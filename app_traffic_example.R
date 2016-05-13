# @knitr data_prep_example
library(dplyr)
library(geosphere)

load("data.RData") # a data frame, d, containing 'long' and 'lat' columns
p <- SpatialPoints(cbind(d$long, d$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
idx1 <- 69 # great circles from coords in all other rows to coords in this row
idx2 <- 2648 # as above

get_paths <- function(x, idx, ...){
  gcInt <- function(x, x1, x2){
    x <- gcIntermediate(x[x1,], x[x2,], ...)
    if(is.list(x)){
      x <- x %>% purrr::map2(c(x1, x1 + 0.5), ~data.frame(.x, .y)) %>%
        bind_rows %>% setnames(c("long", "lat", "group"))
    } else x <- data.frame(x, x1) %>% setnames(c("long", "lat", "group"))
    x
  }
  purrr::map(setdiff(1:length(x), idx), ~gcInt(x, .x, idx)) %>% bind_rows
}

paths1 <- get_paths(p, idx1, addStartEnd=TRUE)
paths2 <- get_paths(p, idx2, addStartEnd=TRUE)

# @knitr setup
library(parallel)
library(gridExtra)
library(raster)
library(data.table)
library(dplyr)
library(ggplot2)

eb <- element_blank()
theme_blank <- theme(axis.line=eb, axis.text.x=eb, axis.text.y=eb,
  axis.ticks=eb, axis.title.x=eb, axis.title.y=eb, legend.position="none",
  panel.background=eb, panel.border=eb, panel.grid.major=eb, panel.grid.minor=eb,
  plot.background=element_rect(colour="transparent", fill="transparent"))

world <- map_data("world")

# @knitr gc_segments
df_segs <- function(d, seg.size, n.frames, replicates=1, direction="fixed"){
  n <- nrow(d)
  if(n < 3) stop("Data not appropriate for this operation.")
  if(seg.size < 3) stop("Segment size too small.")
  z <- round(runif(2, 2, seg.size))
  z[z > n] <- n
  n1 <- ceiling(diff(c((z[1] - z[2]), n))/z[1])
  if(n.frames - n1 < 100) stop("Insufficient frames")
  offset <- sample(0:(n.frames - n1), replicates)
  
  f <- function(k, d, n, n1, z, offset){
    ind2 <- z[1]*k
    ind1 <- max(ind2 - z[2], 1)
    if(ind2 > n) ind2 <- n
    d <- slice(d, ind1:ind2)
    purrr::map(offset, ~mutate(d,
      group=ifelse(replicates==1, group, group + as.numeric(sprintf(".%d", k))),
      frameID=.x + k)) %>% bind_rows
  }
  
  if(direction=="reverse") d <- mutate(d, long=rev(long), lat=rev(lat))
  if(direction=="random" && rnorm(1) < 0) d <- mutate(d, long=rev(long), lat=rev(lat))
  d <- purrr::map(1:n1, ~f(.x, d, n, n1, z, offset)) %>% bind_rows %>% arrange(group, frameID)
  d
}

n.frames <- 900
set.seed(1)
paths2 <- mutate(paths2, group=group + max(paths1$group))
paths <- bind_rows(paths1, paths2) %>% split(.$group)

# @knitr make_segs
paths <- mclapply(paths, df_segs, seg.size=5, n.frames=n.frames, replicates=1, direction="random", mc.cores=32) %>% bind_rows

# @knitr make_segs_alt
paths <- paths %>% split(.$group) %>% purrr::map(~df_segs(.x, 5, n.frames, replicates=1, direction="random")) %>% bind_rows

# @knitr marmap
d.bath <- read.csv("marmap_coord_-180;-90;180;90_res_10.csv") %>%
  data.table %>% setnames(c("long", "lat", "z"))
r <- raster(extent(-180,180,-90,90), res=1/6)
projection(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
r <- setValues(r, d.bath$z)

# @knitr hemisphere
project_to_hemisphere <- function(lat, long, lat0, long0){
    hold <- c(lat, long)
    x <- (pi/180)*c(lat, lat0, long-long0)
    inview <- sin(x[1])*sin(x[2]) + cos(x[1])*cos(x[2])*cos(x[3]) > 0
    data.table(long=hold[2], lat=hold[1], inview=inview)
}

# @knitr plot_setup1
n.period <- 120
lon_seq <- rep(seq(0, 360, length.out=n.period + 1)[-(n.period + 1)], length=n.frames)
lat_seq <- rep(41, length(lon_seq))
paths <- paths %>% split(.$frameID)
d.bath.agg <- r %>% aggregate(2) %>% rasterToPoints %>% data.table %>% setnames(c("long", "lat", "z"))

# @knitr plot_setup2
d.tiles <- mclapply(1:n.period,
  function(i, x, lon, lat){
    left_join(x, project_to_hemisphere(x$lat, x$long, lat[i], lon[i])) %>%
      filter(inview) %>% dplyr::select(-inview) %>% mutate(frameID=i)
    },
  x=d.bath.agg, lat=lat_seq, lon=lon_seq, mc.cores=32)
  
z.range <- purrr::map(d.tiles, ~range(.x$z, na.rm=TRUE)) %>% unlist %>% range
d.world <- purrr::map(1:n.period, ~mutate(world, frameID=.x))

# @knitr plot_func
save_maps <- function(x, lon_seq, lat_seq, col=NULL, type="network", z.range=NULL){
  if(is.null(col)) col <- switch(type,
    network=c("#FFFFFF25", "#1E90FF25", "#FFFFFF", "#1E90FF50"),
    maptiles=c("black", "steelblue4"),
    maplines="white")
  i <- x$frameID[1]
  if(type=="network") x.lead <- group_by(x, group) %>% slice(n())
  g <- ggplot(x, aes(long, lat))
  if(type=="maptiles"){
    if(is.null(z.range)) z.range <- range(x$z, na.rm=TRUE)
    g <- ggplot(x, aes(long, lat, fill=z)) + geom_tile() +
      scale_fill_gradientn(colors=col, limits=z.range)
  } else {
    g <- ggplot(x, aes(long, lat, group=group))
    if(type=="maplines") g <- g + geom_path(colour=col)
    if(type=="network") g <- g + geom_path(colour=col[2]) + geom_path(colour=col[1]) +
      geom_point(data=x.lead, colour=col[3], size=0.6) +
      geom_point(data=x.lead, colour=col[4], size=0.3)
  }
  g <- g + theme_blank + coord_map("ortho", orientation=c(lat_seq[i], lon_seq[i], 23.4))
  dir.create(outDir <- file.path("frames", type), recursive=TRUE, showWarnings=FALSE)
  png(sprintf(paste0(outDir, "/", type, "_%03d.png"), i),
    width=4*1920, height=4*1080, res=300, bg="transparent")
  print(g)
  dev.off()
  NULL
}

mclapply(paths, save_maps, lon_seq, lat_seq, type="network", mc.cores=30)
mclapply(d.world, save_maps, lon_seq, lat_seq, type="maplines", mc.cores=30)
mclapply(d.tiles, save_maps, lon_seq, lat_seq, type="maptiles", z.range=z.range, mc.cores=30)
