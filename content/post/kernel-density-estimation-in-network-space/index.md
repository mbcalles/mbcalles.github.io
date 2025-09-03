---
title: Kernel Density Estimation of Point Processes in Network Space in R
summary: Run KDE estimates on a road network.
date: 2021-01-22
authors:
  - admin
tags:
  - R programming
  - Network analysis
  - Kernel density
math: true
---

# Kernel Density Estimation of Point Processes in Network Space in R

- [Running Code](#running-code)
- [Set Up](#set-up)
- [Getting the example data](#getting-the-example-data)
  - [Study Area Extent](#study-area-extent)
  - [Road Network Data](#road-network-data)
  - [Police Report Data](#police-report-data)
- [Prepare the network data](#prepare-the-network-data)
- [Prepare the network data](#prepare-the-network-data-1)
  - [Lixelize the cleaned network](#lixelize-the-cleaned-network)
- [Prepare the point data](#prepare-the-point-data)
- [Calculate Network Kernel Density](#calculate-network-kernel-density)
- [Mapping out the results](#mapping-out-the-results)

**Update (September 2025):** This was a blog post originally from January 2021 where I implemented a network-based KDE in R based on the algorithm outlined in [Xie & Yan (2008)](https://www.sciencedirect.com/science/article/pii/S0198971508000318). At the time there was not a ton of R functionality for network analysis; packages like `sfnetworks` did not yet exist and neither did `spNetwork` to conduct network-based KDE. A much more thorough examination of network-based KDEs and their implementation in R (current as of 2025) can be found here: https://jeremygelb.github.io/spNetwork/articles/NKDE.html. I have updated this post slightly (including code for extracting police data and using `stplanr` for network cleaning). Most of the code still holds up as it relied on igraph and tidygraph directly. While likely not of particular use anymore given the advanced functionality found in `spNetwork`, this little post can stand as a relic of its time.

The network KDE is a 1-D version of the planar kernel density estimator, with $\tau$ (bandwidth) based on network distances rather than Euclidean distances, and the output is based on *lixels* (a 1-D version of pixels) rather than pixels across 2-D Euclidean space.

To produce kernel density estimates (KDE) of point processes in a linear network:

$$
\lambda(z)= \sum_{i=1}^{n} \frac{1}{\tau}\, k\!\left(\frac{d_{iz}}{\tau}\right) y_i
$$

Where,

- $\lambda(z)$ is the density at location $z$;
- $\tau$ is the bandwidth (linear network distance);
- $k$ is the kernel function, typically a function of the ratio of $d_{iz}$ to $\tau$;
- $d_{iz}$ is the linear network distance from event $i$ to location $z$.

I wanted to implement a network-based KDE in R based on the algorithm outlined in [Xie & Yan (2008)](https://www.sciencedirect.com/science/article/pii/S0198971508000318). The network KDE is a 1-D version of the planar kernel density estimator, with $\tau$ (bandwidth) based on network distances rather than Euclidean distances, and the output is based on *lixels* (a 1-D version of pixels) rather than pixels across 2-D Euclidean space.

In this post I use data from Vancouver, BC as a case study for implementing a network kernel density estimator of points in a network in R.

# Set Up

``` r
library(tidygraph)
library(igraph)
library(stplanr)
library(cancensus)
library(osmdata)
library(dplyr)
library(sf)
library(sfnetworks)
library(ggplot2)
library(stringr)
```

# Getting the example data

The first step is to load example data.

## Study Area Extent

I use the `get_census` function from the [`cancensus`](https://mountainmath.github.io/cancensus/index.html) package to extract a `sf` object with `POLYGON` geometries, representing the study extent: the City of Vancouver.

``` r
study_area <- get_census(dataset='CA16', regions=list(CSD="5915"),
                     level='CSD', use_cache = FALSE,geo_format = "sf") %>%
  filter(name=="Vancouver (CY)") %>%
  st_transform(crs = 26910) 

ggplot(study_area) +
  geom_sf() + 
  coord_sf() +
  theme_void()
```

![](kde_files/figure-commonmark/unnamed-chunk-3-1.png)

## Road Network Data

Next I use the `osmdata` package to download street network files for the City of Vancouver. The `getbb()` function defines a bounding box for the City of Vancouver.

``` r
bbx <- getbb("Vancouver, BC")
```

Next we use the `opq()` and `add_osm_feature` functions to obtain open street map road network data. The `osmdata_sf()`

``` r
streets <- bbx %>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value=c("residential", "living_street",
                            "service","unclassified",
                            "pedestrian", "footway",
                            "track","path","motorway", "trunk",
                          "primary","secondary", 
                          "tertiary","motorway_link",
                          "trunk_link","primary_link",
                          "secondary_link",
                          "tertiary_link")) %>%
  osmdata_sf()
  
  
streets_sf <- st_transform(streets$osm_lines,crs=26910) %>% 
  st_intersection(.,st_geometry(study_area)) |> 
  st_cast(to = "MULTILINESTRING") |> 
  st_cast(to = "LINESTRING")

ggplot() +
  geom_sf(data = streets_sf,
          aes(color=highway),
          size = .4,
          alpha = .65)+
  theme_void()
```

![](kde_files/figure-commonmark/unnamed-chunk-5-1.png)

## Police Report Data

Finally I get point data of police reports in Vancouver from the City of Vancouver. I subset the data to only include vehicle collisions reported to police in 2018 (last full year of data available).

``` r
# core base + tidy read
url  <- "https://geodash.vpd.ca/opendata/crimedata_download/crimedata_csv_all_years.zip"
zipf <- tempfile(fileext = ".zip")
csvd <- tempfile(fileext = ".csv")

download.file(url, destfile = zipf, mode = "wb")
csv_name <- unzip(zipf, list = TRUE)$Name[grepl("\\.csv$", unzip(zipf, list=TRUE)$Name)]
unzip(zipf, files = csv_name, exdir = dirname(csvd))
file.rename(file.path(dirname(csvd), csv_name), csvd)
```

    [1] TRUE

``` r
# fast, friendly read
library(readr)
crime <- readr::read_csv(csvd, show_col_types = FALSE)
dplyr::glimpse(crime)
```

    Rows: 925,938
    Columns: 10
    $ TYPE          <chr> "Break and Enter Commercial", "Break and Enter Commercia…
    $ YEAR          <dbl> 2022, 2022, 2023, 2022, 2022, 2022, 2022, 2022, 2022, 20…
    $ MONTH         <dbl> 1, 1, 9, 6, 3, 3, 2, 2, 4, 8, 9, 2, 4, 4, 5, 2, 4, 10, 8…
    $ DAY           <dbl> 5, 3, 14, 17, 15, 19, 23, 25, 30, 2, 11, 24, 1, 28, 11, …
    $ HOUR          <dbl> 7, 16, 3, 5, 5, 6, 23, 10, 22, 15, 0, 4, 4, 12, 18, 22, …
    $ MINUTE        <dbl> 34, 19, 30, 16, 14, 42, 0, 15, 45, 31, 5, 8, 7, 5, 0, 5,…
    $ HUNDRED_BLOCK <chr> "10XX ALBERNI ST", "10XX ALBERNI ST", "10XX ALBERNI ST",…
    $ NEIGHBOURHOOD <chr> "West End", "West End", "West End", "West End", "West En…
    $ X             <dbl> 491015.9, 491036.1, 491065.3, 491067.3, 491102.2, 491102…
    $ Y             <dbl> 5459166, 5459146, 5459130, 5459115, 5459092, 5459092, 54…

``` r
collisions <- crime %>% 
  filter(YEAR==2025 & str_detect(TYPE,"Vehicle Collision"))%>% 
  filter(X!=0) %>%
  st_as_sf(.,coords = c("X","Y"),crs=26910)
```

Locations of police reported vehicle collisions are below.

``` r
ggplot(collisions) +
  geom_sf() + 
  coord_sf() + 
  theme_void()
```

![](kde_files/figure-commonmark/unnamed-chunk-7-1.png)

Now we have the data we need to try and estimate the density of police-reported collisions or “hotspots” in the City of Vancouver.

NOTE: These data do not differentiate between different types of collisions between road users, and police road collision data consistently miss out on the majority of crashes that occur to active transport users (bicyclists and pedestrians).

# Prepare the network data

In the algorithm outlined in [Xie & Yan (2008)](https://www.sciencedirect.com/science/article/pii/S0198971508000318) the resulting density estimates on a network are based on a spatial unit of analysis they term a *lixel*. *Lixels* are basically the 1-D version of a pixel, and refer to how long the segments in our network are from which we calculate densities of events.

To create a *lixelized* network I combine the `split_lines()` and the `lixelize network()` functions.

The `split_lines()` function has three arguments:

- **input_lines**: `sf` object with `LINSETRING` or `MULTILINESTRING` geometry
- **max_length**: the *lixel* size, defined by the maximum length of an individual linestring in the data.
- **id**: a unique id for each line segment in the `sf` object

The `split_lines()` function is given below:

``` r
split_lines <- function(input_lines, max_length, id) {

  input_lines <- input_lines %>% ungroup()

  geom_column <- attr(input_lines, "sf_column")

  input_crs <- sf::st_crs(input_lines)

  input_lines[["geom_len"]] <- sf::st_length(input_lines[[geom_column]])

  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])

  too_short <- filter(select(all_of(input_lines),all_of(id), all_of(geom_column), geom_len), geom_len < max_length) %>% select(-geom_len)

  too_long <- filter(select(all_of(input_lines),all_of(id), all_of(geom_column), geom_len), geom_len >= max_length)

  rm(input_lines) # just to control memory usage in case this is big.

  too_long <- mutate(too_long,
                     pieces = ceiling(geom_len / max_length),
                     fID = 1:nrow(too_long)) %>%
    select(-geom_len)

  split_points <- sf::st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),] %>%
    select(-pieces)

  split_points <- mutate(split_points, split_fID = row.names(split_points)) %>%
    group_by(fID) %>%
    mutate(piece = 1:n()) %>%
    mutate(start = (piece - 1) / n(),
           end = piece / n()) %>%
    ungroup()

  new_line <- function(i, f, t) {
    lwgeom::st_linesubstring(x = too_long[[geom_column]][i], from = f, to = t)[[1]]
  }

  split_lines <- apply(split_points[c("fID", "start", "end")], 1,
                       function(x) new_line(i = x[["fID"]], f = x[["start"]], t = x[["end"]]))

  rm(too_long)

  split_lines <- st_sf(split_points[c(id)], geometry = st_sfc(split_lines, crs = input_crs))

  lixel <- rbind(split_lines,too_short) %>% mutate(LIXID = row_number())

  return(lixel)
}
```

The second function needed to prepare the network data is the `lixelize_network` function. This takes two arguments, an `sf` object with `LINESTRING` or `MULTILINESTRING` geometry and the lixel_length argument. This function will create two lixelized networks: 1) the output network where the line data are *lixelized* and no segment is larger than the a given pixel size and 2) a shortest distance network which will be used to calculate road network lengths in our network kernel density estimates that we will dive into shortly.

The lixelize network function is given below:

``` r
lixelize_network <- function(sf_network,max_lixel_length,uid){
  print("Splitting input spatial lines by lixel length...")
  target_lixel <- split_lines(input_lines = sf_network,max_length = max_lixel_length,id = uid)
  print("Create corresponding shortest distance network...")
  shortest_distance_network <- split_lines(input_lines = target_lixel,max_length = max_lixel_length/2,id = uid)
  return(list(target_lixel=target_lixel,shortest_distance_network=shortest_distance_network))
}
```

Our approach varies from Xie and Yan (2008) in that we define a “maximum” length for lixels based on the `split_lines()` function, and segments that are larger than the specified length are split into equal sized lixels. Instead of having $n$ lixels with 1 extra “residual lixel” we instead have $n$+1 lixels of equal length, under the lixel limit. For example a 45m segment with a maximum lixel length of 10m would be split into 5 segments of 9m instead of 4 lixels of 10m and a residual lixel of 5m.

Now we have the data we need to try and estimate the density of police-reported collisions or “hotspots” in the City of Vancouver.

NOTE: These data do not differentiate between different types of collisions between road users, and police road collision data consistently miss out on the majority of crashes that occur to active transport users (bicyclists and pedestrians).

# Prepare the network data

In the algorithm outlined in [Xie & Yan (2008)](https://www.sciencedirect.com/science/article/pii/S0198971508000318) the resulting density estimates on a network are based on a spatial unit of analysis they term a *lixel*. *Lixels* are basically the 1-D version of a pixel, and refer to how long the segments in our network are from which we calculate densities of events.

To create a *lixelized* network I combine the `split_lines()` and the `lixelize network()` functions.

The `split_lines()` function has three arguments:

- **input_lines**: `sf` object with `LINSETRING` or `MULTILINESTRING` geometry
- **max_length**: the *lixel* size, defined by the maximum length of an individual linestring in the data.
- **id**: a unique id for each line segment in the `sf` object

The `split_lines()` function is given below:

``` r
split_lines <- function(input_lines, max_length, id) {

  input_lines <- input_lines %>% ungroup()

  geom_column <- attr(input_lines, "sf_column")

  input_crs <- sf::st_crs(input_lines)

  input_lines[["geom_len"]] <- sf::st_length(input_lines[[geom_column]])

  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])

  too_short <- filter(select(input_lines,all_of(id), all_of(geom_column), geom_len), geom_len < max_length) %>% select(-geom_len)

  too_long <- filter(select(input_lines,all_of(id), all_of(geom_column), geom_len), geom_len >= max_length)

  rm(input_lines) # just to control memory usage in case this is big.

  too_long <- mutate(too_long,
                     pieces = ceiling(geom_len / max_length),
                     fID = 1:nrow(too_long)) %>%
    select(-geom_len)

  split_points <- sf::st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),] %>%
    select(-pieces)

  split_points <- mutate(split_points, split_fID = row.names(split_points)) %>%
    group_by(fID) %>%
    mutate(piece = 1:n()) %>%
    mutate(start = (piece - 1) / n(),
           end = piece / n()) %>%
    ungroup()

  new_line <- function(i, f, t) {
    lwgeom::st_linesubstring(x = too_long[[geom_column]][i], from = f, to = t)[[1]]
  }

  split_lines <- apply(split_points[c("fID", "start", "end")], 1,
                       function(x) new_line(i = x[["fID"]], f = x[["start"]], t = x[["end"]]))

  rm(too_long)

  split_lines <- st_sf(split_points[c(id)], geometry = st_sfc(split_lines, crs = input_crs))

  lixel <- rbind(split_lines,too_short) %>% mutate(LIXID = row_number())

  return(lixel)
}
```

The second function needed to prepare the network data is the `lixelize_network` function. This takes two arguments, an `sf` object with `LINESTRING` or `MULTILINESTRING` geometry and the lixel_length argument. This function will create two lixelized networks: 1) the output network where the line data are *lixelized* and no segment is larger than the a given pixel size and 2) a shortest distance network which will be used to calculate road network lengths in our network kernel density estimates that we will dive into shortly.

The lixelize network function is given below:

``` r
lixelize_network <- function(sf_network,max_lixel_length,uid){
  print("Splitting input spatial lines by lixel length...")
  target_lixel <- split_lines(input_lines = sf_network,max_length = max_lixel_length,id = uid)
  print("Create corresponding shortest distance network...")
  shortest_distance_network <- split_lines(input_lines = target_lixel,max_length = max_lixel_length/2,id = uid)
  return(list(target_lixel=target_lixel,shortest_distance_network=shortest_distance_network))
}
```

Our approach varies from Xie and Yan (2008) in that we define a “maximum” length for lixels based on the `split_lines()` function, and segments that are larger than the specified length are split into equal sized lixels. Instead of having $n$ lixels with 1 extra “residual lixel” we instead have $n$+1 lixels of equal length, under the lixel limit. For example a 45m segment with a maximum lixel length of 10m would be split into 5 segments of 9m instead of 4 lixels of 10m and a residual lixel of 5m.

## Lixelize the cleaned network

Next step is to lixelize the cleaned network. I choose a 25m maximum lixel length (e.g. every linear segment in the network will be 25m or under. This function may take some time to run.

``` r
lixel_list_25m <- lixelize_network(
  sf_network = streets_sf,
  max_lixel_length = 25,
  uid = "osm_id"
    
)
```

    [1] "Splitting input spatial lines by lixel length..."
    [1] "Create corresponding shortest distance network..."

``` r
str(lixel_list_25m)
```

    List of 2
     $ target_lixel             : sf [224,904 × 3] (S3: sf/tbl_df/tbl/data.frame)
      ..$ osm_id  : chr [1:224904] "4231647" "4231647" "4231647" "4231647" ...
      ..$ geometry:sfc_LINESTRING of length 224904; first list element:  'XY' num [1:2, 1:2] 487849 487873 5455038 5455039
      ..$ LIXID   : int [1:224904] 1 2 3 4 5 6 7 8 9 10 ...
      ..- attr(*, "sf_column")= chr "geometry"
      ..- attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA
      .. ..- attr(*, "names")= chr [1:2] "osm_id" "LIXID"
     $ shortest_distance_network: sf [436,612 × 3] (S3: sf/tbl_df/tbl/data.frame)
      ..$ osm_id  : chr [1:436612] "4231647" "4231647" "4231647" "4231647" ...
      ..$ geometry:sfc_LINESTRING of length 436612; first list element:  'XY' num [1:2, 1:2] 487849 487861 5455038 5455039
      ..$ LIXID   : int [1:436612] 1 2 3 4 5 6 7 8 9 10 ...
      ..- attr(*, "sf_column")= chr "geometry"
      ..- attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA
      .. ..- attr(*, "names")= chr [1:2] "osm_id" "LIXID"

# Prepare the point data

Here I need to ensure that the point data are lined up with the cleaned network data. I use a function called `st_snap_points` which “snaps” the points to the location of the nearest road segment within a specified distance tolerance. The function is given below:

``` r
st_snap_points = function(x, y, max_dist = 1000) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}
```

We snap the points to our road network with a 20m tolerance:

``` r
coll_snp <- st_snap_points(collisions,streets_sf,max_dist = 20) #this only returns the geometry - doesn't preserve attributes

coll_snp <- st_sf(collisions %>%
                  st_drop_geometry() %>%
                  mutate(geom=coll_snp)) #rerturn the attributes to the snapped locations
```

# Calculate Network Kernel Density

The function I have created here for calculating density of events on a network follows the algorithm outline in Xie and Yan (2008) and takes advantage of functionality from `sf`, `igraph` and `tidygraph` to to calculate shortest distances from each lixel in our network to the point processes we specify. Much of the code in this function is adapted from an excellent resource on using these packages for network analysis: https://www.r-spatial.org/r/2019/09/26/spatial-networks.html

The `network_kde()` function is defined below:

``` r
#function for calculating the centre of 
st_line_midpoints <- function(sf_lines = NULL) {
  
  g <- st_geometry(sf_lines)
  
  g_mids <- lapply(g, function(x) {
    
    coords <- as.matrix(x)
    
    # this is just a copypaste of maptools:::getMidpoints()):
    get_mids <- function (coords) {
      dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid <- sum(dist)/2
      dist_cum <- c(0, cumsum(dist))
      end_index <- which(dist_cum > dist_mid)[1]
      start_index <- end_index - 1
      start <- coords[start_index, ]
      end <- coords[end_index, ]
      dist_remaining <- dist_mid - dist_cum[start_index]
      mid <- start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    
    mids <- st_point(get_mids(coords))
  })
  
  out <- st_sfc(g_mids, crs = st_crs(sf_lines))
  out <- st_sf(out)
  out <- bind_cols(out,sf_lines %>% st_drop_geometry()) %>% rename(geom=out)
}


network_kde <- function(lixel_list,point_process,bandwidth = 100,n_cores=1,attribute=1,point_process_is_lixel_midpoint=FALSE){
  

  #3. Create a network of lixels by establishing the network topology between lixels as well as between lixels and lxnodes.

  #Define topology for calculating shortest distances by converting lixel with half the length to calculate distances
  #The shorter the lixel length the more accurate the calculation of network distance from nearest node of source lixel to nearest node of target lixel
  print("Defining network topology...")
  
  require(parallel)
  require(doParallel)
  require(igraph)
  require(tidygraph)
  require(sf)  
  require(dplyr)


  # Create nodes at the start and end point of each edge

  nodes <- lixel_list$shortest_distance_network %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(LIXID = L1) %>%
    group_by(LIXID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2))

  # Give each node a unique index

  nodes <- nodes %>%
    mutate(xy = paste(.$X, .$Y),
           xy = factor(xy, levels = unique(xy))) %>%
    group_by(xy)%>%
    mutate(nodeID = cur_group_id()) %>%
    ungroup() %>%
    select(-xy)

  # Combine the node indices with the edges

  start_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)

  end_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)

  lixel_list$shortest_distance_network = lixel_list$shortest_distance_network %>%
    mutate(from = start_nodes, to = end_nodes)

  # Remove duplicate nodes
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(LIXID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(lixel_list$shortest_distance_network))

  # Convert to tbl_graph
  graph <- tbl_graph(nodes = nodes, edges = as_tibble(lixel_list$shortest_distance_network), directed = FALSE)

  graph <- graph %>%
    activate(edges) %>%
    mutate(length = st_length(geometry))

  lixel_list$shortest_distance_network <- NULL

  #4. Create the center points of all the lixels for the target lixel (lxcenters)
  print("Calculating lixel midpoints (lxcenters)...")

  lxcenters <- st_line_midpoints(lixel_list$target_lixel)

  #5. Select a point process occurring within the road network

  #input as function parameter

  #6. For each point find its nearest lixel. Count the number of points nearest to a lixel and assigned as property of lixel.

  #Points should be snapped to network within some distance threshold prior to

  if(point_process_is_lixel_midpoint==FALSE){

  print("Counting number of events on each lixel...")

  point_process <- st_join(point_process,lixel_list$target_lixel["LIXID"],
                           join = st_is_within_distance, dist = 0.001) #for each point assign the LIXID that it falls on

  source_lixels <- point_process %>% #summarize the number of points by LIXID (e.g. count the points for each LIXID)
    filter(!is.na(LIXID)) %>%
    group_by(LIXID) %>%
    summarise(n_events = n(),`.groups`="drop") %>%
    st_drop_geometry()

  source_lixels <- inner_join(lxcenters,source_lixels,by="LIXID") #define geometry for source lixels
  }

  if(point_process_is_lixel_midpoint==TRUE){

    source_lixels <- lxcenters %>% mutate(n_events = 1)
  }

  print(paste0(sum(source_lixels$n_events)," events on ",nrow(source_lixels)," source lixels"))

  #7. Define a search bandwidth, measured with the short path network distance

  #input as function parameter

  #8. Calculate the shortest-path network distance from each source lixel to lxcenter of all its neighouring lixels within the seach bandwidth


  nearest_node_to_source_lixel <- st_nearest_feature(source_lixels,nodes)#find nodes from shortest distance network associated with each source lixel
  nearest_node_to_lxcentre <- st_nearest_feature(lxcenters,nodes)#find nodes from shortest distance network associated with each lxcenters

  print("Calculating distances from each source lixel to all other lixel centers... ")
  

  cl <- makeCluster(n_cores)
  registerDoParallel(cores=cl) #parallel computing

  distances <- foreach::foreach(i = 1:length(nearest_node_to_source_lixel),.packages = c("magrittr","igraph","tidygraph","sf")) %dopar% {
    temp <- distances(
      graph = graph,
      weights = graph %>% activate(edges) %>% pull(length),
      v = nearest_node_to_source_lixel[i]
    )[,nearest_node_to_lxcentre]

    data.frame(LIXID = lxcenters[temp<=max(bandwidth),]$LIXID,
               dist = temp[temp<=max(bandwidth)])
  }

 stopCluster(cl)

  rm("graph")


  # gauss <-function(x) {
  #
  #   t1 <- 1/(sqrt(2*pi))
  #   t2 <- exp(-((x^2)/(2*bandwidth^2)))
  #   n <- (1/bandwidth)*(t1*t2)
  #   n <- ifelse(x>bandwidth,0,n)
  #
  #
  #   return(n)
  #
  # }
  #
  # plot(gauss(0:(bandwidth*2)))

  quartic <-function(x,r) {
    K <- 3/pi
    t1 <- 1-(x^2/r^2)
    q <- ifelse(x>r,0,K*t1)
    q <- (1/r)*q

    return(q)
  }

  LIXID <- unlist(lapply(lapply(distances,`[`,"LIXID"),function(x) pull(x)))
  distances_list <- lapply(lapply(distances,`[`,"dist"),function(x) pull(x))

  d_list <- list()

  for(i in 1:length(bandwidth)){

    print(paste0("Applying kernel density estimator with bandwidth of ",bandwidth[i],"m ... "))

    d <- lapply(distances_list,function(x) quartic(x,bandwidth[i]))
    d <- mapply(`*`,d,attribute) #multiply by attribute of source lixel
    d <- mapply(`*`,d,source_lixels$n_events) #sum over number of events that occured on the source lixel

    d_list[[i]] <- data.frame(unlist(d))
  }

  d_cols <- as.data.frame(do.call(cbind,d_list))

  names(d_cols) <- paste0("kde_bw_",bandwidth)

  #sum densities over all lixels
  density <- cbind(LIXID,d_cols) %>%
    group_by(LIXID) %>%
    summarise(across(everything(),sum),`.groups`="drop")

  network_kde <- left_join(lixel_list$target_lixel,density,by = "LIXID") %>%
    mutate(length = st_length(.)) %>%
    replace(., is.na(.), 0)

  return(list(network_kde=network_kde,neighbours = distances))

}
```

As of now I have only set up the algorithm to use a quartic kernel.

Next step is to compute the density of points on a network! The `network_kde()` function implements parallel computing functionality using the `parallel` and `doParallel` packages so we need to specify how many cores are on the computer being used:

``` r
n_cores <- bigstatsr::nb_cores()
```

Finally we specify the lixelized network we are working with, the point process from which we estimate densities of events, the bandwidths we wish to use in our density estimates and the number of cores for parallel computation. The last argument can be ignored for now and set to FALSE as it is something I’m working on to speed up computation if your point process is already associated with the segment (e.g. volumes of road users on a segment).

``` r
system.time({
  collisions_network_kde <- network_kde(lixel_list = lixel_list_25m,
                                 point_process = coll_snp,
                                 bandwidth = c(50,100,150,250),
                                 n_cores = n_cores,
                                 point_process_is_lixel_midpoint = FALSE)
})
```

    [1] "Defining network topology..."

    [1] "Calculating lixel midpoints (lxcenters)..."
    [1] "Counting number of events on each lixel..."
    [1] "733 events on 623 source lixels"
    [1] "Calculating distances from each source lixel to all other lixel centers... "
    [1] "Applying kernel density estimator with bandwidth of 50m ... "
    [1] "Applying kernel density estimator with bandwidth of 100m ... "
    [1] "Applying kernel density estimator with bandwidth of 150m ... "
    [1] "Applying kernel density estimator with bandwidth of 250m ... "

       user  system elapsed 
     105.76    2.36  182.16 

``` r
glimpse(collisions_network_kde$network_kde)
```

    Rows: 224,904
    Columns: 8
    $ osm_id     <chr> "4231647", "4231647", "4231647", "4231647", "4231647", "423…
    $ geometry   <LINESTRING [m]> LINESTRING (487848.5 545503..., LINESTRING (4878…
    $ LIXID      <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, …
    $ kde_bw_50  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    $ kde_bw_100 <dbl> 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.00000…
    $ kde_bw_150 <dbl> 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.00000…
    $ kde_bw_250 <dbl> 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.00000…
    $ length     [m] 24.40429 [m], 24.40429 [m], 24.40429 [m], 24.40429 [m], 24.40…

The output of the `network_kde()` is an `sf` object with `LINESTRING` geometry, where each line segment is a lixel and the columns refer to density of events in network space for the specified bandwidths.

# Mapping out the results

Below we can then map out the results. We will use the 150m bandwidth in this example:

``` r
ggplot() +
  geom_sf(data=study_area,color=NA,fill="gray5")+
  geom_sf(data = collisions_network_kde$network_kde,aes(color=kde_bw_150))+
  scale_color_viridis_c(option = "C",direction = 1,name="Density") + 
  # geom_sf(data=mask,color=NA,fill="white")+
  theme_void()+
  theme(panel.background= element_rect(fill = NA,color=NA),
        legend.position = "top",
        legend.title.align = 0,
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'))+
  coord_sf(xlim = c(483690.9,499329.2),
           ylim = c(5450534.8,5463381.0))
```

![](kde_files/figure-commonmark/unnamed-chunk-18-1.png)

Zooming into downtown:

``` r
ggplot() +
  geom_sf(data=study_area,color=NA,fill="gray5")+
  geom_sf(data = collisions_network_kde$network_kde,aes(color=kde_bw_150))+
  scale_color_viridis_c(option = "C",direction = 1,name="Density") + 
  theme_void()+
  theme(panel.background= element_rect(fill = NA,color=NA),
        legend.position = "top",
        legend.title.align = 0,
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'))+
  coord_sf(xlim = c(491321.39853389-1000,491321.39853389+1100),
           ylim = c(5459015.689874-1000,5459015.689874+1000))
```

![](kde_files/figure-commonmark/unnamed-chunk-19-1.png)

That’s the network KDE function as it is for now. Some challenges going forward are computation time for higher resolution pixels and implementing different kernels other than quartic.

<details>
<summary>
Reproducibility receipt
</summary>

    [1] "2025-09-03 16:15:49 PDT"

    Local:    main C:/Users/micha/Documents/GitHub/mbcalles.github.io
    Remote:   main @ origin (https://github.com/mbcalles/mbcalles.github.io)
    Head:     [2a60de6] 2025-09-03: publication

    R version 4.4.1 (2024-06-14 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26100)

    Matrix products: default


    locale:
    [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8   
    [3] LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=English_Canada.utf8    

    time zone: America/Vancouver
    tzcode source: internal

    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2     readr_2.1.5      
     [5] stringr_1.5.1     ggplot2_3.5.1     sfnetworks_0.6.4  sf_1.0-18        
     [9] dplyr_1.1.4       osmdata_0.2.5     cancensus_0.5.7   stplanr_1.2.2    
    [13] igraph_2.1.1      tidygraph_1.3.1  

    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1       viridisLite_0.4.2      farver_2.1.2          
     [4] fastmap_1.2.0          spatstat.geom_3.3-3    spatstat.explore_3.3-3
     [7] bigassertr_0.1.7       flock_0.7              digest_0.6.37         
    [10] rpart_4.1.23           timechange_0.3.0       lifecycle_1.0.4       
    [13] spatstat.data_3.1-2    magrittr_2.0.3         compiler_4.4.1        
    [16] rlang_1.1.4            tools_4.4.1            utf8_1.2.4            
    [19] yaml_2.3.10            knitr_1.48             labeling_0.4.3        
    [22] bit_4.5.0              classInt_0.4-10        curl_5.2.3            
    [25] xml2_1.3.6             rmio_0.4.0             abind_1.4-8           
    [28] KernSmooth_2.23-24     bigstatsr_1.6.2        withr_3.0.2           
    [31] purrr_1.0.2            grid_4.4.1             polyclip_1.10-7       
    [34] fansi_1.0.6            git2r_0.35.0           e1071_1.7-16          
    [37] colorspace_2.1-1       scales_1.3.0           spatstat.utils_3.1-0  
    [40] cli_3.6.3              rmarkdown_2.28         crayon_1.5.3          
    [43] generics_0.1.3         rstudioapi_0.17.1      httr_1.4.7            
    [46] tzdb_0.4.0             DBI_1.2.3              proxy_0.4-27          
    [49] splines_4.4.1          spatstat.model_3.3-2   vctrs_0.6.5           
    [52] Matrix_1.7-0           jsonlite_1.8.9         geojsonsf_2.0.3       
    [55] hms_1.1.3              bit64_4.5.2            tensor_1.5            
    [58] bigparallelr_0.3.2     spatstat.univar_3.0-1  tidyr_1.3.1           
    [61] units_0.8-5            goftest_1.2-3          parallelly_1.45.1     
    [64] glue_1.8.0             spatstat.random_3.3-2  lwgeom_0.2-14         
    [67] codetools_0.2-20       cowplot_1.1.3          lubridate_1.9.3       
    [70] stringi_1.8.4          gtable_0.3.6           deldir_2.0-4          
    [73] munsell_0.5.1          tibble_3.2.1           pillar_1.9.0          
    [76] rappdirs_0.3.3         htmltools_0.5.8.1      R6_2.5.1              
    [79] httr2_1.0.5            sfheaders_0.4.4        vroom_1.6.5           
    [82] evaluate_1.0.1         lattice_0.22-6         class_7.3-22          
    [85] Rcpp_1.0.13            spatstat.linnet_3.2-2  nlme_3.1-164          
    [88] spatstat.sparse_3.1-0  mgcv_1.9-1             xfun_0.48             
    [91] pkgconfig_2.0.3       

</details>
