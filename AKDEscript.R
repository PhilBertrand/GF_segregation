## Data are available via DRYAD repository: Available upon acceptance of the manuscript by *Scientific Reports*

########################################################################################
############################# ctmm Home range analysis #################################
########################################################################################

## For more in-depth information, please see ctmmweb package usage page https://ctmm-initiative.github.io/ctmmwebdoc/articles/package_usage.html
## And key-papers regarding these analyses:
## Fleming CH, Fagan WF, Mueller T, Olson KA, Leimgruber P, Calabrese JM. 2015. Rigorous home range estimation with movement data: a new autocorrelated kernel density estimator. Ecology. 96(5):1182–1188. doi:10.1890/14-2010.1.
## Calabrese JM, Fleming CH, Gurarie E. 2016. Ctmm: an R package for analyzing animal relocation data as a continuous-time stochastic process. Methods Ecol Evol. 7(9):1124–1132. doi:10.1111/2041-210X.12559.
## Winner K, Noonan MJ, Fleming CH, Olson KA, Mueller T, Sheldon D, Calabrese JM. 2018. Statistical inference for home range overlap. Methods Ecol Evol. 9(7):1679–1691. doi:10.1111/2041-210X.13027.
## Fleming CH, Noonan MJ, Medici EP, Calabrese JM. 2019. Overcoming the challenge of small effective sample sizes in home-range estimation. Methods Ecol Evol. 10(10):1679–1689. doi:10.1111/2041-210X.13270.

require(ctmm)
require(ctmmweb)
require(rgeos)
require(track2kba)
require(move)

tmp <- data ## your data, might be performed on different subgrouping (e.g., population, colony)
ntr <- aggregate(birdTrip ~ ring*colony, tmp, function(x) length(unique(x))) ## Weighing data, here nb of trips per individual 

## Projection shift, ideally metric based
crdref <- CRS("+init=epsg:4326")
pts <- SpatialPointsDataFrame(tmp[,c("Longitude", "Latitude")], tmp,proj4string=crdref)
s <- spTransform(pts, "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
df <- as.data.frame(s)

ex <- extent(120000, 460000, 8600000, 8900000) ## extent specification

## Data formatting
t <- move(x=df$Longitude.2, y=df$Latitude.2, ## Using move package facilitating data format
          time=as.POSIXct(df$datetime,
          format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
          proj=CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"),
          data=df, animal=df$trackID)

data <- as.telemetry(t) ## ctmm format

mtr <- par_try_models(data) ## Parallele model fitting
dfms <- data.frame(summary_tried_models(mtr)) ## Model outputs
df <- subset(dfms, dfms$dAICc == 0) ## Selecting of best model based on AICc
ml <- flatten_models(mtr) ## Conversion
namesSUB <- data.frame(df$model_name, df$identity)
colnames(namesSUB) <- c("model_name", "identity")
mls2 <- ml[namesSUB$model_name]
vls2 <- vario_list[namesSUB$identity]

plot_vario(vls2, mls2, model_color = "purple", fraction=1) ## Model diagnostic; good fit of the function?

tl <- data[namesSUB$identity]
memory.limit(size=800000) ## Updating memory use
hrange <- akde(tl, CTMM = mls2, grid=list(extent = ex, dr = 500)) ## AKDE analysis
names(hrange) <- names(mls2)
la <- lapply(hrange, summary, IC = "LOOCV")

## Population-level home range
mhr <- mean(hrange, weights = ntr$birdTrip)
summary(mhr)
plot(mhr, level.UD = 0.95)

## Sensitivity analysis
tt <- lapply(1:length(hrange), function(x) raster(hrange[[x]], DF = "PDF"))
rt <- stack(tt)
names(rt) <- paste0(df$identity)
repA <- repAssess(s, rt, avgMethod = "mean", nCores = 8, levelUD = 50, bootTable = T, iteration = 999)  

## Visualization          
ggplot() +
  geom_point(repA[[2]], mapping = aes(x=SampleSize, y=InclusionMean), shape = 21, size = 2, 
             fill = "#238A8DFF", alpha = 0.8) +
  geom_path(repA[[2]], mapping = aes(x=SampleSize, y=pred), size = 1, color = "#440154FF", linetype = "dashed", show.legend=FALSE) +
  theme_bw() +
  theme(axis.text=element_text(size=10)) +
  labs(x = "Sample size (n)", y = "Inclusion") +
  scale_fill_viridis_d(begin = 1, end = 0.5) +
  scale_x_continuous(limits = c(0,12.5))
