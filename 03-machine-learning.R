#load the relevant packages
library(powdR)
library(leaflet)
library(ggplot2)
library(Cubist)
library(plotly)



#Load the soil property data
props <- read.csv(file = "data/clusters_and_properties.csv")

#Load the XRPD data

#Get the full file paths
xrpd_paths <-  dir("data/xrd", full.names = TRUE)

#Load the data
xrpd <- read_xy(files = xrpd_paths)

#Make sure the data are interpolated onto the same 2theta scale
#as there are small differences within the dataset
xrpd <- interpolate(xrpd, new_tth = xrpd$icr030336$tth)


#Check that the names of the xrpd data match the SSN column in props
identical(names(xrpd), props$SSN)


leaflet(props) %>%  #create a leaflet object from props
  addTiles() %>% #add the default tiles for the map surface
  addCircleMarkers(~Longitude, ~Latitude)


leaflet(props[props$Sentinel_site == "Bana",]) %>%
  addTiles() %>%
  addCircleMarkers(~Longitude, ~Latitude)


ggplot(data = props, aes(K)) +
  geom_histogram()

summary(props$K)


plot(as_multi_xy(xrpd[props$SSN[props$Sentinel_site == "Didy"]]),
     wavelength = "Cu",
     normalise = TRUE,
     interactive = TRUE)


plot(as_multi_xy(xrpd[props$Sentinel_site == "Didy"]),
     wavelength = "Cu",
     normalise = TRUE,
     xlim = c(26,27))



#-----------------------------------------------------------------------------
#Data pre-treatment

#Extract a quartz pattern from rockjock
quartz <- as_xy(data.frame(rockjock$tth,
                           rockjock$xrd$QUARTZ))

plot(quartz, wavelength = "Cu")


## ---- fig.cap = "Diffractograms from the Didy site aligned to the quartz pattern extracted from the rockjock reference library.", out.width='100%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE----
#Align the xrpd data to this quartz pattern
#using a restricted 2theta range of 10 to 60
xrpd_pt <- align_xy(xrpd, std = quartz,
                    xmin = 10, xmax = 60,
                    xshift = 0.2)

plot(as_multi_xy(xrpd_pt[props$Sentinel_site == "Didy"]),
     wavelength = "Cu",
     normalise = TRUE,
     xlim = c(26,27))


## ---- fig.cap = "Aligned diffractograms subset to remove data below 5 degrees 2theta.", out.width='100%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE----
subset_xrpd <- function (x, xmin, xmax) {

  x <- x[x[[1]] >= xmin & x[[1]] <= xmax, ]

  return(x)

}

xrpd_pt <- lapply(xrpd_pt, subset_xrpd,
               xmin = 5, xmax = 75)

plot(as_multi_xy(xrpd_pt[props$Sentinel_site == "Didy"]),
     wavelength = "Cu")


## ---- fig.cap = "Aligned, subset and mean centred diffractograms.", out.width='100%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE----
#Create a mean centering function
mean_center <- function(x) {

  x[[2]] <- x[[2]] - mean(x[[2]])

  return(x)

}

#apply the function to all patterns
xrpd_pt <- lapply(xrpd_pt, mean_center)

#Inspect the data from Didy
plot(as_multi_xy(xrpd_pt[props$Sentinel_site == "Didy"]),
     wavelength = "Cu")


## ---- message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------------------------------
#Set the seed for random number generation
#so that results are reproducible.
set.seed(10)

#Randomly select 75% of the samples
selection <- sample(1:nrow(props),
                    size = round(nrow(props)*0.75),
                    replace = FALSE)


## ---- message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------------------------------
#double check that the order of xrpd and props
#data match
identical(names(xrpd_pt), props$SSN)

#create a data frame
xrpd_df <- multi_xy_to_df(as_multi_xy(xrpd_pt),
                          tth = TRUE)

#transpose the data so that each sample is a row
cubist_xrpd <- data.frame(t(xrpd_df[-1]))


## ---- message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------------------------------
library(Cubist)

#Create a Cubist model for K
cubist_K <- cubist(x = cubist_xrpd[selection,],
                   y = props$K[selection],
                   committees = 10)


## ---- out.width='80%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE--------------------------
#Predict K
predict_K <- predict(cubist_K,
                     cubist_xrpd[-selection,],
                     neighbours = 9)

#Plot measured K vs predicted K
plot(x = log(props$K[-selection]), log(predict_K),
     xlab = "log(Measured K), ppm", ylab = "log(Predicted K), ppm")
abline(0,1) #Add a 1:1 line

#R2 of measured K vs predicted K
cor(log(props$K[-selection]), log(predict_K))^2


## ----feature-plot, fig.cap = "Features selected by Cubist for the prediction of K. Grey sticks denote the fraction of variable use in the regression models and black sticks denote the fraction of variable use in decision. Red line is the mean diffractogram for the dataset.", out.width='100%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE----
#Extract the usage data from the Cubist model
K_usage <- cubist_K$usage

#Make sure the variable names are numeric so that they can be ordered
K_usage$Variable <- as.numeric(substr(K_usage$Variable,
                                      2,
                                      nchar(K_usage$Variable)))

#Order the K_usage data by Variable now it is numeric
K_usage <- K_usage[order(K_usage$Variable), ]

#Add the 2theta axis
K_usage$tth <- xrpd_df$tth

#Add a mean diffracotgram to the data
K_usage$counts <- rowMeans(xrpd_df[-1])

#Create a function that normalises a vector data to a minimum of 0
#and maximum of 1.
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#Create the plot
ggplot(data = K_usage) +
geom_linerange(aes(x = tth, ymax = Model/100, ymin = 0),
               colour = "grey81", size = 1) +
geom_linerange(aes(x = tth, ymax = Conditions/100, ymin = 0),
               colour = "grey19", size = 1) +
geom_line(aes(x = tth, y = range01(counts)), colour = "red",
          size = 0.5) +
ylab("Scaled counts and fraction of variable use\n") +
xlab("2theta") +
theme_bw()


## ---- warning=FALSE, cache = TRUE-------------------------------------------------------------------------------------------
f1 <- fps(lib = rockjock,
          smpl = xrpd$icr014764,
          std = "QUARTZ",
          refs = c("K-feldspar",
                   "Quartz",
                   "Mica (Tri)",
                   "Organic matter",
                   "Halloysite",
                   "Kaolinite",
                   "Rutile",
                   "Background"),
          align = 0.2)


## ----advanced-feature-plot, fig.cap = "Combining the results from full pattern summation with the features selected by Cubist for prediction of K concentrations.", out.width='100%', fig.asp=.75, fig.align='center', message=FALSE, warning=FALSE, cache = TRUE----
#Create a plot of the full pattern summation results
p1 <- plot(f1, wavelength = "Cu", group = TRUE)

#Define the layers for the sticks that will have to be placed
#BENEATH the data in p1
p2 <- geom_linerange(data = K_usage,
                     aes(x = tth,
                         ymax = (Model/100)*max(f1$measured),
                         ymin = 0),
                     colour = "grey81", size = 1)

p3 <- geom_linerange(data = K_usage,
                     aes(x = tth,
                         ymax = (Conditions/100)*max(f1$measured),
                         ymin = 0),
                     colour = "grey19", size = 1)

#Order the layers so that p2 and p3 are beneath p1
p1$layers <- c(p2, p3, p1$layers)

#Limit the x-axis to between 25 and 30 degrees so that the
#dominant features can be easily examned
p1 <- p1 +
      scale_x_continuous(limits = c(25, 30))

p1


## ----eval = FALSE-----------------------------------------------------------------------------------------------------------
## library(plotly)
##
## ggplotly(p1)

