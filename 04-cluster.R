## ---- echo=FALSE, cache = TRUE---------------------------------------------------------------------
if (knitr:::is_html_output())
{
  downloadthis::download_link(
  link = "https://data.mendeley.com/public-files/datasets/r6g94fgx55/files/9b742b95-c071-4e3a-a227-939efa034ac7/file_downloaded",
  button_label = "Download soil property csv",
  button_type = "danger",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)
}


## ---- echo=FALSE, cache = TRUE---------------------------------------------------------------------
if (knitr:::is_html_output())
{
  downloadthis::download_link(
  link = "https://data.mendeley.com/public-files/datasets/r6g94fgx55/files/ef33d276-e39c-4b64-b9aa-f2d1d69c51b3/file_downloaded",
  button_label = "Download zipped XRPD data",
  button_type = "danger",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)
}


## ---- eval = FALSE---------------------------------------------------------------------------------
## #install packages that haven't been used in the course before
## install.packages(c("reshape2", "e1071"))


## ---- message = FALSE, warning = FALSE-------------------------------------------------------------
#load the relevant packages
library(reshape2)
library(e1071)
library(leaflet)
library(powdR)
library(ggplot2)
library(plotly)
library(gridExtra)
library(plyr)


## ---- eval=FALSE-----------------------------------------------------------------------------------
## #Load the soil property data
## props <- read.csv(file = "path/to/your/file/clusters_and_properties.csv")
##
## #Get the full file paths of the XRPD data
## xrpd_paths <-  dir("path/to/xrd", full.names = TRUE)
##
## #Load the XRPD data
## xrpd <- read_xy(files = xrpd_paths)
##
## #Make sure the data are interpolated onto the same 2theta scale
## #as there are small differences within the dataset
## xrpd <- interpolate(xrpd, new_tth = xrpd$icr030336$tth)


## ---- echo=FALSE, message=FALSE, warning=FALSE-----------------------------------------------------
#This will run in the background because echo = FALSE

#Load the soil property data
props <- read.csv(file = "data/clusters_and_properties.csv")

load("data/afsis_xrpd_df.Rdata")

#Convert it back to a multiXY list
xrpd <- as_multi_xy(xrpd)


## ----message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------
#load the afsis reference library
data(afsis)

#Extract a quartz pattern from it
quartz <- data.frame(afsis$tth,
                     afsis$xrd$QUARTZ_1_AFSIS)

#Align the xrpd data to this quartz pattern
#using a restricted 2theta range of 10 to 60
xrpd_aligned <- align_xy(xrpd, std = quartz,
                         xmin = 10, xmax = 60,
                         xshift = 0.2)


## ----message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------
xrpd_aligned_sub <- lapply(xrpd_aligned, subset,
                           tth >= 6)


## ----message=FALSE, warning=FALSE, cache = TRUE----------------------------------------------------
pca <- xrpd_pca(xrpd_aligned_sub,
                mean_center = TRUE,
                root_transform = 2,
                components = 5)

#View the variance explained by the first 5 PCs
pca$eig[1:5,]


## ----pc-fig, fig.cap = "The first 5 PCs plotted against one-another", out.width='100%', fig.align='center', fig.asp=2.3, message=FALSE, warning=FALSE, cache = TRUE----
#Define the x-axis components
x <- c(1, 1, 1, 1,
       2, 2, 2,
       3, 3,
       4)

#Define the y-axis components
y <- c(2, 3, 4, 5,
       3, 4, 5,
       4, 5,
       5)

#Create and empty list
p <- list()

#Populate each item in the list using the dimension defined
#in x and y
for (i in 1:length(x)) {

   p[[i]] <- ggplot(data = pca$coords) +
             geom_point(aes_string(x = paste0("Dim.", x[i]),
                                   y = paste0("Dim.", y[i])),
                        shape = 21,
                        size = 3)

}

grid.arrange(grobs = p,
             ncol = 2)


## ----fig.cap = "The loading of Dim.1.", out.width='100%', fig.align='center', fig.asp=0.75, message=FALSE, warning=FALSE, cache = TRUE----
ggplot(data = pca$loadings) +
  geom_line(aes(x = tth, y = Dim.1)) +
  geom_hline(yintercept = 0)


## ---- fig.cap = "Full pattern summation applied to the loading of Dim.1.", out.width='100%', fig.align='center', fig.asp=0.75, cache = TRUE----
#Load the rockjock library
data(rockjock)

#Merge the rockjock and afsis libraries
rockjock_afsis <- merge(rockjock, afsis)

#All patterns in the library need to be square root transformed
#because this transformation was applied to the soil data
#during the use of xrpd_pca(). In order to avoid errors with the square
#root transforms, any reference pattern with negative counts
#must be removed from the library

remove_index <- which(unlist(lapply(rockjock$xrd, min)) < 0)

rockjock_afsis <- subset(rockjock_afsis,
                         refs = names(rockjock_afsis$xrd)[remove_index],
                         mode = "remove")

#Square root transform the counts
rockjock_afsis_sqrt <- rockjock_afsis
rockjock_afsis_sqrt$xrd <- sqrt(rockjock_afsis_sqrt$xrd)

#Produce a fit using a subset of common soil minerals
dim1_fit <- fps_lm(rockjock_afsis_sqrt,
                   smpl = data.frame(pca$loadings$tth,
                                     pca$loadings$Dim.1),
                   refs = c("Quartz",
                            "Organic matter",
                            "Plagioclase",
                            "K-feldspar",
                            "Goethite",
                            "Illite",
                            "Mica (Di)",
                            "Kaolinite",
                            "Halloysite",
                            "Dickite",
                            "Smectite (Di)",
                            "Smectite (ML)",
                            "Goethite",
                            "Gibbsite",
                            "Amphibole",
                            "Calcite",
                            "Ferrihydrite"),
                     std = "QUARTZ_1_AFSIS",
                     align = 0, #No alignment needed
                     p = 0.01)

plot(dim1_fit, wavelength = "Cu", group = TRUE)

dim1_fit$phases_grouped[order(dim1_fit$phases_grouped$coefficient),]


## ---- fig.cap = "Full pattern summation applied to the loading of Dim.2.", out.width='100%', fig.align='center', fig.asp=0.75, cache = TRUE----
dim2_fit <- fps_lm(rockjock_afsis_sqrt,
                   smpl = data.frame(pca$loadings$tth,
                                     pca$loadings$Dim.2),
                   refs = c("Quartz",
                            "Organic matter",
                            "Plagioclase",
                            "K-feldspar",
                            "Goethite",
                            "Illite",
                            "Mica (Di)",
                            "Kaolinite",
                            "Halloysite",
                            "Dickite",
                            "Smectite (Di)",
                            "Smectite (ML)",
                            "Goethite",
                            "Gibbsite",
                            "Amphibole",
                            "Calcite",
                            "Ferrihydrite"),
                     std = "QUARTZ_1_AFSIS",
                     align = 0, #No alignment needed
                     p = 0.01)

plot(dim2_fit, wavelength = "Cu", group = TRUE)

dim2_fit$phases_grouped[order(dim2_fit$phases_grouped$coefficient),]











#--------------------------------------------------------------------
# Fuzzy clustering
#--------------------------------------------------------------------

#Apply the fuzzy-c-means algorithm to the PCs
fcm <- cmeans(pca$coords[-1],
              center = 9)

#check the data are in the same order
identical(names(fcm$cluster), pca$coords$sample_id)

clusters <- data.frame("SSN" = names(fcm$cluster),
                       "CLUSTER" = paste0("C", unname(fcm$cluster)),
                       pca$coords[-1])

#Reorder the clusters based on Dim.1
#Lowest mean Dim.1 will be Cluster 1
#Highest mean Dim.1 will be Cluster 9
dim1_mean <- aggregate(Dim.1 ~ CLUSTER,
                       data = clusters,
                       FUN = mean)

#Order so that the Dim.1 mean is ascending
dim1_mean <- dim1_mean[order(dim1_mean$Dim.1),]
dim1_mean$NEW_CLUSTER <- paste0("C", 1:nrow(dim1_mean))

#Create a named vector that will be used to revalue cluster names
#Values of the vector are the new values, and old values are the names
rv <- setNames(dim1_mean$NEW_CLUSTER, # the vector values
               dim1_mean$CLUSTER) #the vector names

#use the revalue function to create a new cluster column
clusters$NEW_CLUSTER <- revalue(clusters$CLUSTER,
                                rv)

#Create an empty list
p <- list()

#Populate each item in the list using the dimension already
#defined in x and y above
for (i in 1:length(x)) {

   p[[i]] <- ggplot(data = clusters) +
             geom_point(aes_string(x = paste0("Dim.", x[i]),
                                   y = paste0("Dim.", y[i]),
                                   fill = "NEW_CLUSTER"),
                        shape = 21,
                        size = 3,
                        alpha = 0.5) +
             guides(fill = guide_legend(title="Cluster"))

}

grid.arrange(grobs = p,
             ncol = 2)

p[[1]]

# Fuzzy memberships

#Extract the membership coefficients
members <- data.frame(fcm$membership, check.names = FALSE)

#Add C to all the names to match that used above.
names(members) <- paste0("C", names(members))

#revalue the names based on the new clustering order
names(members) <- revalue(names(members), rv)

#Add membership to the name
names(members) <- paste0("Membership_", names(members))

#Join clusters and members
members <- data.frame(clusters,
                      members)

#Create and empty list for the plots
p <- list()

#Populate each item in the list using the dimension defined
#in x and y
for (i in 1:9) {

   p[[i]] <- ggplot(data = members) +
             geom_point(aes_string(x = "Dim.1",
                                   y = "Dim.2",
                                   fill = paste0("Membership_C",
                                                 i)),
                        shape = 21,
                        size = 3,
                        alpha = 0.5) +
             ggtitle(paste("Cluster", i)) +
             theme(legend.position = "None") +
             scale_fill_gradient(low = "blue",
                                 high = "red")

}

grid.arrange(grobs = p,
             ncol = 3)



#Subsetting samples based on memberships

#Create a blank name to populate with the unique SSNs
member_ssn <- list()

#A loop to omit samples from each cluster with low membership coefficient
for (i in 1:9) {

  memberships <- members[which(members$NEW_CLUSTER == paste0("C", i)),
                         c("SSN", paste0("Membership_C", i))]


  #Extract the samples with top 25 % of membership coefficient for each cluster
  memberships_75 <- which(memberships[[2]] > quantile(memberships[[2]],
                                                      probs = 0.75))

  member_ssn[[i]] <- memberships$SSN[memberships_75]

  names(member_ssn)[i] <- paste0("Cluster.", i)

}

#Unlist the indexes
member_ssn <- unname(unlist(member_ssn))

members_sub <- members[which(members$SSN %in% member_ssn),]

#Plot the results
#Create and empty list
p <- list()

#Populate each item in the list using the dimension defined
#in x and y
for (i in 1:length(x)) {

   p[[i]] <- ggplot(data = members_sub) +
             geom_point(aes_string(x = paste0("Dim.", x[i]),
                                   y = paste0("Dim.", y[i]),
                                   fill = "NEW_CLUSTER"),
                        shape = 21,
                        size = 3,
                        alpha = 0.5) +
             guides(fill = guide_legend(title="Cluster"))

}

grid.arrange(grobs = p,
             ncol = 2)













#---------------------------------------------------------
# Exploring results
#---------------------------------------------------------

#Join clusters to the soil property data
clust_sub <- join(members_sub[c("SSN", "NEW_CLUSTER")],
                  props, by = "SSN")

#Plot cluster 9
leaflet(clust_sub[which(clust_sub$NEW_CLUSTER == "C4"), ]) %>%
  addTiles() %>%
  addCircleMarkers(~Longitude, ~Latitude)



#Create a blank list to populate
cluster_xrpd <- list()

for (i in 1:9) {
cluster_xrpd[[i]] <- xrpd_aligned[clust_sub$SSN[clust_sub$NEW_CLUSTER == paste0("C",i)]]
names(cluster_xrpd)[i] <- paste0("C", i)
}

#Plot cluster 9
plot(as_multi_xy(cluster_xrpd$C4), wavelength = "Cu", normalise = TRUE, int = TRUE)



# Barplots of total nutrient concentrations nutrient
#--------------------------------------------------------------
#Subset the data to the variables of interest
total_nutrients <- subset(clust_sub,
                          select = c(NEW_CLUSTER,
                                     TOC, K, Ca,
                                     Mn, Fe, Ni,
                                     Cu, Zn))

#Create geometric mean function
gmean <- function(x) {exp(mean(log(x)))}

#Compute overall geometric means for each nutrient (i.e. column)
total_nutrients_gmeans <- apply(total_nutrients[-1], 2, gmean)

#Compute the geometric mean for each nutrient by cluster
cluster_gmeans <- by(total_nutrients[-1],
                     as.factor(total_nutrients$NEW_CLUSTER),
                     function(x) apply(x, 2, gmean))

#Bind the data by row
cluster_gmeans <- do.call(rbind, cluster_gmeans)

#Calculate log-ratios by dividing by the geometric mean
#of each nutrient
bp <- apply(cluster_gmeans,1,
            function(x) log(x/total_nutrients_gmeans))


#Use melt from the reshape2 package so that the data are
#in the right format for ggplot2:
bpm <- melt(bp)

#Create a barplot using ggplot2
g1 <- ggplot(bpm,
             aes(fill=Var1, y=value, x=Var2)) +
      geom_bar(position="dodge", stat="identity") +
      theme(legend.title = element_blank()) +
      xlab("Cluster") +
      ylab("")

#Make the barplot interactive
ggplotly(g1)


## ---- fig.cap = "The mean diffractogram of Cluster 9.", out.width='100%', fig.align='center', fig.asp=0.75, cache = TRUE, message = FALSE, warning = FALSE----
#Create a blank list to be populated
xrpd_clusters <- list()

#Calculate the mean diffractogram of each cluster
for (i in 1:9) {

  #Extract the SSN for the cluster
  ssns <- clust_sub$SSN[which(clust_sub$NEW_CLUSTER == paste0("C", i))]

  #Extract the aligned xrpd data and make sure it is a multiXY object
  xrpd_clusters[[i]] <- as_multi_xy(xrpd_aligned[ssns])

  names(xrpd_clusters)[i] <- paste0("C", i)

  #Convert to a data frame
  xrpd_clusters[[i]] <- multi_xy_to_df(xrpd_clusters[[i]],
                                       tth = TRUE)

  #Calculate mean diffractogram
  xrpd_clusters[[i]] <- as_xy(data.frame("tth" = xrpd_clusters[[i]]$tth,
                                         "counts" = rowMeans(xrpd_clusters[[i]][-1])))

}

#Plot the mean diffractogram of Cluster 9
plot(xrpd_clusters$C3, wavelength = "Cu")


#Full pattern summation of each mean pattern

usuals <- c("Quartz", "Kaolinite",
            "Halloysite", "Dickite",
            "Smectite (Di)", "Smectite (ML)",
            "Illite", "Mica (Di)", "Muscovite",
            "K-feldspar", "Plagioclase",
            "Goethite", "Maghemite", "Ilmenite",
            "Hematite", "Gibbsite", "Magnetite",
            "Anatase", "Amphibole", "Pyroxene",
            "Calcite", "Gypsum", "Organic matter",
            "HUMIC_ACID", "FERRIHYDRITE_HUMBUG_CREEK",
            "FERRIHYDRITE",
            "BACK_POS")

#Define the amorphous phases
amorph <- c("ORGANIC_MATTER", "ORGANIC_AFSIS",
            "HUMIC_ACID", "FERRIHYDRITE_HUMBUG_CREEK",
            "FERRIHYDRITE")

clusters_quant <- lapply(xrpd_clusters, afps,
                         lib = rockjock_afsis,
                         std = "QUARTZ_1_AFSIS",
                         refs = usuals,
                         amorphous = amorph,
                         align = 0.2,
                         lod = 0.05,
                         amorphous_lod = 0,
                         force = "BACK_POS")

#Load the regrouping structures for rockjock and afsis
data(rockjock_regroup)
data(afsis_regroup)

#lapply the regrouping structure
clusters_quant_rg <- lapply(clusters_quant,
                            regroup,
                            y = rbind(rockjock_regroup[1:2],
                                      afsis_regroup[1:2]))

#Extract the quantitative data
quant_table <- summarise_mineralogy(clusters_quant_rg,
                                    type = "grouped",
                                    order = TRUE)

#Reduce to the 10 most common phases for a barplot
quant_table <- quant_table[2:11]


#Rename clusters C1 to C9
rownames(quant_table) <- paste0("C", rownames(quant_table))

#Use melt from the reshape2 package so that the data are
#in the right format for ggplot2:
quant_table_m <- melt(as.matrix(quant_table))

#Create a barplot using ggplot2
g2 <- ggplot(quant_table_m,
             aes(fill=Var2, y=value, x=Var1)) +
      geom_bar(position="dodge", stat="identity") +
      theme(legend.title = element_blank()) +
      xlab("Cluster") +
      ylab("")

#Make the barplot interactive
ggplotly(g2)

