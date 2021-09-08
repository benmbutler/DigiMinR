Build instructions for the Digital Soil Mineralogy with R course
================

All of the documentation from the Digital Soil Mineralogy with R course that's hosted online can be built from the source files hosted on GitHub. This `README` file will detail how to go about doing this build.

# Pre-requisites

It is assumed that you have both R and RStudio already installed on your
machine.

# Step 1: Clone the DigiMinR repository from GitHub

Use Git to clone the DigiMinR repository from GitHub to your machine via
<https://github.com/benmbutler/DigiMinR>, or simply download the
repository in ZIP format from the same link.

# Step 2: Install packages from CRAN

The following packages are required to build this **bookdown** project:

``` r
install.packages(c("bookdown",
                   "powdR",
                   "leaflet",
                   "Cubist",
                   "reshape2",
                   "e1071",
                   "gridExtra",
                   "devtools",
                   "downloadthis"))
``` 

# Step 3: Install an archived version of **baseline**

Unfortunately the current version of the **baseline** package (which is used by **powdR**) frequently crashes when fitting backgrounds to XRPD data. The current solution is to install the archived version 1.2.1 from CRAN:

``` r
install.packages("http://cran.r-project.org/src/contrib/Archive/baseline/baseline_1.2-1.tar.gz",
                 repos = NULL,
                 type = "source")
```


# Step 4: Install the **mars2mull** package from GitHub

Now that **devtools** is installed it will be used to download the
**mars2mull** package from GitHub, which contains data required for
Chapter 5 of the Digital Soil Mineralogy with R course

``` r
devtools::install_github("benmbutler/mars2mull")
```

# Step 5: Open the DigiMinR R project

In the cloned or downloaded DigiMinR folder you will find a file called
`DigiMinR.Rproj`. Opening this file should result in the DigiMinR
project being opened in RStudio, which sets your working directory to
that of all of `DigiMinR.Rproj`.

# Step 6: Initial build

Once the project is loaded into RStudio you can build the course
documentation by running:

``` r
bookdown::render_book()
```

The first time you run the build it will take some time (&gt; 1 hour) to
run all of the code required to create the html files that are stored in
`/_book`. Upon building these files you will notice the creation of a
`/_bookdown_files` repository that stores the cache for all data derived
within the course content, allowing subsequent builds can be a lot
faster.

# Step 7: Edit the `.Rmd` files and re-build

The `.Rmd` files used to create the course content are:

1.  `index.Rmd`: The homepage including information on the
    prerequisites, what to expect, acknowledgements etc. This file also
    include various YAML headings used to format the resulting
    documentation.
2.  `01-intro.Rmd`: Chapter 1 - An introduction to loading, plotting and
    manipulating XRPD data in R.
3.  `02-quant.Rmd`: Chapter 2 - Quantitative phase analysis using full
    pattern summation of XRPD data.
4.  `03-machine-learning.Rmd`: Chapter 3 - Using machine learning to
    predict and interpret soil properties from XRPD data.
5.  `04-cluster.Rmd`: Chapter 4 - Using Cluster analysis of XRPD data to
    understand mineral-nutrient relationships.
6.  `05-mars.Rmd`: Chapter 5 - Comparing XRPD data from Mars and
    Scotland.
7.  `999-references.Rmd`: References (no need to edit)

Editing any of the plain text within the `.Rmd` files should result in a
fast re-build of the documentation using `bookdown::render_book()`. However,
if you were to edit any of code within the R chunks, then rebuilds may
take longer because various bits of analysis will have to be recomputed
and cached.
