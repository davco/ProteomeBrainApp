# ProteomeBrainApp 

ProteomeBrainApp is an interactive app built in R, showing protein expression across human brain tissues from baseline and different stages of Alzheimer disease.


## Features

* Currently there are 3 reanalyzed proteomics experiments embedded within the app.
* It shows the protein expression values across baseline human brain tissues and from different Braak stages of Alzheimer disease.
* The protein expression values are shown in color over a MRI scan using the R package EveTemplate.
* Filters the protein expression by selecting the gene ID and the biological condition.

![Capture](https://user-images.githubusercontent.com/15140798/161384015-4606784c-88bd-4957-9e40-657da498b483.PNG)


## Installation

To install the required packages to run the app

```r
install.packages(c("shiny", "shinythemes", "dplyr", "readr", "stringi", "data.table"))

source("https://neuroconductor.org/neurocLite.R")
neuro_install('EveTemplate')

```

## Getting Started

Once installed, load the library and run the app with:

```r
library(shiny)
runGitHub("ProteomeBrainApp", "davco")
```


## License

The shiny package as a whole is licensed under the GPLv3.

