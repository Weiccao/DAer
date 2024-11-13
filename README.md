# DAer
 Data Augmentation for censored expectile regression

 
## Installation
You can install the **development** version from
[Github](https://github.com/Weiccao/DAer)

```s
# install.packages("remotes")
remotes::install_github("Weiccao/DAer")
```

## Usage

```s
library(dirttee)
library(dplyr)
data("colcancer")
dat <- list()
dat[['delta']] <- colcancer$death
dat[['Y']] <- colcancer$logfollowup
dat[['X']] <- as.matrix(dplyr::select(colcancer, LNE, age))
f <- 'y ~ LNE+age'
res <- DAer(dat, censortype = 'right', formula = f)

```

## License

This package is free and open source software, licensed under GPL-3.
