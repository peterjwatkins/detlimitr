# detlimitr

The **detlimitr** (estimating *det*ection *limit*s using *R*) package contains functions for use with calibration data in analytical chemistry. In this package, the detection limit can be calculated for either linear, quadratic or power regression.


## Installation

If required, the *detlimiter* package can be downloaded using

``` R
#install.packages("devtools") 
devtools::install_github("peterjwatkins/detlimitr", force=TRUE)
```

## Usage

Some example code:
```R 
library(detlimitr)
data(mtbe)
calcDL(mtbe)  # or delimitr::calcDL(mtbe)
tabulateDL(mtbe) # or delimitr::tabulateDL(mtbe)
plotDL(mtbe)     # or delimitr::plotDL(mtbe)

data("chloromethane")
plotDL(chloromethane, "q")
```

Further information and examples can be found in 'detlimitr.Rmd'.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)
