##################################################
### Code to prepare datasets provided with package
##################################################

### Load the ebmt dataset from mstate package and save it do /data
###
data("ebmt4")
ebmt <- ebmt4

### Create the predicted risks....

### Use the p9 programs to get in the correct form....

### Use this data
usethis::use_data(ebmt, overwrite = TRUE)
usethis::use_data(msebmcal, overwrite = TRUE)
usethis::use_data(ebmcal, overwrite = TRUE)
usethis::use_data(tps0, overwrite = TRUE)
usethis::use_data(tps100, overwrite = TRUE)
