###
### Testing package functionality
###
a <- 3
b <- 4
rm(list=ls())
load_all()
library(devtools)
data("ebmtcal")
data("msebmtcal")
str(msebmtcal)
data("tps0")
data("tps100")
str(class(msebmtcal))
typeof(msebmtcal)

##########################
### STAB WEIGHTS STUFF ###
##########################

dat.calib.blr.unstab <- calc_calib_blr(data.mstate = msebmtcal,
                                data.raw = ebmtcal[1:500, ],
                                j=1,
                                s=0,
                                t.eval = t.eval,
                                tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))) %>% slice(1:500),
                                curve.type = "rcs",
                                rcs.nk = 3,
                                w.covs = c("year", "agecl", "proph", "match"),
                                CI = 95,
                                CI.R.boot = 200)

dat.calib.blr.stab <- calc_calib_blr(data.mstate = msebmtcal,
                                 data.raw = ebmtcal[1:500, ],
                                 j=1,
                                 s=0,
                                 t.eval = t.eval,
                                 tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))) %>% slice(1:500),
                                 curve.type = "rcs",
                                 rcs.nk = 3,
                                 w.covs = c("year", "agecl", "proph", "match"),
                                 w.stabilised= TRUE,
                                 CI = 95,
                                 CI.R.boot = 200)



plot(dat.calib.blr.unstab[["plotdata"]][[6]]$obs.upper - dat.calib.blr.unstab[["plotdata"]][[6]]$obs.lower, dat.calib.blr.stab[["plotdata"]][[6]]$obs.upper - dat.calib.blr.stab[["plotdata"]][[6]]$obs.lower,
     xlab = "UNSTABILISED", ylab = "STABILISED")
abline(0,1)


testtt3 <- calc_weights(data.mstate = msebmtcal,
                       data.raw = ebmtcal,
                       covs = c("agecl", "year"),
                       j = 1,
                       landmark.type = "all",
                       s = 0,
                       t.eval = 1826,
                       max.weight = 10,
                       stabilised = FALSE)

weights1 <- calc_weights(data.mstate = msebmtcal,
                        data.raw = ebmtcal,
                        covs = c("agecl", "year", "proph", "match"),
                        t.eval = 1826,
                        s = 0,
                        landmark.type = "all",
                        j = 1,
                        max.weight = 10,
                        stabilised = FALSE)
sum(is.na(weights1))

weights2 <- calc_weights(data.mstate = msebmtcal,
                         data.raw = ebmtcal,
                         covs = c("agecl", "year", "proph", "match"),
                         t.eval = t.eval,
                         s = 0,
                         landmark.type = "all",
                         j = 1,
                         max.weight = 10,
                         stabilised = FALSE)


weights3 <- calc_weights(data.mstate = data.mstate,
                        data.raw = data.boot,
                        covs = w.covs,
                        t.eval = t.eval,
                        s = s,
                        landmark.type = w.landmark.type,
                        j = j,
                        max.weight = w.max,
                        stabilised = w.stabilised)


getwd()
t.eval <- 1826
str(ebmtcal)

a <- seq(0.1, 0.9, 0.1)
log(a/(1-a))

load_all()
dat.calib.blr <- calc_calib_blr(data.mstate = msebmtcal,
                                data.raw = ebmtcal,
                                j=1,
                                s=0,
                                t.eval = 1826,
                                tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                curve.type = "rcs",
                                rcs.nk = 3,
                                w.covs = c("year", "agecl", "proph", "match"))
plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
plot.calib_blr(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
ftype(plot)
class(dat.calib.blr)
library(devtools)
install.packages("devtools")
install.packages("usethis")
install.packages("rlang")
install.packages("Rtools")
install.packages("stringr")
attributes(msebmtcal)
rm(list=ls())
devtools::install()
sessionInfo()
mypaths <- .libPaths()
mypaths
.libPaths()
.libPaths(mypaths[2])
install.packages("rlang")
library(calibmsm)

install.packages("rlang")
devtools::document()
rm(list=ls())
load_all()
calc_calib_mlr
devtools::check(vignettes = FALSE)
devtools::install()

load_all()
dat.calib.mlr.j1.s0.s <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      smoother.type = "s",
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")

dat.calib.mlr.j1.s0.sm.ps <- calc_calib_mlr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=0,
                                        t.eval = 1826,
                                        tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        smoother.type = "sm.ps",
                                        ps.int = 3, degree = 2,
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")

dat.calib.mlr.j1.s0.sm.os <- calc_calib_mlr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=0,
                                        t.eval = 1826,
                                        tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        smoother.type = "sm.os",
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")


plot.calib_mlr(dat.calib.mlr.j1.s0.s)
plot.calib_mlr(dat.calib.mlr.j1.s0.sm.ps)
plot.calib_mlr(dat.calib.mlr.j1.s0.sm.os)


dat.calib.mlr.j1.s0 <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")
plot.calib_mlr(dat.calib.mlr.j1.s0)


dat.calib.mlr.j1.s0 <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")
plot.calib_mlr(dat.calib.mlr.j1.s0)




dat.calib.blr.j1.s100 <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=100,
                                        t.eval = 1826,
                                        tp.pred = tps100 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")

dat.calib.blr.j2.s100 <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=2,
                                        s=100,
                                        t.eval = t.eval,
                                        tp.pred = tps100 %>% filter(j == 2) %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                        w.covs = c("year", "agecl", "proph", "match"))

dat.calib.blr.j3.s100 <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=3,
                                        s=100,
                                        t.eval = t.eval,
                                        tp.pred = tps100 %>% filter(j == 3) %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                        w.covs = c("year", "agecl", "proph", "match"))

plot(dat.calib.blr.j1.s100, combine = TRUE, nrow = 2, ncol = 3)
plot(dat.calib.blr.j2.s100, combine = TRUE, nrow = 2, ncol = 3)
plot(dat.calib.blr.j3.s100, combine = TRUE, nrow = 2, ncol = 3)
data("ebmt")
data("ebmtcal")

expect_error(1 / "a")

#   data.mstate <- msebmt
#   data.raw <- ebmt
#   covs <- NULL
#   j
#   landmark.type <- "state"
#   s
#   t.eval <- t.eval
#   max.weight <- 10
#
