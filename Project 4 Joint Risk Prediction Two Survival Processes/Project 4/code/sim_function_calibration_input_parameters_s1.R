### input for calibration plots
n.knot <- 3
x.lim <- c(0,0.1)
y.lim <- c(0,0.1)
pred.eval <- seq(0.001, 0.062, 0.0001)
font.size <- 10
linewidths <- c(0.05,0.7)
linecolours <- c(alpha("black", 0.5), "cyan", "red")

### input for ggsave/combined graphs/output files
grid.labels.in <- c("product", "joint-o", "msm", "c-clay", "c-gumb", "c-frank", "f-norm", "f-gam")
label.x.in <- 0.12
label.y.in <- 0.99
font.label.size.in <- 10
dpi.in <- 300
width.in <- 10.8
height.in <- 5

### input for summary/horizontal plots
lwd.horiz <- 0.6
lwd <- 0.6
x.lim <- c(0, max(pred.eval))
x.lim.horiz <- c(0, 0.06)
y.lim.horiz.med <- c(-0.005, 0.015)
y.lim.horiz.p.range <- c(0, 0.1)
linecolours.horiz <- 1:8
xlab <- "Predicted probability"
ylab.med <- "Median calibration"
ylab.p.range <- "Calibration variation"
ylab <- "Observed probability"


### Define analysis methods we want to summarise for
anal.methods.sim <- c("msm", "product", "joint", "c.clay", "c.gumb", "c.frank", "f.normal_weib", "f.gamma_weib")

### The name of the folder contains the results for this scenario (to load data frames from where the data is stored)
scen <- "s1"