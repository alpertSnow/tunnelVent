## boyancy-driven ventilation
f <- function(u){

rho <- 1.205 # density at 20C, kg/m3
g <- 9.8
H <- 5   # height difference, m
L <- 200  # tunnel length, m
Q <- 139  # heat generation rate, kW/s
U.out <- 10  # wind velocity, m/s
a <- 2.7 # size of duct, m
b <- 2.8 # size of duct, m
mu <- 1.88e-5 # dynamic (absolute) viscosity, kg/(m s)
k <- 0.1 # roughness, m
T.out <- 273+15  # outdoor temp, K
cp <- 1.005 # kJ/(kg K)
T.in <- T.out + Q/(cp * u * rho * a * b)

dT <- abs(T.in - T.out)   # temperature difference, K
d.e <- 1.30 * (a * b)^0.625 / (a + b)^0.25  # effective diameter, m
d.h <- 2 * a * b / (a + b)     # hydraulic diameter, m
Re <- rho * u * d.h / mu # Reynolds number, turbulent: Re > 4000. Turbulent: TRUE

# solve minor pressure loss coefficient, lambda
f1 <- function(lambda){1 / lambda^0.5 + 2 * log(2.51 / (Re * lambda^0.5) + (k / d.h) / 3.72)}
root <- uniroot(f1, c(0,10))
lambda <- root$root
#

P.major.loss <- lambda * L / d.h * (rho * u^2)/2 # major pressure loss, Pa

xi <- 50 # total minor pressure loss coefficient
P.minor.loss <- xi * (rho * u^2)/2

P.loss <- P.major.loss + P.minor.loss

Ps <- rho * g * H * dT / T.in   # boyancy-driven pressure, Pa
Pw <- 0.5 * rho * U.out  # dynamic pressure, Pa


ACH <- u/L*3600 # air change rate per hour
Q.air <- a * b * u * 3600  # air volume rate per hour, m3/h
Ps - P.loss}
print(data.frame(dP, ACH, Q.air, P.loss, dT))


