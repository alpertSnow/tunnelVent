## boyancy-driven ventilation

# key parameters!!
H <- 3   # height difference, m
L <- 400  # tunnel length, m
Q <- 50  # heat generation rate, kW
a <- 2 # size of duct, m
b <- 2.4 # size of duct, m
T.out <- 273+15  # outdoor temp, K
xi <- 100 # total minor pressure loss coefficient
U.out <- 2  # wind velocity, m/s
#

# properties 
rho <- 1.205 # density at 20C, kg/m3
g <- 9.8
mu <- 1.88e-5 # dynamic (absolute) viscosity, kg/(m s)
k <- 0.1 # roughness, m
cp <- 1.005 # kJ/(kg K)
#

f <- function(u){
        # calculated parameters
        T.in <- T.out + Q/(cp * u * rho * a * b)
        dT <- abs(T.in - T.out)   # temperature difference, K
        d.e <- 1.30 * (a * b)^0.625 / (a + b)^0.25  # effective diameter, m
        d.h <- 2 * a * b / (a + b)     # hydraulic diameter, m
        Re <- rho * u * d.h / mu # Reynolds number, turbulent: Re > 4000. Turbulent: TRUE
        #
        
        # solve minor pressure loss coefficient, lambda
        f1 <- function(lambda){1 / lambda^0.5 + 2 * log(2.51 / (Re * lambda^0.5) + (k / d.h) / 3.72)}
        lambda <- uniroot(f1, c(0,10))$root
        #

        P.major.loss <- lambda * L / d.h * (rho * u^2)/2 # major pressure loss, Pa
        P.minor.loss <- xi * (rho * u^2)/2
        
        P.loss <- P.major.loss + P.minor.loss
        
        Ps <- rho * g * H * dT / T.in   # boyancy-driven pressure, Pa
        Pw <- 0.5 * rho * U.out  # dynamic pressure, Pa
        Ps - P.loss
}

u <- uniroot(f, c(1e-5, 100))$root
T.in <- T.out + Q/(cp * u * rho * a * b)
dT <- abs(T.in - T.out)   # temperature difference, K
Ps <- rho * g * H * abs(T.in - T.out) / T.in   # boyancy-driven pressure, Pa
ACH <- u/L*3600 # air change rate per hour
Q.air <- a * b * u * 3600  # air volume rate per hour, m3/h
print(data.frame(u, ACH, Q.air, Ps, dT))
