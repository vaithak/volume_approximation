
library(volesti)

d=40
m=80
P=GenCube(3,'H')

system.time({ vol = volume(P, nn=d, mm=m, rounding = TRUE) })
