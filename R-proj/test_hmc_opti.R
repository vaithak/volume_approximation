
library(volesti)
rets = opti_hmc(60,60,70,500 + 4*3600)


save(rets, file = "sdp_hmc_hnr3.RData")
