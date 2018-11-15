library(ballgown)
data_directory ="./"
HCC_expr = ballgown(dataDir=data_directory, samplePattern='HK', meas='all')
save(HCC_expr, file='HCC_expr.rdata')

