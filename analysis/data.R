hpc <- read.csv("hpc.csv")
pdo <- read.csv("pdo.csv")
ctemp <- read.csv("ctemp.csv")

colnames(hpc) <- c("year", "total")
hpc.sum.year <- with(hpc,aggregate(total, by=list(year), FUN = sum))
colnames(hpc.sum.year) <- c("Year", "Total HPC Release")


pdo.mean.year <- with(pdo, aggregate(Value, by=list(Year), FUN= mean))
colnames(pdo.mean.year) <- c("Year", "Mean PDO")
pdo.mean.year <- pdo.mean.year[-length(pdo.mean.year$Year),]


if(!require(lubridate)) install.packages('lubridate',repos = "http://cran.us.r-project.org")

ctemp <- na.omit(ctemp)
colnames(ctemp) <- c("date", "temp")
ctemp$date <- mdy(ctemp$date)
ctemp$year <- year(ctemp$date)
ctemp.mean.year <- with(ctemp, aggregate(temp, by=list(year), FUN= mean))
colnames(ctemp.mean.year) <- c("Year", "Mean Creek Temp")
ctemp.mean.year <- ctemp.mean.year[-length(ctemp.mean.year$Year),]


dat <- cbind(hpc.sum.year, pdo.mean.year[,2], ctemp.mean.year[,2])
colnames(dat) <- c("Year","Total HPC Release", "Mean PDO", "Mean Creek Temp" )

if(!require(xlsx)) install.packages('xlsx',repos = "http://cran.us.r-project.org")

write.xlsx(dat, "CLEANDATA.xlsx")

save(dat, file = "CLEANDATA.rdata")

