cardinals <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Arizona Cardinals.csv", header=T, sep = ",")
falcons <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Atlanta Falcons.csv", header=T, sep = ",")
ravens <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Baltimore Ravens.csv", header=T, sep = ",")
bills <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Buffalo Bills.csv", header=T, sep = ",")
panthers <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Carolina Panthers.csv", header=T, sep = ",")
bears <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Chicago Bears.csv", header=T, sep = ",")
bengals <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Cincinnati Bengals.csv", header=T, sep = ",")
browns <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Cleveland Browns.csv", header=T, sep = ",")
cowboys <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Dallas Cowboys.csv", header=T, sep = ",")
broncos <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Denver Broncos.csv", header=T, sep = ",")
lions <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Detroit Lions.csv", header=T, sep = ",")
packers <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Green Bay Packers.csv", header=T, sep = ",")
texans <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Houston Texans.csv", header=T, sep = ",")
colts <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Indianapolis Colts.csv", header=T, sep = ",")
jaguars <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Jacksonville Jaguars.csv", header=T, sep = ",")
chiefs <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Kansas City Chiefs.csv", header=T, sep = ",")
dolphins <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Miami Dolphins.csv", header=T, sep = ",")
vikings <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Minnesota Vikings.csv", header=T, sep = ",")
patriots <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\New England Patriots.csv", header=T, sep = ",")
saints <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\New Orleans Saints.csv", header=T, sep = ",")
giants <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\New York Giants.csv", header=T, sep = ",")
jets <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\New York Jets.csv", header=T, sep = ",")
raiders <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Oakland Raiders.csv", header=T, sep = ",")
eagles <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Philadelphia Eagles.csv", header=T, sep = ",")
steelers <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Pittsburgh Steelers.csv", header=T, sep = ",")
chargers <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\San Diego Chargers.csv", header=T, sep = ",")
niners <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\San Francisco 49ers.csv", header=T, sep = ",")
seahawks <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Seattle Seahawks.csv", header=T, sep = ",")
rams <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\St. Louis Rams.csv", header=T, sep = ",")
bucs <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Tampa Bay Buccaneers.csv", header=T, sep = ",")
titans <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Tennessee Titans.csv", header=T, sep = ",")
redskins <- read.table("C:\\Users\\Beans\\Desktop\\BRIDGE MAN\\Junior Year\\Spring Semester\\Design of Experiments\\Written Project\\Washington Redskins.csv", header=T, sep = ",")
nfl <- rbind(cardinals, falcons, ravens, bills, panthers, bears, bengals, browns, cowboys, broncos, lions, packers, texans, colts, jaguars, chiefs, dolphins, vikings, patriots, saints, giants, jets, raiders, eagles, steelers, chargers, niners, seahawks, rams, bucs, titans, redskins)
attach(nfl)
nfl <- subset(nfl, Year >= 1970)
Pts <- factor(Pts)
Pts.1 <- factor(Pts.1)
nfl.aov1 <- aov(W ~ Pts * Pts.1)
summary(nfl.aov1)
TukeyHSD(nfl.aov1, "Pts")
TukeyHSD(nfl.aov1, "Pts.1")
interaction.plot(Pts,Pts.1,W, ylim=c(0.7*min(W),1.3*max(W)))
Yds <- factor(Yds)
Yds.1 <- factor(Yds.1)
nfl.aov2 <- aov(W ~ Yds * Yds.1)
summary(nfl.aov2)
TukeyHSD(nfl.aov2, "Yds")
TukeyHSD(nfl.aov2, "Yds.1")
interaction.plot(Yds,Yds.1,W, ylim=c(0.7*min(W),1.3*max(W)))
nfl.aov3 <- aov(W ~ Pts * Pts.1 + Yds*Yds.1)
summary(nfl.aov3)
anova(nfl.aov1, nfl.aov3)
anova(nfl.aov2, nfl.aov3)
