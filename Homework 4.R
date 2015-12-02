rating <- read.table("http://www.stat.ufl.edu/~winner/data/product_rating.dat", header = F, col.names = c("product", "score"))
attach(rating)
(ybar <- mean(score))
(ybar.trt <- as.vector(tapply(score,product,mean)))
(var.trt <- as.vector(tapply(score,product,var)))
(rep.trt <- as.vector(tapply(score,product,length)))
(num.trt <- length(ybar.trt))
sse <- 0;    sstr <- 0
for (i1 in 1:num.trt) {
sse <- sse + (rep.trt[i1]-1)*var.trt[i1]
sstr <- sstr + rep.trt[i1]*((ybar.trt[i1]-ybar)^2)
}
df.trt <- num.trt-1
df.err <- sum(rep.trt)-num.trt
(mse <- sse/df.err)
(mstr <- sstr/df.trt)
(f_star <- mstr/mse)
(p_f_star <- 1-pf(f_star,df.trt,df.err))
(ybar + c(qt(.025,df.trt)*sqrt(mstr/sum(rep.trt)),qt(.975,df.trt)*sqrt(mstr/sum(rep.trt))))
(c(sse/qchisq(.975,df.err),sse/qchisq(.025,df.err)))
L <- (1/rep.trt[1])*((mstr/mse)*(1/qf(.975,df.trt,df.err))-1)
U <- (1/rep.trt[1])*((mstr/mse)*(1/qf(.025,df.trt,df.err))-1)
(c(L/(1+L),U/(1+U)))
s2.mu <- (mstr-mse)/rep.trt[1]
df.num <- (rep.trt[1]*s2.mu)^2
df.den <- (mstr^2/df.trt) + (mse^2/df.err)
(df <- round(df.num/df.den))
(c(df*s2.mu/qchisq(.975,df),df*s2.mu/qchisq(.025,df)))
horse <- read.table("http://www.stat.ufl.edu/~winner/data/equine_seat.dat", header = F, col.names = c("position", "combo", "stride"))
attach(horse)
position <- factor(position)
combo <- factor(combo)
horse.aov1 <- aov(stride ~ position*combo)
summary(horse.aov1)
(ybar <- mean(stride))
(ybar.pos <- as.vector(tapply(stride,position,mean)))
(ybar.combo <- as.vector(tapply(stride,combo,mean)))
(ybar.trt <- as.vector(tapply(stride,list(position,combo),mean)))
(var.trt <- as.vector(tapply(stride,list(position,combo),var)))
(rep.trt <- as.vector(tapply(stride,list(position,combo),length)))
(b <- length(ybar.combo))
(a <- length(ybar.pos))
(n <- rep.trt[1])
sse <- 0;    sstr <- 0; ssa <- 0;   ssb <- 0
for (i1 in 1:(a*b)) {
sse <- sse + (n-1)*var.trt[i1]
sstr <- sstr + n*((ybar.trt[i1]-ybar)^2)
}
for(i1 in 1:a) {
ssa <- ssa + b*n*(ybar.pos[i1]-ybar)^2
}
for(i1 in 1:b) {
ssb <- ssb + a*n*(ybar.combo[i1]-ybar)^2
}
ssab <- sstr-ssa-ssb
df.a <- a-1; df.b <- b-1;   df.ab <- (a-1)*(b-1)
df.err <- a*b*(n-1)
(mse <- sse/df.err)
(msab <- ssab/df.ab)
(msa <- ssa/df.a)
(msb <- ssb/df.b)
(f_star_ab <- msab/mse)
(p_f_star_ab <- 1-pf(f_star_ab,df.ab,df.err))
(f_star_a <- msa/msab)
(p_f_star_a <- 1-pf(f_star_a,df.a,df.ab))
(f_star_b <- msb/mse)
(p_f_star_b <- 1-pf(f_star_b,df.b,df.err))
(s2b <- (msb/(a*n)) - (mse/(a*n)))
(s2b_df <- round((s2b)^2/(((msb/(a*n))^2/(b-1)) + ((mse/(a*n))^2/(a*b*(n-1))))))
(c(s2b_df*s2b/qchisq(.975,s2b_df),s2b_df*s2b/qchisq(.025,s2b_df)))
q_crit <- 3.635
se_trt_diff <- sqrt(2 * msab/(b*n))
(tukey_hsd <- (q_crit/sqrt(2))*se_trt_diff)
(bon_msd <- qt(1-(0.05/(2*a*(a-1)/2)),(a-1)*(b-1))*se_trt_diff)
for(i1 in 1:(a-1)) {
for (i2 in (i1+1):a) {
print(cbind(i1,i2,(ybar.pos[i1]-ybar.pos[i2])+c(-tukey_hsd,tukey_hsd)))
print(cbind(i1,i2,(ybar.pos[i1]-ybar.pos[i2])+c(-bon_msd,bon_msd)))
}}
clouds <- read.csv("http://www.stat.ufl.edu/~winner/sta4211/clouds.csv",header=TRUE)
attach(clouds)
seeded1 <- subset(clouds, seeded == 1)
unseeded <- subset(clouds, seeded == 0)
t.test(seeded1$testarea, unseeded$testarea, alternative = "two.sided", var.equal= TRUE)
plot(testarea[seeded==1],cntrarea[seeded==1],xlab="Rainfall Amounts in Test Area",ylab="Rainfall Amounts in Control Area",xlim=c(min(testarea[seeded==1]) - .5,max(testarea[seeded==1]) + .5),
      ylim=c(min(cntrarea[seeded==1]) - .5,max(cntrarea[seeded == 1]) + .5),pch=1,main="Analysis of Covariance - Additive Model")
points(testarea[seeded==0],cntrarea[seeded==0],pch=2)
legend("bottomright",,c("Seeded","Unseeded"),pch=c(1,2))
seeded <- factor(seeded)
results <- lm(testarea ~ cntrarea + seeded)
anova(results)
