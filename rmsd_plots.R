##A script for plotting root mean square deviation plots
##By Dorothy Nyamai
setwd('F:\\PhD_2018Data\\MD_26_5_2018\\MD_results_4_6_2018\\PfTrpRs_16_10_2018\\PfTrpRs_RMSD')


options(stringsAsFactors=FALSE)
###PfProRs####
PfTrpRs_TYM_table = read.table('PfTrpRs_Holo_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_TYM_dataframe = data.frame(TIME=PfTrpRs_TYM_table$TIME, PfTrpRs_TYM=PfTrpRs_TYM_table$PfTrpRs_TYM, TYM_rmsd=PfTrpRs_TYM_table$TYM_rmsd)
colnames(PfTrpRs_TYM_table) = c("TIME","PfTrpRs_TYM","TYM_rmsd")

time_PfTrpRs_TYM <- PfTrpRs_TYM_dataframe$TIME
rmsd_PfTrpRs_bb <- PfTrpRs_TYM_dataframe$PfTrpRs_TYM
rmsd_TYM <- PfTrpRs_TYM_dataframe$TYM_rmsd 

png(filename ="PfTrpRs_rmsd_comb_ed.png",width = 18, height = 10, units = 'in', res=300)

attach(mtcars)
par(mfrow=c(3,3))
par(mar=c(5.1, 5.1, 4.1, 2.1))


par(mar=c(4.1, 5.1, 1.6, 0.4))
plot(time_PfTrpRs_TYM,rmsd_PfTrpRs_bb, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.6, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_TYM,rmsd_TYM, col="blue", lwd=1.5, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM")),col=c("black","blue"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns


##PfTrpRs_195


#PfTrpRs_195_table = read.table('PfTrpRs_195_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


#PfTrpRs_195_dataframe = data.frame(TIME=PfTrpRs_195_table$TIME, PfTrpRs_195=PfTrpRs_195_table$PfTrpRs_195, TYM_rmsd=PfTrpRs_195_table$TYM_rmsd, LIG_rmsd=PfTrpRs_195_table$LIG_rmsd)
#colnames(PfTrpRs_195_table) = c("TIME","PfTrpRs_195","TYM_rmsd", "LIG_rmsd")

#time_PfTrpRs_195 <- PfTrpRs_195_dataframe$TIME
#rmsd_PfTrpRs_195 <- PfTrpRs_195_dataframe$PfTrpRs_195
#rmsd_TYM <- PfTrpRs_195_dataframe$TYM_rmsd 
#rmsd_LIG <- PfTrpRs_195_dataframe$LIG_rmsd 

#par(mar=c(4.1, 5.1, 1.6, 0.4))
#plot(time_PfTrpRs_195,rmsd_PfTrpRs_195, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.6, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
#axis(side = 1, seq(0,200,by=10), tck = -0.02)
#lines(time_PfTrpRs_195,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
#lines(time_PfTrpRs_195,rmsd_LIG, col="red", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
#legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("SANC195")),col=c("black","blue","red"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns


##PfTrpRs_257


PfTrpRs_257_table = read.table('PfTrpRs_257_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_257_dataframe = data.frame(TIME=PfTrpRs_257_table$TIME, PfTrpRs_257=PfTrpRs_257_table$PfTrpRs_257, TYM_rmsd=PfTrpRs_257_table$TYM_rmsd, LIG_rmsd=PfTrpRs_257_table$LIG_rmsd)
colnames(PfTrpRs_257_table) = c("TIME","PfTrpRs_257","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_257 <- PfTrpRs_257_dataframe$TIME
rmsd_PfTrpRs_257 <- PfTrpRs_257_dataframe$PfTrpRs_257
rmsd_TYM <- PfTrpRs_257_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_257_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.6, 0.4))
plot(time_PfTrpRs_257,rmsd_PfTrpRs_257, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.6, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_257,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_257,rmsd_LIG, col="green", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("SANC257")),col=c("black","blue","green"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns


##PfTrpRs_438
PfTrpRs_438_table = read.table('PfTrpRs_438_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_438_dataframe = data.frame(TIME=PfTrpRs_438_table$TIME, PfTrpRs_438=PfTrpRs_438_table$PfTrpRs_438, TYM_rmsd=PfTrpRs_438_table$TYM_rmsd, LIG_rmsd=PfTrpRs_438_table$LIG_rmsd)
colnames(PfTrpRs_438_table) = c("TIME","PfTrpRs_438","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_438 <- PfTrpRs_438_dataframe$TIME
rmsd_PfTrpRs_438 <- PfTrpRs_438_dataframe$PfTrpRs_438
rmsd_TYM <- PfTrpRs_438_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_438_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.6, 0.4))
plot(time_PfTrpRs_438,rmsd_PfTrpRs_438, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.6, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_438,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_438,rmsd_LIG, col="magenta", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("SANC438")),col=c("black","blue","magenta"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns

##PfTrpRs_465
PfTrpRs_465_table = read.table('PfTrpRs_465_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_465_dataframe = data.frame(TIME=PfTrpRs_465_table$TIME, PfTrpRs_465=PfTrpRs_465_table$PfTrpRs_465, TYM_rmsd=PfTrpRs_465_table$TYM_rmsd, LIG_rmsd=PfTrpRs_465_table$LIG_rmsd)
colnames(PfTrpRs_465_table) = c("TIME","PfTrpRs_465","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_465 <- PfTrpRs_465_dataframe$TIME
rmsd_PfTrpRs_465 <- PfTrpRs_465_dataframe$PfTrpRs_465
rmsd_TYM <- PfTrpRs_465_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_465_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.6, 0.4))
plot(time_PfTrpRs_465,rmsd_PfTrpRs_465, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.6, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_465,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_465,rmsd_LIG, col="chocolate", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("SANC465")),col=c("black","blue","chocolate"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns

##ZINC

setwd('F:\\PhD_2018Data\\MD_26_5_2018\\MD_results_4_6_2018\\MD_Zinc_17_9_2019\\MD_Zinc_3_10_2018\\PfTrpRs_Zinc\\PfTrpRs_RMSD')



##PfTrpRs_108


PfTrpRs_108_table = read.table('PfTrpRs_108_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_108_dataframe = data.frame(TIME=PfTrpRs_108_table$TIME, PfTrpRs_108=PfTrpRs_108_table$PfTrpRs_108, TYM_rmsd=PfTrpRs_108_table$TYM_rmsd, LIG_rmsd=PfTrpRs_108_table$LIG_rmsd)
colnames(PfTrpRs_108_table) = c("TIME","PfTrpRs_108","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_108 <- PfTrpRs_108_dataframe$TIME
rmsd_PfTrpRs_108 <- PfTrpRs_108_dataframe$PfTrpRs_108
rmsd_TYM <- PfTrpRs_108_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_108_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.8, 0.4))
plot(time_PfTrpRs_108,rmsd_PfTrpRs_108, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.4, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_108,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_108,rmsd_LIG, col="darkviolet", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("PC124921071")),col=c("black","blue","darkviolet"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns




##PfTrpRs_161
PfTrpRs_161_table = read.table('PfTrpRs_161_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_161_dataframe = data.frame(TIME=PfTrpRs_161_table$TIME, PfTrpRs_161=PfTrpRs_161_table$PfTrpRs_161, TYM_rmsd=PfTrpRs_161_table$TYM_rmsd, LIG_rmsd=PfTrpRs_161_table$LIG_rmsd)
colnames(PfTrpRs_161_table) = c("TIME","PfTrpRs_161","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_161 <- PfTrpRs_161_dataframe$TIME
rmsd_PfTrpRs_161 <- PfTrpRs_161_dataframe$PfTrpRs_161
rmsd_TYM <- PfTrpRs_161_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_161_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.8, 0.4))
plot(time_PfTrpRs_161,rmsd_PfTrpRs_161, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.4, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_161,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_161,rmsd_LIG, col="darkgreen", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("PC124836483")),col=c("black","blue","darkgreen"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns

##PfTrpRs_209
PfTrpRs_209_table = read.table('PfTrpRs_209_rmsd.csv',sep=',', header = TRUE, quote="'", stringsAsFactors=FALSE, fill=TRUE)


PfTrpRs_209_dataframe = data.frame(TIME=PfTrpRs_209_table$TIME, PfTrpRs_209=PfTrpRs_209_table$PfTrpRs_209, TYM_rmsd=PfTrpRs_209_table$TYM_rmsd, LIG_rmsd=PfTrpRs_209_table$LIG_rmsd)
colnames(PfTrpRs_209_table) = c("TIME","PfTrpRs_209","TYM_rmsd", "LIG_rmsd")

time_PfTrpRs_209 <- PfTrpRs_209_dataframe$TIME
rmsd_PfTrpRs_209 <- PfTrpRs_209_dataframe$PfTrpRs_209
rmsd_TYM <- PfTrpRs_209_dataframe$TYM_rmsd 
rmsd_LIG <- PfTrpRs_209_dataframe$LIG_rmsd 

par(mar=c(4.1, 5.1, 1.8, 0.4))
plot(time_PfTrpRs_209,rmsd_PfTrpRs_209, type="l", xaxt='n', cex.lab=1.8, cex.axis=1.4, ylim=c(0, 0.5), lwd=1.3, pch=19,col='black',xlab='Time (ns)',ylab=expression(paste("Distance (nm)")))
axis(side = 1, seq(0,200,by=10), tck = -0.02)
lines(time_PfTrpRs_209,rmsd_TYM, col="blue", lwd=1.3, type="l", pch=19)
lines(time_PfTrpRs_209,rmsd_LIG, col="deeppink4", lwd=1.3, type="l", pch=19)
##abline(v=207, col="azure3",lwd=1, lty=2)
##abline(v=215, col="azure3",lwd=1, lty=2)
##abline(v=225, col="azure3",lwd=1, lty=2)
legend("topright",legend=c(expression("Protein backbone"), expression("TYM"), expression("PC99601290")),col=c("black","blue","deeppink4"), lty=1:1, cex=1.2, lwd=5, bty="n") ##, ncol=2 use this if you want the legent to be in two columns


dev.off()

