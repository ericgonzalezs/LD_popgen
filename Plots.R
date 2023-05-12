library(data.table)
library(ggplot2)
library(ggpubr)

        centromers <- read.table("Centromeres_LC075745_HA412_no_NA.bed", head=F)
        centromers <- centromers[order(centromers[,1], centromers[,2]), ]
        centromers$start <- pmin(centromers$V2, centromers$V3)
        centromers$end <- pmax(centromers$V2, centromers$V3)

        centromers_1 <- subset(centromers, centromers$V1=="Ha412HOChr01")
        centromers_2 <- subset(centromers, centromers$V1=="Ha412HOChr02")
        centromers_3 <- subset(centromers, centromers$V1=="Ha412HOChr03")
        centromers_4 <- subset(centromers, centromers$V1=="Ha412HOChr04")
        centromers_5 <- subset(centromers, centromers$V1=="Ha412HOChr05")
        centromers_6 <- subset(centromers, centromers$V1=="Ha412HOChr06")
        centromers_7 <- subset(centromers, centromers$V1=="Ha412HOChr07")
        centromers_8 <- subset(centromers, centromers$V1=="Ha412HOChr08")
        centromers_9 <- subset(centromers, centromers$V1=="Ha412HOChr09")
        centromers_10 <- subset(centromers, centromers$V1=="Ha412HOChr10")
        centromers_11 <- subset(centromers, centromers$V1=="Ha412HOChr11")
        centromers_12 <- subset(centromers, centromers$V1=="Ha412HOChr12")
        centromers_13 <- subset(centromers, centromers$V1=="Ha412HOChr13")
        centromers_14 <- subset(centromers, centromers$V1=="Ha412HOChr14")
        centromers_15 <- subset(centromers, centromers$V1=="Ha412HOChr15")
        centromers_16 <- subset(centromers, centromers$V1=="Ha412HOChr16")
        centromers_17 <- subset(centromers, centromers$V1=="Ha412HOChr17")

       #inv
Inv1 <- data.frame(start=c(6000000), end=c(10000000))
Inv2 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv3 <- data.frame(start=c(1000000000000), end=c(100000000000000))
Inv4 <- data.frame(start=c(1000000000000), end=c(100000000000000))
Inv5 <- data.frame(start=c(148000000), end=c(177000000))
Inv6 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv7 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv8 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv9 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv10 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv11 <- data.frame(start=c(19000000), end=c(54000000))
Inv12 <- data.frame(start=c(1000000000000), end=c(10000000000000))
Inv13 <- data.frame(start=c(9000000, 139000000), end=c(110000000, 158000000))
Inv14 <- data.frame(start=c(101000000, 129000000), end=c(129000000, 135000000))
Inv15 <- data.frame(start=c(104000000), end=c(176000000))
Inv16 <- data.frame(start=c(10000000, 147000000), end=c(23000000, 159000000))
Inv17 <- data.frame(start=c(188000000), end=c(197000000))



        a = 1
        args <- commandArgs(trailingOnly = TRUE)

        for(f in args) {

                    L <- sub("#", "", readLines(f))

                    chr <- read.table(text = L, skip = 2, header = TRUE)
                    chr$lengt <- chr$right_snp - chr$left_snp
                    chr$Pperbase <- chr$p0.500/chr$lengt
                    chr$POS <- (chr$right_snp + chr$left_snp)/2
                    chr$kb_lengt <- chr$lengt/1000
                    chr$Pperkb  <- chr$p0.500/chr$kb_lengt
                    chr$log10_Pperbase <- log10(chr$Pperbase)
                    chr$log10_kb <- log10(chr$Pperkb)

                    name <- gsub(".txt", "", f)
                    centro <- paste("centromers_", a, sep="")
                    INV <- paste("Inv", a, sep="")

                    PL <- ggplot(chr, aes(x=POS, y=Pperbase)) + geom_point() +  geom_vline(xintercept= eval(as.name(centro))$start, color = "purple", alpha = .05) +
                    geom_vline(xintercept= eval(as.name(INV))$start, color = "green") + geom_vline(xintercept= eval(as.name(INV))$end, color = "green") +
                    xlim(0, max(chr$POS)) + ylim(0, max(chr$Pperbase)) +
                    annotate("rect", xmin = eval(as.name(INV))$start, xmax = eval(as.name(INV))$end, ymin = 0, ymax = max(chr$Pperbase), alpha = .2, fill = "green") +
                    ggtitle(paste("Chr", a, sep=" ")) +
                    theme(axis.text=element_text(size=16,face="bold") , axis.title=element_text(size=16,face="bold"), plot.title = element_text(size = 20, face = "bold") )

                    #geom_rect(aes(xmin=eval(as.name(INV))$start, xmax=eval(as.name(INV))$end, ymin=0, ymax=Inf, color = "green", alpha = .05))

                    assign(paste("CAL_plot_N", a, sep=""),PL)

                    assign(paste("chr_CAL_N", a, sep=""),chr)

                    spline <- smooth.spline(chr$POS, chr$Pperbase, spar=0.35)

                    assign(paste("spline_CAL_N", a, sep=""),spline)

                    print(a)
                   #print(mean(chr$Pperbase))
                    print(summary(chr$Pperbase))
                    a=a + 1
                }


#jpeg("All_CHR.jpg", width=2000, height= 2000)
#ggarrange(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7, plot_8, plot_9, plot_10, plot_11, plot_12, plot_13, plot_14, plot_15, plot_16, plot_17, ncol = 2, nrow =9)
#dev.off()

#jpeg("spline_example.jpg", width=2000, height= 2000)
#plot(chr_newmut_ancmatrix_1$POS, chr_newmut_ancmatrix_1$Pperbase, type="n")
#lines(spline_newmut_ancmatrix_1, col='blue', lwd=2)
#dev.off()

#jpeg("log10p_per_kb_example.jpg", width=2000, height= 2000)
#plot(chr_1$POS, chr_1$log10_kb)
#lines(chr_1$POS, chr_1$log10_kb, col='blue', lwd=2)
#dev.off()

#jpeg("log10p_per_b_example.jpg", width=2000, height= 2000)
#plot(chr_1$POS, chr_1$log10_Pperbase)
#dev.off()

save.image("Plots_cal_2M_map_TE_N.Rdata")


############

