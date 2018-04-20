# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
#------------------------------------------------------------------------------
# Designed for use with R version R/3.3.2

main <-function() {

    options(scipen=10000)

    # Read input arguments
    args = commandArgs(trailingOnly = TRUE)
    data_file = args[1]
    prefix = args[2]
    
    # Read data
    customTheme = setupRun()
    print("Reading data file..........")
    data = read_delim(data_file, delim = "\t", escape_double = FALSE, 
                      col_names = TRUE, na = "NA")
    #-------------------------------Filter
    filtered_data = filter_transcripts(data)
    gene_and_transcript_venns(filtered_data, "PB36", "PB38", prefix)
    plot_3prime_difference_hist(filtered_data, prefix, customTheme)
    plot_5prime_difference_hist(filtered_data, prefix, customTheme)
    print_known_novel(filtered_data)
    plot_gene_tpm_correlation(filtered_data, "PB36", "PB38", prefix)
    plot_transcript_tpm_correlation(filtered_data, "PB36", "PB38", prefix)
}

# ------------------Functions------------------------------------------

setupRun <- function() {

    # Check for packages and install if not found
    #print("Checking for R packages...")
    #all_packages <- c("ggplot2", "readr", "grid")
    #new_packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]
    #if(length(new_packages)) install.packages(new_packages)

    # Load packages
    library(ggplot2)
    library(readr)
    library(grid)
    library(gridExtra)
    library(reshape2)
    library(gtable)
    library(VennDiagram)
    library(dplyr)

    # Create custom theme for plots
    # axis.text controls tick mark labels
    customTheme = suppressMessages(theme_bw(base_family = "Helvetica", base_size = 14) +
        theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
        theme(plot.title = element_text(lineheight=1, size= 13.5, margin=margin(-10,1,1,1))) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +

        theme(axis.title.x = element_text(color="black", size=18, vjust = -2, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(color="black", vjust=0.75, size=16),
        axis.title.y = element_text(color="black", size=18, margin=margin(0,10,0,0)),
        axis.text.y  = element_text(color="black", vjust=0.75, size=16))  +
        theme(legend.text = element_text(color="black", size = 14), legend.title = element_text(color="black", size=11), legend.key.size = unit(0.5, "cm")))

    return(customTheme)
}

getMedMaxLabel <- function(v) {
    # Returns a print-ready label of the median and max of vector v
    v = as.numeric(v)
    medianV = round(median(v), 3)
    meanV = round(mean(v), 3)
    maxV = max(v)
    minV = min(v)
    if (abs(minV) > maxV) {
        maxV = minV
    }
    print(maxV)
    medianLabel = paste("Median dist =", medianV, sep = " ")
    meanLabel = paste("Mean dist =", meanV, sep = " ")
    maxLabel = paste("Max dist =", maxV, sep = " ")
    plotLabel = paste( medianLabel, meanLabel, maxLabel, sep = "\n")
    return(plotLabel)
}

plot_3prime_difference_hist <- function(data, outprefix, customTheme) {
    # Makes a histogram showing the distribution of 3' end difference
    # sizes (for transcripts with permissive matches)

    filtered = subset(data, annotation_status == "known")
    maxCount = max(table(factor(filtered$diff_3)))
    med_dist = median(filtered$diff_3)
    lab = getMedMaxLabel(filtered$diff_3)
    #q = quantile(abs(filtered$diff_3), probs = 0.75)[1]
    q = 1000

    plotname = paste(outprefix, "3prime_hist.png", sep = "_")
    png(filename = plotname,
     width = 2000, height = 2000, units = "px",
    bg = "white",  res = 300)

    p = ggplot(filtered, aes(diff_3)) + geom_bar(stat="count", fill="dodgerblue4") +
        xlab("Difference at 3' end (bp)") + ylab("Count") + customTheme +
        coord_cartesian(xlim = c(-q,q))  +
        annotate("text", x = q/2, y = maxCount*0.75, label = lab, 
                 color = "black", size = 6)
    print(p)
    dev.off()
}

plot_5prime_difference_hist <- function(data, outprefix, customTheme) {
    # Makes a histogram showing the distribution of 3' end difference
    # sizes (for transcripts with permissive matches)

    filtered = subset(data, annotation_status == "known")
    maxCount = max(table(factor(filtered$diff_5)))
    med_dist = median(filtered$diff_5)
    lab = getMedMaxLabel(filtered$diff_5)
    #q = quantile(abs(filtered$diff_5), probs = 0.75)[1]
    q = 1000

    plotname = paste(outprefix, "5prime_hist.png", sep = "_")
    png(filename = plotname,
     width = 2000, height = 2000, units = "px",
    bg = "white",  res = 300)

    p = ggplot(filtered, aes(diff_5)) + geom_bar(stat="count", fill="dodgerblue4") +
        xlab("Difference at 5' end (bp)") + ylab("Count") + customTheme +
        coord_cartesian(xlim = c(-q,q))  +
        annotate("text", x = q/2, y = maxCount*0.75, label = lab,
                 color = "black", size = 6)
    print(p)
    dev.off()
}

filter_transcripts <- function(data) {
    # Removes transcripts that were seen only once
    transcript_freqs = table(data$transcript_id)
    transcript_freqs_filter = transcript_freqs[transcript_freqs >= 2]
    keep = names(transcript_freqs_filter)
    data_filter = subset(data, transcript_id %in% keep)
    return(data_filter)
}

print_known_novel <- function(data) {
    # Count the number of unique known and novel transcripts in the data
    data = data[!duplicated(data$transcript_id), ]
    print(paste("Known: ", nrow(subset(data, annotation_status == "known")), sep=""))
    print(paste("Novel: ", nrow(subset(data, annotation_status == "novel")), sep=""))

}

gene_and_transcript_venns <- function(data, name1, name2, outprefix) {
    # Makes a set of gene and transcript venn diagrams to compare what was
    # detected in the two input samples

    # Gene venn diagram
    genes1 = unique(subset(data, dataset == name1)$gene_id)
    genes2 = unique(subset(data, dataset == name2)$gene_id)

    plotname = paste(outprefix, "gene_venn.png", sep = "_")
    png(filename = plotname,
     width = 3000, height = 3000, units = "px",
    bg = "white",  res = 300)

    p1 = draw.pairwise.venn(area1=length(genes1), area2=length(genes2), cross.area=length(intersect(genes1,genes2)), category = c(name1, name2), fill = c("navy", "cornflowerblue"), lty = "blank", euler.d = T,scaled = T, label.col="black", fontfamily = rep("Helvetica", 3), cat.fontfamily = rep("Helvetica", 2), alpha = 0.65, cat.pos = c(350,10), cat.dist = c(0.052, 0.052), cat.cex = c(1.5,1.5), cex = rep(1.5,3))

    print(p1)
    dev.off()
    print(paste("Number of unique genes found across datasets: ", length(unique(c(genes1, genes2))), sep = ""))
    #print(length(intersect(genes1,genes2))*100/length(genes1))

    plotname = paste(outprefix, "transcript_venn.png", sep = "_")
    png(filename = plotname,
     width = 3000, height =3000, units = "px",
    bg = "white",  res = 300)

    # Transcript venn diagram
    isoforms1 = unique(subset(data, dataset == name1)$transcript_id)
    isoforms2 = unique(subset(data, dataset == name2)$transcript_id)
    p2 = draw.pairwise.venn(area1=length(isoforms1), area2=length(isoforms2), cross.area=length(intersect(isoforms1, isoforms2)), category = c(name1, name2), fill = c("navy", "cornflowerblue"), lty = "blank", euler.d = T,scaled = T, label.col="black", fontfamily = rep("Helvetica", 3), cat.fontfamily = rep("Helvetica", 2), alpha = 0.65, cat.pos = c(350,10), cat.dist = c(0.052, 0.052), cat.cex = c(1.5,1.5), cex = rep(1.5,3), inverted=T, ext.text = T)
    print(p2)
    dev.off()
    print(length(intersect(isoforms1, isoforms2))*100/length(isoforms1))
    print(paste("Number of unique transcripts found across datasets: ", length(unique(c(isoforms1, isoforms2))), sep = ""))
}

plot_gene_tpm_correlation <- function(data, name1, name2, outprefix) {
    # Compute TPMs of datasets and plot their correlation
    data = subset(data, is.na(transcript_id) == F)
    d1 = subset(data, dataset == name1)
    d2 = subset(data, dataset == name2)

    norm_d1 = 1000000/nrow(d1)
    norm_d2 = 1000000/nrow(d2)

    d1$TPM = norm_d1
    d2$TPM = norm_d2

    tpm1 = d1[, c("gene_id", "transcript_id", "TPM")]
    tpm2 = d2[, c("gene_id", "transcript_id", "TPM")]

    # Sum TPMs by transcript ID
    TPM.per.gene.1 = aggregate(tpm1$TPM, by = list(tpm1$gene_id), sum)
    TPM.per.gene.2 = aggregate(tpm2$TPM, by = list(tpm2$gene_id), sum)
    colnames(TPM.per.gene.1) = c("gene_id", "TPM.1")
    colnames(TPM.per.gene.2) = c("gene_id", "TPM.2")

    # Merge datasets
    m = as.data.frame(merge(TPM.per.gene.1, TPM.per.gene.2, by = "gene_id", all=T))
    m[is.na(m)] = 0
    m$logTPM.1 = log(m$TPM.1 + 1, base=2)
    m$logTPM.2 = log(m$TPM.2 + 1, base=2)
    
    # Correlations: beware- Pearson's correlation changes a lot depending on whether it is run on the looged data or not.
    pearsonCorr = cor.test(~TPM.1 + TPM.2, data=m, method = "pearson", continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~TPM.1 + TPM.2, data=m, method = "spearman", continuity = FALSE, conf.level = 0.95, exact=FALSE)$estimate
    #pearsonCorr = cor.test(m$TPM.1, m$TPM.2, method = "pearson")$estimate
    #spearmanCorr = cor.test(m$TPM.1, m$TPM.2, method = "spearman", exact=FALSE)$estimate

    # Plot 
    plotname = paste(outprefix, "gene_TPM_corr.png", sep = "_")
    xlabel = paste("log2(TPM+1) of gene in ", name1, sep="")
    ylabel = paste("log2(TPM+1) of gene in ", name2, sep="")
    png(filename = plotname,
    width = 2000, height = 2000, units = "px",
        bg = "white",  res = 300)

    p = ggplot(m, aes(x = logTPM.1, y = logTPM.2)) + geom_jitter(alpha = 0.5, col="navy") + theme_bw() + xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=16)) + theme(axis.text.x = element_text(color = "black", size=16), axis.text.y = element_text(color = "black", size=16)) + annotate("text", x = 4, y = 15, label = paste("Pearson r: ", round(pearsonCorr, 3), "\nSpearman rho: ", round(spearmanCorr, 3), sep=""),  color="black", size = 7) + coord_cartesian(xlim=c(0,17), ylim=c(0,17))
    print(p)
    dev.off()
    
}

plot_transcript_tpm_correlation <- function(data, name1, name2, outprefix) {
    # Compute TPMs of datasets and plot their correlation
    data = subset(data, is.na(transcript_id) == F)
    d1 = subset(data, dataset == name1)
    d2 = subset(data, dataset == name2)

    norm_d1 = 1000000/nrow(d1)
    norm_d2 = 1000000/nrow(d2)

    d1$TPM = norm_d1
    d2$TPM = norm_d2

    tpm1 = d1[, c("gene_id", "transcript_id", "TPM")]
    tpm2 = d2[, c("gene_id", "transcript_id", "TPM")]

    # Sum TPMs by transcript ID
    TPM.per.t.1 = aggregate(tpm1$TPM, by = list(tpm1$transcript_id), sum)
    TPM.per.t.2 = aggregate(tpm2$TPM, by = list(tpm2$transcript_id), sum)
    colnames(TPM.per.t.1) = c("transcript_id", "TPM.1")
    colnames(TPM.per.t.2) = c("transcript_id", "TPM.2")

    # Merge datasets
    m = as.data.frame(merge(TPM.per.t.1, TPM.per.t.2, by = "transcript_id", all=T))
    m[is.na(m)] = 0
    m$logTPM.1 = log(m$TPM.1 + 1, base=2)
    m$logTPM.2 = log(m$TPM.2 + 1, base=2)

    # Correlations: beware- Pearson's correlation changes a lot depending on whether it is run on the looged data or not.
    pearsonCorr = cor.test(~TPM.1 + TPM.2, data=m, method = "pearson", continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~TPM.1 + TPM.2, data=m, method = "spearman", continuity = FALSE, conf.level = 0.95, exact=FALSE)$estimate
    #pearsonCorr = cor.test(m$TPM.1, m$TPM.2, method = "pearson")$estimate
    #spearmanCorr = cor.test(m$TPM.1, m$TPM.2, method = "spearman", exact=FALSE)$estimate

    # Plot
    plotname = paste(outprefix, "transcript_TPM_corr.png", sep = "_")
    xlabel = paste("log2(TPM+1) of gene in ", name1, sep="")
    ylabel = paste("log2(TPM+1) of gene in ", name2, sep="")
    png(filename = plotname,
    width = 2000, height = 2000, units = "px",
        bg = "white",  res = 300)
 
    p = ggplot(m, aes(x = logTPM.1, y = logTPM.2)) + geom_jitter(alpha = 0.5, col="navy") + theme_bw() + xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=16)) + theme(axis.text.x = element_text(color = "black", size=16), axis.text.y = element_text(color = "black", size=16)) + annotate("text", x = 4, y = 15, label = paste("Pearson r: ", round(pearsonCorr, 3), "\nSpearman rho: ", round(spearmanCorr, 3), sep=""),  color="black", size = 7) + coord_cartesian(xlim=c(0,17), ylim=c(0,17))
    print(p)
    dev.off()
}


main()
