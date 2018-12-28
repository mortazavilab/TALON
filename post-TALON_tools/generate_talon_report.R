# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# generate_talon_report.R is a utility that creates an R markdown and PDF
# report illustrating the results of a single TALON dataset run, or alternately,
# a pair of datasets (i.e. biological replicates).

main <- function() {

    load_packages()
    custom_theme <- custom_theme()
    opt <- parse_options()

    database <- opt$database
    outdir <- opt$outdir
    datasets <- str_split(opt$datasets, ",")[[1]]
    whitelists <- process_whitelists(database, datasets, opt$whitelists)
    message("--------------------------") 

    message("Number of reads successfully annotated in dataset ", datasets[1],
            ": ", count_reads_processed(database, datasets[1]), sep="")
    message("Number of known genes detected at least once in dataset ", 
            datasets[1], ": ", count_known_genes(database, datasets[1]), sep="")
    message("Number of known transcripts detected at least once in dataset ",
            datasets[1], ": ", count_known_transcripts(database, datasets[1]), sep="")
    message("--------------------------")

    plot_read_length_distribution(database, datasets[1], 
                                  whitelists[[1]]$transcript_ID,
                                  custom_theme, outdir)

    plot_isoforms_per_gene(whitelists[[1]], custom_theme, outdir)

    plot_exons_per_read(database, datasets[1],
                                   whitelists[[1]]$transcript_ID,
                                   custom_theme, outdir)

    plot_TPM_distributions(database, datasets[1], whitelists[[1]]$transcript_ID,
                           custom_theme, outdir)

}

process_whitelists <- function(database, datasets, whitelist_input) {
    # If the user provided whitelists, read into a list data structure.
    # If whitelists are missing, then create them by pulling all of the 
    # transcript and gene IDs from the database.

    if (is.null(whitelist_input)) {
        whitelists <- lapply(datasets, fetch_all_genes_and_transcripts, database)           

    } else {
        whitelist_files <- str_split(whitelist_input, ",")[[1]]
        whitelists <- lapply(whitelist_files, read_whitelists)
 
    }

    if (length(datasets) != length(whitelists)) {
           stop("Error: If providing whitelists, make sure that the quantity matches the number of dataset names you are providing.")
    }

    return(whitelists)
}

read_whitelists <- function(whitelist_file) {
    # Read a whitelist from a file
    w <- as.data.frame(read_delim(whitelist_file, ",", escape_double = FALSE,
                            col_names = FALSE, trim_ws = TRUE, na = "NA"))
    colnames(w) <- c("gene_ID", "transcript_ID")   
    return(w)
}

fetch_all_genes_and_transcripts <- function(dataset, database) {
    # Given a dataset, get a table of all genes and transcripts detected
    # (no filtering)

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Count number of reads on record for this dataset
    query_text <- paste("SELECT t.gene_ID,
                                abd.transcript_ID
                         FROM abundance AS abd
                         LEFT JOIN transcripts as t
                             ON t.transcript_ID = abd.transcript_ID
                             AND dataset =", dataset, "", sep="'")

    query <- dbSendQuery(con, query_text)
    genes_and_transcripts <- dbFetch(query, n = -1)

    dbClearResult(query)
    dbDisconnect(con)

    return(genes_and_transcripts)


}

count_reads_processed <- function(database, dataset) {
    # Use the 'observed' table in the database to find out how many sam reads
    # from the given dataset made it into the database. No filtering.

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Count number of reads on record for this dataset
    query_text <- paste("SELECT COUNT(read_name) FROM observed WHERE dataset =", 
                        dataset, "", sep="'")
                       
    query <- dbSendQuery(con, query_text)
    n_reads <- dbFetch(query, n = -1)

    dbClearResult(query)
    dbDisconnect(con)

    return(n_reads)
}

count_known_genes <- function(database, datasets) {
    # Count the number of unique known genes detected in the listed datasets.
    # Only known full splice match transcripts get counted towards this total.

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch known genes
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
    query_text <- paste("SELECT t.gene_ID,
                                ga.value AS gene_annot,
                                ta.value AS transcript_annot
                           FROM abundance AS abd
                           LEFT JOIN transcripts as t 
                               ON t.transcript_ID = abd.transcript_ID
                               AND abd.dataset IN ", dataset_str, " ",
                          "LEFT JOIN gene_annotations as ga  
                               ON t.gene_ID = ga.ID
                               AND ga.attribute = 'gene_status'
                               AND ga.value = 'KNOWN'
                           LEFT JOIN transcript_annotations as ta
                               ON t.transcript_ID = ta.ID
                               AND ta.attribute = 'transcript_status'
                               AND ta.value = 'KNOWN'",
                           sep="")
   
    query <- dbSendQuery(con, query_text)
    genes <- as.data.frame(dbFetch(query, n = -1))
   
    # Filter to remove rows that have a null gene or transcript status col.
    genes <- unique(genes[complete.cases(genes), "gene_ID"])

    dbClearResult(query)
    dbDisconnect(con)
 
    return(length(genes))

}

count_known_transcripts <- function(database, datasets) {
    # Count the number of unique known transcripts detected in the listed datasets

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch known genes
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
    query_text <- paste("SELECT abd.transcript_ID,
                                ta.value AS transcript_annot
                           FROM abundance AS abd
                           LEFT JOIN transcript_annotations as ta
                               ON abd.transcript_ID = ta.ID
                               AND abd.dataset IN ", dataset_str, " ",
                              "AND ta.attribute = 'transcript_status'
                               AND ta.value = 'KNOWN'",
                           sep="")

    query <- dbSendQuery(con, query_text)
    transcripts <- as.data.frame(dbFetch(query, n = -1))

    # Filter to remove rows that have a null gene or transcript status col.
    transcripts <- unique(transcripts[complete.cases(transcripts), "transcript_ID"])

    dbClearResult(query)
    dbDisconnect(con)

    return(length(transcripts))

}


plot_read_length_distribution <- function(database, datasets, whitelist, 
                                          custom_theme, outdir) {
    # Use the read lengths recorded in the 'observed' table to plot a read length
    # distribution. In the query, grab the annotation status so that I can use 
    # it in future versions if I want 

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch read lengths of observed transcripts, along with the 
    # annotation status
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
    whitelist_str <- paste("('", paste(whitelist, collapse = "','"), "')", sep="")
    query_text <- paste("SELECT obs.read_name,
                                obs.transcript_ID,
                                obs.dataset,
                                ta.value AS annot,
                                obs.read_length
                           FROM observed AS obs
                           LEFT JOIN transcript_annotations as ta
                               ON obs.transcript_ID = ta.ID
                               AND obs.dataset IN ", dataset_str, " ",
                              "AND ta.attribute = 'transcript_status'
                               AND ta.value = 'KNOWN'
                               AND obs.transcript_ID IN", whitelist_str,
                           sep="")

    query <- dbSendQuery(con, query_text)
    reads <- as.data.frame(dbFetch(query, n = -1))
    reads$annot[is.na(reads$annot)] <- "NOVEL"
    reads$read_length <- reads$read_length/1000
    reads$annot <- as.factor(reads$annot)

    dbClearResult(query)
    dbDisconnect(con)

    # Plotting and output settings
    fname <- paste(outdir, "/transcript_length_by_annotation_status.png", sep="")
    xlabel <- "Transcript length (kilobases)"
    ylabel <- "Count"
    colors <- c("skyblue", "olivedrab3")

    png(filename = fname,
        width = 3000, height = 2000, units = "px",
        bg = "white",  res = 300)    

    g = ggplot(reads, aes(read_length, color = annot, fill = annot)) + 
               geom_histogram(alpha = 0.5, position="identity") + 
               xlab(xlabel) + ylab(ylabel) + 
               scale_colour_manual(values = c(colors)) +
               scale_fill_manual(values = c(colors)) +
               custom_theme
    print(g)   
 
    dev.off()
}

plot_exons_per_read <- function(database, datasets, whitelist,
                                          custom_theme, outdir) {
    # Fetch the reads in the 'observed' table and join to transcript table in
    # order to plot the distribution of exons per read. In the query, grab the 
    # annotation status as well.

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch n_exons of observed transcripts, along with the
    # annotation status
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
    whitelist_str <- paste("('", paste(whitelist, collapse = "','"), "')", sep="")

    query_text <- paste("SELECT obs.read_name,
                                obs.transcript_ID,
                                obs.dataset,
                                ta.value AS annot,
                                t.n_exons
                           FROM observed AS obs
                           LEFT JOIN transcript_annotations as ta
                               ON obs.transcript_ID = ta.ID
                               AND obs.dataset IN ", dataset_str, " ",
                              "AND ta.attribute = 'transcript_status'
                               AND ta.value = 'KNOWN'
                               AND obs.transcript_ID IN", whitelist_str,
                          "LEFT JOIN transcripts as t
                               ON t.transcript_ID = obs.transcript_ID",
                           sep="")

    query <- dbSendQuery(con, query_text)
    reads <- as.data.frame(dbFetch(query, n = -1))

    reads$annot[is.na(reads$annot)] <- "NOVEL"
    reads$annot <- as.factor(reads$annot)

    dbClearResult(query)
    dbDisconnect(con)

    # Plotting and output settings
    fname <- paste(outdir, "/exons_per_read_by_annotation_status.png", sep="")
    xlabel <- "Number of exons in read"
    ylabel <- "Count"
    colors <- c("skyblue", "olivedrab3")

    png(filename = fname,
        width = 3000, height = 2000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(reads, aes(n_exons, color = annot, fill = annot)) +
               geom_histogram(alpha = 0.5, position="identity") +
               xlab(xlabel) + ylab(ylabel) +
               scale_colour_manual(values = c(colors)) +
               scale_fill_manual(values = c(colors)) +
               custom_theme
    print(g)

    dev.off()
}

plot_TPM_distributions <- function(database, datasets, whitelist,
                                          custom_theme, outdir) {

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch information from abundance table
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
    whitelist_str <- paste("('", paste(whitelist, collapse = "','"), "')", sep="")

    query_text <- paste("SELECT t.gene_ID,
                                abd.transcript_ID,
                                abd.count
                           FROM abundance AS abd
                           LEFT JOIN transcripts as t
                               ON t.transcript_ID = abd.transcript_ID
                               AND abd.dataset IN ", dataset_str, " ",
                              "AND abd.transcript_ID IN", whitelist_str,
                           sep="")

    query <- dbSendQuery(con, query_text)
    abundance <- as.data.frame(dbFetch(query, n = -1))

    dbClearResult(query)
    dbDisconnect(con)

    # Compute gene and transcript TPMs
    transcript_TPMs <- abundance
    transcript_TPMs$TPM <- (transcript_TPMs$count)*1000000/sum(transcript_TPMs$count)

    gene_TPMs <- aggregate(transcript_TPMs$TPM, by=list(transcript_TPMs$gene_ID), FUN=sum)
    colnames(gene_TPMs) <- c("gene_ID", "TPM")
    transcript_TPMs$TPM <- log(transcript_TPMs$TPM + 1, base = 2)
    gene_TPMs$TPM <- log(gene_TPMs$TPM + 1, base = 2)

    # Plotting and output settings for genes
    fname <- paste(outdir, "/gene_TPM_distribution.png", sep="")
    xlabel <- "log2(TPM + 1)"
    ylabel <- "Count"
    colors <- c("skyblue", "olivedrab3")

    png(filename = fname,
        width = 3000, height = 3000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(gene_TPMs, aes(TPM)) +
               geom_histogram(alpha = 0.5, position="identity", 
                              color = colors[1], fill = colors[1]) +
               xlab(xlabel) + ylab(ylabel) + custom_theme
    print(g)
    dev.off() 
 
    # Plotting and output settings for transcripts
    fname <- paste(outdir, "/transcript_TPM_distribution.png", sep="")
    xlabel <- "log2(TPM + 1)"
    ylabel <- "Count"

    png(filename = fname,
        width = 3000, height = 3000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(transcript_TPMs, aes(TPM)) +
               geom_histogram(alpha = 0.5, position="identity",
                              color = colors[2], fill = colors[2]) +
               xlab(xlabel) + ylab(ylabel) + custom_theme 

    print(g)
    dev.off()

}

plot_isoforms_per_gene <- function(whitelist, custom_theme, outdir) {
    # Plot the number of unique isoforms detected per gene

    isoforms_per_genes <- as.data.frame(table(whitelist$gene_ID))

    # Plotting and output settings
    fname <- paste(outdir, "/isoforms_per_gene.png", sep="")
    xlabel <- "Number of unique isoforms detected per gene"
    ylabel <- "Count"
    col <- "skyblue"

    png(filename = fname,
        width = 3000, height = 2000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(isoforms_per_genes, aes(Freq)) +
               geom_histogram(alpha = 0.8, position="identity", 
                              fill = col, color = col) +
               xlab(xlabel) + ylab(ylabel) + custom_theme
    print(g)

    dev.off()

}



custom_theme <- function() {
    # Create custom GGPlot2 theme
    customTheme <- suppressMessages(theme_bw(base_family = "Helvetica", base_size = 18) +
                       theme(axis.line.x = element_line(color="black", size = 0.5),
                             axis.line.y = element_line(color="black", size = 0.5),
                             axis.text.x = element_text(color="black", size = 16),
                             axis.text.y = element_text(color="black", size = 16)))
        #theme(axis.title.x = element_text(color="black", size=18, vjust = -2, margin=margin(5,0,0,0)),
        #axis.text.x  = element_text(color="black", vjust=0.75, size=16),
        #axis.title.y = element_text(color="black", size=18, margin=margin(0,10,0,0)),
        #axis.text.y  = element_text(color="black", vjust=0.75, size=16))  +
        #theme(legend.text = element_text(color="black", size = 14), legend.title = element_text(color="black", size=11), legend.key.size = unit(0.5, "cm")))

    return(customTheme)
}

load_packages <- function() {

    message("Loading required packages...")

    tryCatch({ suppressPackageStartupMessages(library("DBI"))
               suppressPackageStartupMessages(library("ggplot2"))
               suppressPackageStartupMessages(library("knitr"))
               suppressPackageStartupMessages(library("optparse"))
               suppressPackageStartupMessages(library("RSQLite"))
               suppressPackageStartupMessages(library("readr"))
               suppressPackageStartupMessages(library("rmarkdown"))
               suppressPackageStartupMessages(library("stringr"))
               suppressPackageStartupMessages(library("tinytex"))
 
              }, error = function(cond) {
                   message("One of the required packages is missing.")
                   message("Here's the original error message:")
                   message(cond)
              })

    return
}

parse_options <- function() {

    message("Processing input options...")

    option_list <- list(
        make_option(c("--db"), action = "store", dest = "database",
                    help = "TALON database file"),
        make_option(c("--whitelists", "-w"), action = "store", dest = "whitelists",
                    default = NULL, help = "Comma-delimited list of up to two filesnames. File contains whitelisted transcripts for the datasets provided. If this option is omitted, all genes and transcripts are included in the analyses."),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    help = "Comma-delimited list of up to two dataset names to include in the analysis."),
        #make_option(c("--color"), action = "store", dest = "color_scheme",
        #            default = NULL, help = "blue, red, or green"),
        make_option(c("-o", "--outdir"), action = "store", dest = "outdir",
                    help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
