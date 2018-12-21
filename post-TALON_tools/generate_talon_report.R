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
 
    #print(count_known_genes(database, datasets))
    #print(count_known_transcripts(database, datasets))

    plot_read_length_distribution(database, datasets[1], custom_theme, outdir)

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

plot_read_length_distribution <- function(database, datasets, custom_theme, outdir) {
    # Use the read lengths recorded in the 'observed' table to plot a read length
    # distribution. In the query, grab the annotation status so that I can use 
    # it in future versions if I want 

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Query to fetch read lengths of observed transcripts, along with the 
    # annotation status
    dataset_str <- paste("('", paste(datasets, collapse = "','"), "')", sep="")
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
                               AND ta.value = 'KNOWN'",
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

custom_theme <- function() {
    # Create custom GGPlot2 theme
    customTheme <- suppressMessages(theme_bw(base_family = "Helvetica", base_size = 18) +
                       theme(axis.line.x = element_line(color="black", size = 0.5),
                             axis.line.y = element_line(color="black", size = 0.5)))
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
                    default = NULL, help = "TALON database file"),
        #make_option(c("--whitelist"), action = "store", dest = "whitelist",
        #            default = NULL, help = "File of whitelisted transcripts for the Pacbio data"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
        #make_option(c("--ik"), action = "store", dest = "illumina_kallisto",
        #            default = NULL, help = "Illumina Kallisto file."),
        #make_option(c("--color"), action = "store", dest = "color_scheme",
        #            default = NULL, help = "blue, red, or green"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
