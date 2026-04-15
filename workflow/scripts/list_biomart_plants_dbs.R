library(biomaRt)

output <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]

con <- file(log_file, open = "wt")
sink(con, type = "message")

sink(output)

cat("=== Plants host (plants.ensembl.org) ===\n")
tryCatch(
    print(listMarts(host = "https://plants.ensembl.org")),
    error = function(e) cat("ERROR:", conditionMessage(e), "\n")
)

cat("\n=== Default host (ensembl.org) ===\n")
tryCatch(
    print(listMarts()),
    error = function(e) cat("ERROR:", conditionMessage(e), "\n")
)

sink()
sink(type = "message")
close(con)
