arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) != 5) {
  stop(
    "Usage: prepare_module16_phenotypes.R <covariates.txt> <bfile.fam> ",
    "<transformed_methylation_adjusted_pcs_prefix> <female_pheno_dir> <male_pheno_dir>",
    call. = FALSE
  )
}

covariates_file <- arguments[1]
fam_file <- arguments[2]
methylation_prefix <- arguments[3]
female_pheno_dir <- arguments[4]
male_pheno_dir <- arguments[5]

autosomal_csv <- paste0(methylation_prefix, ".csv")
female_chrx_csv <- paste0(methylation_prefix, ".Female.chrX.csv")
male_chrx_csv <- paste0(methylation_prefix, ".Male.chrX.csv")
male_chry_csv <- paste0(methylation_prefix, ".Male.chrY.csv")

female_output <- file.path(female_pheno_dir, "methylation_data.csv")
male_output <- file.path(male_pheno_dir, "methylation_data.csv")

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

check_file <- function(path, label) {
  if (!file.exists(path)) {
    fail("Missing ", label, ": ", path)
  }
}

check_file(covariates_file, "covariates file")
check_file(fam_file, "PLINK fam file")
check_file(autosomal_csv, "autosomal methylation CSV")

dir.create(female_pheno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(male_pheno_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(female_output)) {
  unlink(female_output)
}
if (file.exists(male_output)) {
  unlink(male_output)
}

covariates <- read.table(
  covariates_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  colClasses = "character",
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

if (!"Sex_factor" %in% names(covariates)) {
  fail("Cannot find Sex_factor column in ", covariates_file)
}

sample_column <- names(covariates)[1]
covariate_ids <- as.character(covariates[[sample_column]])
sex_factor <- as.character(covariates[["Sex_factor"]])

female_covariate_ids <- covariate_ids[sex_factor == "F"]
male_covariate_ids <- covariate_ids[sex_factor == "M"]

fam <- read.table(
  fam_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  colClasses = "character",
  comment.char = "",
  quote = ""
)

if (ncol(fam) < 2) {
  fail("PLINK fam file has fewer than two columns: ", fam_file)
}

fam_iids <- as.character(fam[[2]])

ordered_fam_ids <- function(covariate_ids_for_sex, sex_label) {
  ordered_ids <- fam_iids[fam_iids %in% covariate_ids_for_sex]
  missing_from_fam <- setdiff(covariate_ids_for_sex, fam_iids)

  message(
    "Sex_factor ", sex_label, " count in covariates: ",
    length(covariate_ids_for_sex)
  )
  message(
    sex_label, " sample count in fam order: ",
    length(ordered_ids)
  )

  if (length(missing_from_fam) > 0) {
    message(
      "WARNING: ", length(missing_from_fam), " ", sex_label,
      " covariate sample(s) were not found in the fam file"
    )
  }

  ordered_ids
}

female_ids <- ordered_fam_ids(female_covariate_ids, "female")
male_ids <- ordered_fam_ids(male_covariate_ids, "male")

do_female <- length(female_ids) > 0
do_male <- length(male_ids) > 0

if (!do_female && !do_male) {
  fail("No female or male samples from covariates were found in ", fam_file)
}

if (do_female && do_male) {
  message("Cohort contains both female and male samples")
} else if (do_female) {
  message("Cohort female only")
} else {
  message("Cohort male only")
}

open_output_connections <- function(outputs, mode) {
  lapply(outputs, function(output) {
    output$con <- file(output$file, open = mode)
    output
  })
}

close_output_connections <- function(outputs) {
  for (output in outputs) {
    if (!is.null(output$con) && isOpen(output$con)) {
      close(output$con)
    }
  }
}

find_columns <- function(header, ids, sex_label, input_file) {
  columns <- match(ids, header)
  missing <- ids[is.na(columns)]

  if (length(missing) > 0) {
    fail(
      "Missing ", sex_label, " sample(s) in methylation CSV header: ",
      input_file, "\nFirst missing sample(s): ",
      paste(head(missing, 20), collapse = ", ")
    )
  }

  columns
}

write_selected_columns <- function(input_file, outputs, write_header) {
  check_file(input_file, "methylation CSV")

  input_con <- file(input_file, open = "rt")
  on.exit(close(input_con), add = TRUE)

  header_line <- readLines(input_con, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    fail("Methylation CSV is empty: ", input_file)
  }

  header <- strsplit(header_line, ",", fixed = TRUE)[[1]]
  if (length(header) < 2) {
    fail("Methylation CSV has fewer than two columns: ", input_file)
  }

  for (i in seq_along(outputs)) {
    outputs[[i]]$columns <- find_columns(
      header,
      outputs[[i]]$ids,
      outputs[[i]]$label,
      input_file
    )
  }

  output_mode <- if (write_header) "wt" else "at"
  outputs <- open_output_connections(outputs, output_mode)
  on.exit(close_output_connections(outputs), add = TRUE)

  if (write_header) {
    for (output in outputs) {
      writeLines(
        paste(c(header[1], output$ids), collapse = ","),
        con = output$con
      )
    }
  }

  line_number <- 1
  repeat {
    lines <- readLines(input_con, n = 10000, warn = FALSE)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      line_number <- line_number + 1
      fields <- strsplit(line, ",", fixed = TRUE)[[1]]

      for (output in outputs) {
        if (length(fields) < max(output$columns)) {
          fail(
            "Line ", line_number, " in ", input_file,
            " has fewer columns than expected"
          )
        }

        writeLines(
          paste(c(fields[1], fields[output$columns]), collapse = ","),
          con = output$con
        )
      }
    }
  }
}

autosomal_outputs <- list()
if (do_female) {
  autosomal_outputs[[length(autosomal_outputs) + 1]] <- list(
    label = "female",
    ids = female_ids,
    file = female_output
  )
}
if (do_male) {
  autosomal_outputs[[length(autosomal_outputs) + 1]] <- list(
    label = "male",
    ids = male_ids,
    file = male_output
  )
}

message("Reading autosomal methylation CSV once: ", autosomal_csv)
write_selected_columns(autosomal_csv, autosomal_outputs, write_header = TRUE)

if (do_female) {
  message("Appending female chrX methylation probes: ", female_chrx_csv)
  write_selected_columns(
    female_chrx_csv,
    list(list(label = "female", ids = female_ids, file = female_output)),
    write_header = FALSE
  )
  message("Wrote female phenotype: ", female_output)
}

if (do_male) {
  message("Appending male chrX methylation probes: ", male_chrx_csv)
  write_selected_columns(
    male_chrx_csv,
    list(list(label = "male", ids = male_ids, file = male_output)),
    write_header = FALSE
  )

  message("Appending male chrY methylation probes: ", male_chry_csv)
  write_selected_columns(
    male_chry_csv,
    list(list(label = "male", ids = male_ids, file = male_output)),
    write_header = FALSE
  )
  message("Wrote male phenotype: ", male_output)
}

message("Module 16 phenotype preparation completed")
