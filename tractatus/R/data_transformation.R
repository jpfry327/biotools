# R/data_transformation.R

#' Convert read counts to TPM (Transcripts Per Million)
#'
#' This function converts raw read counts to TPM values, normalizing for both
#' sequencing depth and gene length. TPM is calculated as:
#' 1. Divide counts by gene length (in kb) to get RPK (Reads Per Kilobase)
#' 2. Divide RPK by the sum of all RPK values and multiply by 1 million
#'
#' @param data A data frame containing gene information and count data
#' @param lengths Column name containing gene lengths (unquoted)
#' @param samples Column selection for count columns (supports tidyselect)
#' @return A data frame with count columns replaced by TPM values
#' @export
#' @import dplyr
#' @import purrr
#' @import rlang
#' @examples
#' \dontrun{
#' # Convert counts to TPM for multiple samples
#' tpm_data <- counts_to_tpm(my_data, gene_length, starts_with("sample_"))
#'
#' # Or specify exact columns
#' tpm_data <- counts_to_tpm(my_data, length, c(sample1, sample2, sample3))
#' }
counts_to_tpm <- function(data, lengths, samples) {
  # Capture the column selections
  lengths_col <- rlang::enquo(lengths)
  samples_cols <- rlang::enquo(samples)

  # Get sample column names
  sample_names <- data %>%
    dplyr::select(!!samples_cols) %>%
    names()

  # Extract lengths and counts
  gene_lengths <- data %>% dplyr::pull(!!lengths_col)

  # Calculate TPM for each sample
  tpm_data <- data %>%
    dplyr::select(!!samples_cols) %>%
    purrr::map_dfc(function(counts) {
      # Step 1: Divide counts by gene length (in kb) -> RPK
      rpk <- counts / (gene_lengths / 1000)
      # Step 2: Divide by sum of RPK and multiply by 1 million -> TPM
      tpm <- rpk / (sum(rpk, na.rm = TRUE) / 1e6)
      return(tpm)
    })

  # Replace original count columns with TPM values
  data %>%
    dplyr::select(-!!samples_cols) %>%
    dplyr::bind_cols(tpm_data)
}
