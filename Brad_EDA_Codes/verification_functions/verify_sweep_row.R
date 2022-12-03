verify_sweep_row <- function(input_sweep_row)
{
  assert_that(!missing(input_sweep_row),msg = "input_sweep_row is missing.")
  assert_that(!is.null(input_sweep_row),msg = "input_sweep_row cannot be NULL.")
  
  assert_that(is_tibble(input_sweep_row),msg = "input_sweep_row must be a tibble.")
  assert_that(nrow(input_sweep_row) == 1L,msg = "input_sweep_row must have exactly 1L rows.")
  
  has_column_named_model <- "model" %in% names(input_sweep_row)
  assert_that(has_column_named_model,msg = "input_sweep_row must have a column named 'model'.")
  
  return(input_sweep_row)
}