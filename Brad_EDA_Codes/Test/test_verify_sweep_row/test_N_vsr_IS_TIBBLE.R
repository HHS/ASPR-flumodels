test_that(desc = "",
code =
{
  error_regex <- "input_sweep_row must be a tibble."
  expect_error(object = verify_sweep_row(input_sweep_row = list(model= 1,b = 3),regex = error_regex))
  expect_error(object = verify_sweep_row(input_sweep_row = data.frame(model= 1,b = 3),regex = error_regex))
})