test_that(desc = "",
code =
{
  error_regex <- "^input_sweep_row must have a column named \\'model\\'\\.$"
  expect_error(object = verify_sweep_row(input_sweep_row = tibble(a = 1, b = 3),regex = error_regex))
  expect_error(object = verify_sweep_row(input_sweep_row = tibble(k = "model"),regex = error_regex))
})
