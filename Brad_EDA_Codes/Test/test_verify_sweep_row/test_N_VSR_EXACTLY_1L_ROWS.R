test_that(desc = "",
code =
{
  error_regex <- "^input_sweep_row must have exactly 1L rows\\.$"
  expect_error(object = verify_sweep_row(input_sweep_row = tibble(model = c(1,2), b = c(1,3)),regex = error_regex))
  expect_error(object = verify_sweep_row(input_sweep_row = tibble(a =  integer(), model= character()),regex = error_regex))
})