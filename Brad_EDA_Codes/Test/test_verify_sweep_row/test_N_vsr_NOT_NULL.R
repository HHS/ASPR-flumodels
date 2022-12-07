test_that(desc = "",
code =
{
  expect_error(object = verify_sweep_row(input_sweep_row = NULL),regex = "^input_sweep_row cannot be NULL.")
})
