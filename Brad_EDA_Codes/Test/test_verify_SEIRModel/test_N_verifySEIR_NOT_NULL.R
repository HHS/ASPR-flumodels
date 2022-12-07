test_that(desc = "",
code =
{
  expect_error(object = verify_SEIRModel(input_sweep_row = NULL),regex = "^input_model cannot be NULL.")
})
