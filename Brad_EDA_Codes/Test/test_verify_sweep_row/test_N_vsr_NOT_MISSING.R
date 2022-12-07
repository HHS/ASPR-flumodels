test_that(desc = "",
code =
{
  expect_error(object = verify_sweep_row(),regex = "^input_sweep_row is missing\\.$")
})
