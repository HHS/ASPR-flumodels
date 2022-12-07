test_that(desc = "",
code =
{
  expect_error(object = verify_SEIRModel(),regex = "^input_model is missing\\.$")
})
