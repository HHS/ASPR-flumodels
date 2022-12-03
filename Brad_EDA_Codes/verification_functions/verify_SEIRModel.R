verify_SEIRModel <- function(input_model)
{
  assert_that(!missing(input_model),msg = "input_model is missing.")
  assert_that(!is.null(input_model),msg = "input_model cannot be NULL.")
  
  print(input_model)
  print((input_model %>% class))
  assert_that("SEIRModel" %in% (input_model %>% class),msg = "model must be at least of type 'SEIRModel'.")
  
  model_field_names <- names(input_model)
  
  has_correct_model_field_names <- setdiff(model_field_names,c("parameters","rawOutput")) %>% is_empty()
  assert_that(has_correct_model_field_names,msg = "input_model must have and can only have the field names 'parameters' and 'rawOutput.")
  
  return(input_model)
}