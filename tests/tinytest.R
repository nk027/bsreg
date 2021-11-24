
# if(requireNamespace("tinytest", quietly = TRUE)) {
#   set.seed(42)
#   home <- length(unclass(packageVersion("bsreg"))[[1]]) == 4 # 0.0.0.9000
#   # home <- TRUE
#   if(home) { # Only run locally, let CRAN test examples and the vignette
#     tinytest::test_package("bsreg", at_home = home, pattern = "^.*\\.[rR]$")
#   }
# }
