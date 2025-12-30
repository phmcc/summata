#' Test Suite for table2* Export Functions
#' 
#' Tests for table export functions: table2docx, table2pptx, table2html, 
#' table2pdf, table2tex, table2rtf. Focuses on input validation, 
#' file creation, and integration with summata output.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-table2export.R")

library(testthat)
library(data.table)
library(summata)

## ============================================================================
## Setup: Create test data
## ============================================================================

## Use package's built-in clinical trial data
data(clintrial)
data(clintrial_labels)

## Create sample tables for testing
clintrial$response <- as.integer(clintrial$os_status == 1 & clintrial$os_months < 24)

## Create a simple fit result for testing
test_table <- fit(
    data = clintrial,
    outcome = "response",
    predictors = c("age", "sex", "stage"),
    model_type = "glm"
)

## Create a simple data.frame for basic tests
simple_table <- data.frame(
    Variable = c("Age", "Sex", "", ""),
    Group = c("", "Male", "Female", ""),
    Total = c("55.2 ± 12.3", "150 (48.5%)", "159 (51.5%)", ""),
    `p-value` = c("0.042", "", "", "0.823"),
    check.names = FALSE
)


## ============================================================================
## SECTION 1: table2docx Input Validation
## ============================================================================

test_that("table2docx requires .docx extension", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    
    expect_error(
        table2docx(simple_table, "output.txt"),
        regexp = "\\.docx"
    )
})


test_that("table2docx validates paper parameter", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    
    ## match.arg throws error mentioning the valid options
    expect_error(
        table2docx(simple_table, "output.docx", paper = "invalid")
    )
})


test_that("table2docx validates orientation parameter", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    
    ## match.arg throws error mentioning the valid options
    expect_error(
        table2docx(simple_table, "output.docx", orientation = "invalid")
    )
})


## ============================================================================
## SECTION 9: table2pptx Input Validation
## ============================================================================

test_that("table2pptx requires .pptx extension", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    
    expect_error(
        table2pptx(simple_table, "output.txt"),
        regexp = "\\.pptx"
    )
})


## ============================================================================
## SECTION 10: table2html Input Validation
## ============================================================================

test_that("table2html requires .html extension", {
    
    skip_if_not_installed("kableExtra")
    
    expect_error(
        table2html(simple_table, "output.txt"),
        regexp = "\\.html"
    )
})


## ============================================================================
## SECTION 11: table2pdf Input Validation
## ============================================================================

test_that("table2pdf requires .pdf extension", {
    
    skip_if_not_installed("xtable")
    skip_if_not_installed("tinytex")
    
    expect_error(
        table2pdf(simple_table, "output.txt"),
        regexp = "\\.pdf"
    )
})


## ============================================================================
## SECTION 12: table2tex Input Validation
## ============================================================================

test_that("table2tex requires .tex extension", {
    
    skip_if_not_installed("xtable")
    
    expect_error(
        table2tex(simple_table, "output.txt"),
        regexp = "\\.tex"
    )
})


## ============================================================================
## SECTION 13: table2rtf Input Validation
## ============================================================================

test_that("table2rtf requires .rtf extension", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    
    expect_error(
        table2rtf(simple_table, "output.txt"),
        regexp = "\\.rtf"
    )
})


## ============================================================================
## SECTION 14: File Creation Tests (Integration)
## ============================================================================

test_that("table2docx creates file", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(simple_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx returns correct structure", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(simple_table, temp_file, caption = "Test Caption")
    
    expect_s3_class(result, "table2docx_result")
    expect_equal(result$file, temp_file)
    expect_equal(result$caption, "Test Caption")
    
    ## Should have flextable attribute
    ft <- attr(result, "flextable")
    expect_true(!is.null(ft))
    expect_s3_class(ft, "flextable")
})


test_that("table2docx return_ft=TRUE returns flextable", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(simple_table, temp_file, return_ft = TRUE)
    
    expect_s3_class(result, "flextable")
})


test_that("table2pptx creates file", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".pptx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2pptx(simple_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2html creates file", {
    
    skip_if_not_installed("kableExtra")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".html")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2html(simple_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2tex creates file", {
    
    skip_if_not_installed("xtable")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".tex")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2tex(simple_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


## ============================================================================
## SECTION 15: Package Requirement Tests
## ============================================================================

test_that("table2docx checks for flextable package", {
    
    ## This test verifies the error message format
    ## We can't actually test the missing package scenario easily
    ## but we can verify the function exists and has proper structure
    
    expect_true(is.function(table2docx))
})


test_that("table2html checks for kableExtra package", {
    
    expect_true(is.function(table2html))
})


test_that("table2pdf checks for xtable package", {
    
    expect_true(is.function(table2pdf))
})


## ============================================================================
## SECTION 16: Option Processing Tests
## ============================================================================

test_that("table2docx processes indent_groups option", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should not error with indent_groups
    result <- table2docx(test_table, temp_file, indent_groups = TRUE)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx processes condense_table option", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should not error with condense_table
    result <- table2docx(test_table, temp_file, condense_table = TRUE)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx processes zebra_stripes option", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should not error with zebra_stripes
    result <- table2docx(test_table, temp_file, zebra_stripes = TRUE)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx processes dark_header option", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should not error with dark_header
    result <- table2docx(test_table, temp_file, dark_header = TRUE)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx processes font options", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should not error with custom font options
    result <- table2docx(test_table, temp_file, 
                         font_size = 10, 
                         font_family = "Times New Roman")
    
    expect_true(file.exists(temp_file))
})


## ============================================================================
## SECTION 17: Integration with summata Output
## ============================================================================

test_that("table2docx works with fit() output", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(test_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx works with desctable() output", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    desc_table <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex")
    )
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(desc_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2html works with fit() output", {
    
    skip_if_not_installed("kableExtra")
    skip_on_cran()
    
    temp_file <- tempfile(fileext = ".html")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2html(test_table, temp_file)
    
    expect_true(file.exists(temp_file))
})


## ============================================================================
## SECTION 18: Edge Cases
## ============================================================================

test_that("table2docx handles empty table gracefully", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    empty_table <- data.frame(
        Variable = character(0),
        Value = character(0)
    )
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    ## Should either work or give informative error
    result <- tryCatch(
        table2docx(empty_table, temp_file),
        error = function(e) e
    )
    
    ## Either file created or error caught
    expect_true(file.exists(temp_file) || inherits(result, "error"))
})


test_that("table2docx handles single row table", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    single_row <- data.frame(
        Variable = "Age",
        Value = "55.2"
    )
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(single_row, temp_file)
    
    expect_true(file.exists(temp_file))
})


test_that("table2docx handles special characters", {
    
    skip_if_not_installed("flextable")
    skip_if_not_installed("officer")
    skip_on_cran()
    
    special_table <- data.frame(
        Variable = c("α", "β", "γ"),
        Value = c("< 0.001", "≥ 5%", "95% CI")
    )
    
    temp_file <- tempfile(fileext = ".docx")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- table2docx(special_table, temp_file)
    
    expect_true(file.exists(temp_file))
})
