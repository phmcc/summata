withr::defer({
    test_files <- list.files(pattern = "\\.(html|tex|pdf|docx|pptx|rtf)$")
    if (length(test_files) > 0) {
        unlink(test_files)
    }
}, teardown_env())
