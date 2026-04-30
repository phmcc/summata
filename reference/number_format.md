# Number Formatting Utilities

Internal utilities for locale-aware number formatting across all summata
output functions. Supports preset locales (US, European, SI/ISO, plain)
and fully custom separator definitions.

## Global Option

The default number format can be set once per session:

      options(summata.number_format = "eu")

This avoids passing `number_format` to every function call.
