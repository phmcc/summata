# Package Imports

Declares the imports that are called without a namespace prefix, so that
R CMD check does not report them as undefined globals. Functions called
as `package::function()` elsewhere in the package need no entry here,
but their package must still appear under Imports in DESCRIPTION.
