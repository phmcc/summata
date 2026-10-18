# Add the coded events of a survival model to its analysis rows

The group counts of
[`m2dt()`](https://phmcc.codefloe.page/summata/reference/m2dt.md) sum an
event column by group. The events are coded from the model's `Surv()`
response (see
[`get_event_indicator()`](https://phmcc.codefloe.page/summata/reference/get_event_indicator.md))
into a column of a copy of the model's variables, so that the data
supplied are not modified. Observations whose event `Surv()` codes as
missing, such as a status of 9 among statuses of 0 and 1, did not enter
the model and are not counted.

## Usage

``` r
add_event_column(model, model_class, data)
```

## Arguments

- model:

  Fitted survival model.

- model_class:

  Character string of the model class.

- data:

  Data table of the analysis rows.

## Value

List with `data`, the data with an event column where the events could
be coded, and `event_var`, the name of that column or `NULL`.
