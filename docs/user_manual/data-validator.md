<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>

SPDX-License-Identifier: MPL-2.0
-->

# Data Validator

The power-grid-model repository provides an input/update data validator.
For performance reasons, the input/update data is not automatically validated.
However, it would be advisable to validate your data before constructing a PowerGridModel instance.
An alternative approach would be to validate only when an exception is raised, but be aware that not all data errors
will raise exceptions, most of them wil just yield invalid results without warning.
The main validation functions and classes can be included from `power_grid_model.validation`.

```{seealso}
Check the [example](examples/Validation%20Examples.ipynb) for an example of function applications.
```

## Validation Errors

Each validation error is an object which can be converted to a compact human-readable message using `str(error)`.
It contains three member variables `component`, `field` and `ids`, which can be used to gather more specific information
about the validation error, e.g. which object IDs are involved.

```py
class ValidationError:
    # Component(s): e.g. ComponentType.node or [ComponentType.node, ComponentType.line]
    component: ComponentType | list[ComponentType] | None = None

    # Field(s): e.g. AttributeType.id or [AttributeType.line_from, AttributeType.line_to] or
    # [(ComponentType.node, AttributeType.id), (ComponentType.line, AttributeType.id)]
    field: AttributeType | list[AttributeType] | list[tuple[ComponentType, AttributeType]] | None = None

    # IDs: e.g. [1, 2, 3] or [(ComponentType.node, 1), (ComponentType.line, 1)]
    ids: list[int] | list[tuple[ComponentType, int]] | None = None
```

## Validation warnings

Some input is valid, but probably not what was intended.
For such input, the validation functions emit a `ValidationWarning` through the standard Python `warnings` module.
A warning does not change the result: the data is still valid, `validate_input_data` and `validate_batch_data` still
return `None`, and the assertion functions do not raise.
Each kind of warning is emitted at most once per side per call (for batch data it is not repeated per scenario), and
its message names the IDs of the components.

For example, a `ValidationWarning` is emitted when `u1` or `u2` of a `single_phase_transformer` differs by more than
20% from the nominal winding voltage implied by the `u_rated` of its node (see
[Reduced-phase networks](reduced-phase-networks.md#validation)).

Use the `warnings` module to turn warnings into errors, or to record them:

```py
import warnings

from power_grid_model.validation import ValidationWarning, validate_input_data

with warnings.catch_warnings():
    warnings.simplefilter("error", ValidationWarning)  # raise instead of warn
    errors = validate_input_data(input_data, symmetric=False)
```

## Validation functions

### Manual validation

The validation functions below can be used to validate input/batch data manually.
The functions require `input_data` of {py:class}`SingleDataset <power_grid_model.data_types.SingleDataset>` type,
which is power-grid-model input data, and `symmetric: bool`, stating if
the data will be used for symmetric or asymmetric calculations.
`calculation_type: CalculationType` is optional and can be supplied to allow missing values for unused fields; see the
[API reference](../api_reference/python-api-reference.md#enum) for more information.
To validate update/batch data `update_data` of {py:class}`BatchDataset <power_grid_model.data_types.BatchDataset>` type,
power-grid-model update data, should also be supplied.

- `validate_input_data(input_data, calculation_type, symmetric) -> list[ValidationError]` validates input_data.
- `validate_batch_data(input_data, update_data, calculation_type, symmetric) -> dict[int, list[ValidationError]]`
  validates input_data in combination with batch/update data.

```{note}
When doing [batch calculations](./calculations.md#batch-calculations), the input data set by itself is not required to
be valid, as long as all missing and invalid values are overwritten in all scenarios in the batch data set.
This means, that `validate_input_data(input_data, **kwargs)` may fail, while
`validate_batch_data(input_data, update_data, **kwargs)` may succeed.
In such cases, the latter is leading when only running batch calculations.
Running single calculations on an incomplete input data set is, of course, unsupported.
```

Validating a Cartesian product of datasets used in PGM's `update_data` via providing it with `list[BatchDataset]`is
done by validating each individual `BatchDataset` that conforms the Cartesian product.

### Assertions

Instead of manual validation it is possible to use the assertion functions below to assert valid data.
In the case of invalid data a `ValidationException` will be raised, which is an exception storing the name of the
validated data, a list/dict of errors and a convenient conversion to string to display a summary of all the errors when
printing the exception.

- `assert_valid_input_data(input_data, calculation_type, symmetric)` raises `ValidationException`
- `assert_valid_batch_data(input_data, calculation_type, update_data, symmetric)` raises `ValidationException`

Other than the result/exception type, the behaviour of these functions is the same as for the
[manual validation](#manual-validation) functions.

### Utilities

- `errors_to_string(errors, name, details)` converts a set of errors to a human-readable (multi-line) string
  representation

```{seealso}
For more information on the validation functions please see the
[API reference](../api_reference/python-api-reference.md#validation).
```

## Reduced-phase networks

The validation functions check the attributes of one- and two-phase lines (`asym_line` with `phases`),
`single_phase_transformer` and the `connection` of loads and generators.
Validate such networks with `symmetric=False`.
The rules that need the phases of the nodes, which are derived from the topology, are only checked by the calculation:
for example, whether all phases of a transformer winding are present at its node, or whether a load has power on an
absent phase.
See [Reduced-phase networks](reduced-phase-networks.md#validation) for the full list.
