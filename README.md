# General EAF Material and Heat Balance

This repository contains a configurable fixed-form Fortran program for electric
arc furnace (EAF) material and heat balance calculations.

The solver is intended for engineering studies where plant data may be
incomplete or inconsistent. It uses a component-balance matrix and solves
unknown stream masses by weighted least squares, then evaluates sensible heat
and user-defined heat terms.

## Contents

| File | Purpose |
|---|---|
| `eaf_general_balance.f` | Fixed-form Fortran source code. |
| `general_eaf_input.txt` | Commented example input file. |
| `README.md` | This guide. |
| `Makefile` | Optional build/run helper. |
| `.gitignore` | Excludes generated binaries and output files. |

Generated files:

| File | Purpose |
|---|---|
| `eaf_general_balance` or `eaf_general_balance.exe` | Compiled executable. |
| `general_eaf_output.txt` | Solver output report. |

Generated files are ignored by Git.

## Features

The model supports:

- Steel mass and target steel chemistry.
- Scrap, DRI, HBI, hot metal, FeMn, and other charge streams.
- Fluxes and slag-former compositions.
- Oxygen, air, fuel, moisture, carbon injection, and electrodes.
- Slag chemistry and optional slag mass targets.
- Off-gas chemistry and optional measured gas-flow targets.
- Dust and other loss streams.
- User-defined component lists.
- Weighted least-squares material-balance solving.
- Matrix rank, degree-of-freedom, and redundant-equation diagnostics.
- Stream fraction closure and untracked-mass diagnostics.
- Sensible heat from user-supplied coefficients.
- Electrical input, burner/fuel terms, reaction heat, cooling losses, and
  other user-defined heat terms.

The program does not include a thermodynamic database. Heat-content
coefficients and reaction/fuel/loss terms must be supplied in the input file.

## Requirements

GNU Fortran is recommended:

```text
gfortran
```

The source is fixed-form Fortran and should be compiled in legacy mode.

## Quick Start

Compile:

```text
gfortran -std=legacy -ffixed-form -ffixed-line-length-none -O2 \
  eaf_general_balance.f -o eaf_general_balance
```

Run on Linux/macOS:

```text
./eaf_general_balance
```

Run on Windows PowerShell:

```text
.\eaf_general_balance.exe
```

The program reads:

```text
general_eaf_input.txt
```

and writes:

```text
general_eaf_output.txt
```

If `make` is available:

```text
make
make run
```

## Input File Structure

`general_eaf_input.txt` is line-oriented. Blank lines and lines beginning with
`#` or `!` are ignored.

Sections must appear in this order:

1. Case title.
2. Number of conserved components or elements.
3. Component names and balance weights.
4. Number of streams.
5. Stream definitions.
6. Number of soft stream-mass targets or measurements.
7. Soft target definitions.
8. Number of user heat terms.
9. Heat term definitions.

Names must not contain spaces. Use underscores.

## Component Definitions

Each component line has this format:

```text
component_name balance_weight
```

Example:

```text
C    1.0
Mn   1.0
Si   1.0
Fe   1.0
O    1.0
N    1.0
H    1.0
CaO  1.0
```

The listed quantities are usually chemical elements. They can also be
pseudo-components, such as `CaO`, when that is how the plant balance is
specified.

The balance weight controls how strongly each component balance is enforced in
the least-squares objective.

## Stream Definitions

Each stream line has this format:

```text
name side fixed mass hsgn temp_K HA HB HC HD frac_1 frac_2 ...
```

There must be one mass fraction for every component listed in the component
section.

| Field | Meaning |
|---|---|
| `name` | Stream name. |
| `side` | `1` for input, `-1` for output. |
| `fixed` | `1` for known mass, `0` for solved mass. |
| `mass` | Stream mass in kg. Use `0.0` for solved streams. |
| `hsgn` | Sensible-heat sign: `1` demand, `-1` credit, `0` ignored. |
| `temp_K` | Stream temperature in kelvin. |
| `HA HB HC HD` | Sensible-heat coefficients. |
| `frac_1...` | Component mass fractions. |

For example:

```text
SCRAP 1 0 0.0 0.0 298.15 0.0 0.0 0.0 0.0 ...
STEEL -1 1 1000.0 1.0 1873.0 -13456.8 747.6 0.0 0.0 ...
```

## Soft Targets and Measurements

Soft target lines have this format:

```text
name stream_index target_mass_kg weight
```

Example:

```text
GAS_MEAS 12 106.25 0.20
```

This adds a weighted equation:

```text
mass(stream 12) = 106.25 kg
```

A practical interpretation is:

```text
weight = 1 / sigma
```

where `sigma` is the expected uncertainty in kg. A weight of `0.20` corresponds
to about `5 kg` uncertainty.

## Heat Terms

Heat term lines have this format:

```text
name stream_index factor coeff_J_per_kg const_kWh
```

The calculated term is:

```text
Q_k = const_kWh + factor * coeff_J_per_kg * mass(stream_index) / 3.6E6
```

If `stream_index = 0`, the term is a pure constant:

```text
Q_k = const_kWh
```

Sign convention:

```text
positive kWh = heat demand or heat loss
negative kWh = heat supply or exothermic credit
```

Examples:

```text
ELECTRIC_INPUT 0 0.0 0.0 -360.0
PANEL_LOSS     0 0.0 0.0   35.0
FUEL_HEATING_VALUE 8 1.0 -50000000.0 0.0
```

## Governing Equations

### Component Material Balance

For stream `i` and component `j`:

| Symbol | Meaning |
|---|---|
| `m_i` | Stream mass, kg. |
| `x_ij` | Mass fraction of component `j` in stream `i`. |
| `s_i` | Stream side, `+1` for input and `-1` for output. |

The residual for component `j` is:

```text
r_j = sum_i s_i * m_i * x_ij
```

For an exact component balance:

```text
r_j = 0
```

Plant data is often inconsistent, so the program reports residuals rather than
assuming they vanish.

### Matrix Form

Unknown stream masses are collected in vector `x`. The program builds:

```text
A*x = b
```

Rows of `A` come from:

- Component balances.
- Soft stream-mass targets.

### Weighted Least Squares

The program minimizes:

```text
Phi = sum_k (w_k * (A_k*x - b_k))^2
```

The normal equations solved are:

```text
(A^T W^2 A) * x = A^T W^2 * b
```

The code solves these equations by Gaussian elimination with partial pivoting.

### Degrees of Freedom

The program estimates the weighted matrix rank and reports:

```text
DOF = number_of_unknown_stream_masses - rank(A)
```

Interpretation:

- `DOF > 0`: underspecified model.
- `DOF = 0`: independent equations determine the unknown stream masses.
- `number_of_equations - rank(A) > 0`: redundant or extra equations are
  present.

### Stream Closure

For each stream:

```text
FRACSUM_i = sum_j x_ij
UNTRACKED_KG_i = m_i * (1 - FRACSUM_i)
```

If `FRACSUM` differs from one, the component list does not account for all
stream mass. Add an `OTHER` component if closure is required.

### Sensible Heat

For stream `i`, the supplied heat-content function is:

```text
h_i(T) = HA_i + HB_i*T + HC_i*T^2 + HD_i/T
```

The sensible heat contribution is:

```text
Q_sensible,i = hsgn_i * m_i * h_i(T_i) / 3.6E6
```

### Net Heat Balance

The total heat result is:

```text
Q_net = sum_i Q_sensible,i + sum_k Q_k
```

Sign convention:

```text
positive Q_net = additional heat demand
negative Q_net = net heat surplus
```

## Example Case

The included example solves six unknown stream masses:

| Stream | Approximate mass, kg |
|---|---:|
| `SCRAP` | 1013.73771 |
| `FEMN` | 7.52594 |
| `AIR` | 97.35209 |
| `LIME_CAO` | 40.35914 |
| `SLAG` | 63.15415 |
| `GAS` | 107.25349 |

Expected summary:

```text
WEIGHTED MATRIX RANK = 6
DEGREES OF FREEDOM AFTER RANK TEST = 0
REDUNDANT/EXTRA EQUATIONS AFTER RANK TEST = 6
WEIGHTED RMS RESIDUAL = 0.22324 kg
NET HEAT BALANCE = 137.97013 kWh
```

The example is a template, not a universal EAF model. Replace the stream
compositions, measurements, target masses, and heat terms with data for your
case.

## Validation Checklist

After editing the input file, check:

- Every stream has one fraction per component.
- Input streams use `side = 1`.
- Output streams use `side = -1`.
- Known stream masses use `fixed = 1`.
- Solved stream masses use `fixed = 0`.
- Soft targets refer to the correct stream index.
- Heat terms refer to the correct stream index.
- `FRACSUM` and `UNTRACKED_KG` values are understood.
- Rank, DOF, and residuals are physically acceptable.
- Heat-term signs match the convention.

## Limits

The compiled dimensions are:

| Limit | Value |
|---|---:|
| Components/elements | 20 |
| Streams | 60 |
| Unknown stream masses | 40 |
| Equations | 140 |
| Heat terms | 120 |

These are set near the top of `eaf_general_balance.f`:

```fortran
      PARAMETER (MAXE=20,MAXS=60,MAXU=40,MAXQ=140,MAXT=120)
```

Increase them and recompile if needed.

