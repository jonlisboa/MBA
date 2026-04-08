# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an Agent-Based Model (ABM) simulating a spatial economy of Ilhéus, Brazil. It models households, firms, labor markets, government policy (Bolsa Família), and macroeconomic dynamics using pure R.

## Running the Model

**Requirements:** R with `sf` and `dplyr` packages installed.

```r
install.packages(c("sf", "dplyr"))
```

**Run:**
```r
source("MODEL V5")
```

No build step — the file runs end-to-end as a simulation script.

**Important:** Line ~173 hardcodes a Windows path for the geobr CSV:
```r
csv_path <- "C:\\Users\\Usuário\\OneDrive\\Stock Flow Consistent\\br_geobr_mapas_municipio.csv"
```
Update this path when running on a different machine.

**Outputs** (written to working directory, prefixed with a random ID):
- `{id}_abm_v5_tfp_cyclic_macro_pop.csv` — time-series macro panel
- `{id}_abm_v5_tfp_cyclic_hh_final_pop.csv` — final household cross-section
- `{id}_abm_v5_tfp_cyclic_overview.png` — 15-panel diagnostic plot

## Key Parameters (lines 18–99)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `N_HH` | 2000 | Initial number of households |
| `MAXT` | 300 | Simulation periods |
| `BURN_IN` | 90 | Periods excluded from post-sim analysis |
| `POP_GROWTH_RATE` | 0.002 | Births per period |
| `STARTUP_CAPITAL` | 75 | Initial firm capital for entrepreneurs |

## Architecture

The file (`MODEL V5`, ~1,351 lines) is structured as five sequential sections:

### 1. Geography Setup (lines 161–336)
Loads municipality boundary from geobr CSV, converts to UTM, defines 25+ named road networks with waypoints, and builds a combined spatial density kernel (HDI + road proximity) used to place agents realistically.

### 2. Agent Initialization (lines 344–449)
Creates household and firm arrays, computes a full N_HH × N_HH Euclidean distance matrix, and seeds 2% of population as initial firms.

### 3. Main Simulation Loop (lines 480–1205)
Each of `MAXT` periods executes these phases in order:

| Phase | Lines | Description |
|-------|-------|-------------|
| A. Closure | 494–511 | Insolvent firms exit; owner covers deficit |
| A2. Population | 534–628 | Births added; inherit wealth |
| B. Entry | 630–719 | Logit-based firm entry by entrepreneurs |
| C. Planning | 742–771 | Firms set prices, plan production and hiring |
| D. Labour market | 773–816 | Spatial job search and matching |
| E. Layoffs | 822–836 | Over-employed firms shed workers |
| F. Production | 846+ | Output, wages, taxes |
| G–K. Consumption/Finance/Gov | — | Spending, firm P&L, government fiscal accounts |
| L. Record | 1117–1203 | Log 48-column macro/distributional statistics |

### 4. Output & Diagnostics (lines 1207–1351)
Post-simulation TFP trend analysis, summary statistics printed to console, CSV exports, and PNG generation.

## Key Data Structures

**Household arrays** (length grows with births):
- `hh_M` — wealth; `hh_L` — debt; `hh_Employed` — employer firm index (0 = unemployed)
- `hh_Wage`, `hh_firm_owner`, `hh_unemp_spells`, `hh_emp_spells`

**Firm arrays** (length = N_SLOTS, grows with births):
- `f_status` — `"Open"` or `"Vacant"`
- `f_M`, `f_L`, `f_L_max` — money, loans, credit limit
- `f_N`, `f_N_des` — actual and desired employment
- `f_Y`, `f_Inv`, `f_p`, `f_w`, `f_a` — output, inventory, price, wage, TFP
- `f_Profits`, `f_Sold` — period accounting

**Government** (`lg` list): `lg$M`, `lg$L`, `lg$tax_rev`, `lg$transfer`, `lg$deficit`

**History** (`hist` dataframe): 48 columns, one row per period — unemployment, wages, price index (Laspeyres), Gini, firm counts, government accounts, etc.

## Design Notes

- All indices are 1-based (R convention).
- Spatial matching uses exponential decay on a precomputed Euclidean distance matrix — no graph library.
- Downward nominal wage stickiness is enforced explicitly.
- The model uses EWMA for firms' sales expectations and softmax for stochastic job matching.
- Population arrays (households and firms) grow dynamically each period via births.
