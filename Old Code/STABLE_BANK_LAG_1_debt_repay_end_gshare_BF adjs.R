# =============================================================================
# ABM v5 stable2 --� R port
# Ilhéus spatial economy: logistic entry, owner-liability closure,
# Bolsa Família transfers, pseudo-capital (a as depreciating asset),
# Laspeyres price index, real boundary from geobr CSV.
#
# Requires: sf, dplyr  (for boundary parsing)
# All spatial search done via Euclidean distance matrix (no dodgr needed).
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

#setwd("C:\\Users\\r3420266\\OneDrive\\TCC MBA")
set.seed(rnorm(1,30,3))

# -�--�- Parameters -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
N_HH          <- 1500
MAXT          <- 300

# -�--�- Population growth -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
# Demographic engine: each period, ceil(N_HH * POP_GROWTH_RATE) new households
# are born inside the municipality boundary, inheriting a small wealth endowment
# drawn from the lower tail of the current HH wealth distribution (mimicking
# entry with modest initial assets, consistent with lifecycle models;
# see Huggett 1996, Aiyagari 1994).
# New entrants are unemployed and have zero firm ownership.
# N_HH is updated dynamically; all HH arrays grow accordingly.
# Set POP_GROWTH_RATE <- 0 to disable without changing anything else.
GOV_TRANSF    <- 5
# ----- Bolsa Família inflation correction rule -----------------------------
# GOV_TRANSF is fixed in nominal terms each period (no wage-tracking).
# Every BF_INFL_FREQ periods, if cumulative inflation since the last
# correction exceeds zero, the transfer is partially uprated to restore
# real purchasing power. The nominal transfer NEVER falls (downward
# rigidity): deflation leaves GOV_TRANSF unchanged, improving real value.
# GOV_TRANSF_BASE is a hard floor — the transfer never drops below it.
#
# Glossary:
#   GOV_TRANSF_BASE  : nominal floor (never breached downward)
#   BF_INFL_PASSTHRU : θ, fraction of cumulative inflation passed through (0–1)
#   BF_INFL_FREQ     : κ, periods between corrections (12 = annual)
#   bf_last_pi       : price index at last correction (state variable)
GOV_TRANSF_BASE  <- 5       # hard nominal floor
BF_INFL_PASSTHRU <- 0.50    # 50% pass-through of cumulative inflation
BF_INFL_FREQ     <- 12L     # correct every 12 periods
POP_GROWTH_RATE  <- 0.00        # ~0.1% per period
# SFC birth endowment: each parent transfers HERITAGE_FRAC of their wealth.
# Total endowment = 2 * HERITAGE_FRAC * mean(parental wealth).
# Government tops up to POP_WEALTH_FLOOR if combined transfer is too small.
HERITAGE_FRAC    <- 0.10         # 10% of each parent's wealth → 20% joint
POP_WEALTH_FLOOR <- GOV_TRANSF   # minimum endowment (government tops up if needed)
MIN_WAGE      <- 2
AS_SCALE       <- 1.0   # animal spirits scale: floor = AS_SCALE * mean(wage)
ANIMAL_SPIRITS_MIN <- MIN_WAGE  # hard floor of last resort if wages collapse
THETA         <- 0.6
ALPHA_S       <- 0.35   # was 0.6 --� slower sales expectation update damps demand-signal amplification
LAMBDA_N      <- 0.50   # was 0.95 --� gradual hiring spreads entry waves over ~4 periods
# instead of snapping to target immediately; key fix for oscillations
MU_PHI        <- 0.02
MU_MIN        <- 0.01
MU_MAX        <- 1000
PHI_W         <- 0.15   # was 0.25 --� slower wage adjustment reduces wage-price feedback

MAX_WAGE      <- 1000
LAMBDA0       <- 1e-3
P_FLOOR       <- 0.1
GAMMA_MD      <- 0.03     # was 0.05 --� slower markdown, reduces deflationary spikes
# Consumption function: hh_des_C = A1*hM + A2*wage
# Equilibrium nominal hM* = (1-A2)/A1 × wage × N_HH
# Real hM grows over time as LBD raises productivity and prices deflate.
# A1=0.15, A2=0.85 produces a realistic precautionary savings buffer
# and keeps demand from over-reacting to wealth fluctuations.
A1            <- 0.1 #from 0.04
A2            <- 0.85
RL_FIRM       <- 0.01
# ----- Banking sector (endogenous money, Godley & Lavoie 2007, ch.7-11) --
# Single bank, endogenous money creation. Credit creates deposits.
# SFC mechanics:
#   Lending:    f_M ↑, f_L ↑           (money CREATED — both sides expand)
#   Repayment:  f_M ↓, f_L ↓           (money DESTROYED — both sides shrink)
#   Loan int:   f_M ↓ (via profit chain), bank_NW ↑  (spread income)
#   Dep int:    hh_M/f_M/lg$M ↑, bank_NW ↓           (interest expense)
#   Bad debt:   bank_NW ↓              (equity absorbs loss, no money vanishes)
#   Dividends:  hh_M ↑, bank_NW ↓     (every period, from net interest income)
# CAR constraint: total loans ≤ bank_NW / CAR_MIN  (Basel-like)
# Total M changes via: credit creation/destruction + government deficits.
# References: Lavoie (2014) ch.4; Moore (1988); McLeay et al. (2014) BoE QB;
#             Minsky (1986) ch.9; Nikolaidi & Stockhammer (2017).
RD               <- 0.001    # deposit interest rate (per period)
CAR_MIN          <- 0.08     # minimum capital adequacy ratio (Basel-like)
LOAN_MATURITY    <- 50L      # amortisation horizon (periods to full repayment)
STARTUP_LEVERAGE <- 1.0      # bank matches owner equity 1:1 at entry
BANK_NW_INIT_MULT <- 20    # initial bank_NW = mult * N_SEED * STARTUP_CAPITAL
BANK_DIV_RATE    <- 0.50     # fraction of period net interest income distributed
# Minimum bank equity before dividends are suspended.
# Below this threshold the bank retains all profit to rebuild capital,
# mimicking Basel III Pillar 2 supervisory dividend restrictions.
# At BANK_NW_FLOOR=0, dividends are always paid (original behaviour).
# BANK_NW_FLOOR is a clean Monte Carlo treatment variable: higher values
# make the bank more conservative, build equity faster, and widen the
# CAR lending headroom — at the cost of lower household dividend income.
# Rule of thumb: set to ~10% of mean outstanding firm loans.
BANK_NW_FLOOR    <- 2500     # suspend dividends if bank_NW falls below this
PAYROLL_TAX   <- 0.10
PROFIT_TAX    <- 0.27
WEALTH_FLOOR  <- 10        # overridden dynamically each period
INCOME_FRAC   <- 0.8
RL_GOV        <- 0.005
BF_MEDIAN_FRAC  <- 0.45    # transfer = 45% of median employed wage
BF_WEALTH_FRAC  <- 4.00    # wealth floor = 4x median wage
# ----- Fiscal rule toggle -----------------------------------------------
# FISCAL_STABILIZER = TRUE  : Keynesian stabiliser mode.
#   G = G_SHARE * EWMA(tax_rev), smoothed over ~FISCAL_EWMA_TAU periods.
#   Government spends against trend; deficits expand in downturns.
#   The smoothed revenue base (tax_rev_smooth) is a state variable
#   initialised to 0 and updated every period regardless of the toggle.
#
# FISCAL_STABILIZER = FALSE : Pro-cyclical mode (original behaviour).
#   G = G_SHARE * max(tax_rev - transfer_bill, 0).
#   Procurement tracks current net revenue; collapses in recessions.
#
# Both modes: funded from lg$M first; deficit-financed via lg$L if short.
# G_SHARE and FISCAL_EWMA_TAU are shared across both modes.
# Toggle is a clean treatment variable for Monte Carlo factorial design.
FISCAL_STABILIZER <- TRUE
FISCAL_EWMA_TAU   <- 20    # smoothing horizon (periods); alpha = 1/tau
# ----- Endogenous G_SHARE (fiscal buffer targeting) ------------------------
# G_SHARE adjusts each period toward a desired lg$M cash buffer.
# When lg$M > LG_M_TARGET the government has excess cash → G_SHARE rises
# → more procurement → cash drawn down. When lg$M < LG_M_TARGET the
# government is cash-constrained → G_SHARE falls → less borrowing needed.
# This replaces the fixed G_SHARE and resolves the cash-hoarding failure
# mode (G30) and the perpetual-deficit failure mode (G80) simultaneously.
#
# Mechanism (stock-flow norm, Godley & Lavoie 2007 ch.3):
#   adj_t = G_SHARE_PHI * (lg$M_t - LG_M_TARGET) / (tax_rev_smooth_t + 1)
#   G_SHARE_t = clamp(G_SHARE_{t-1} + adj_t, G_SHARE_MIN, G_SHARE_MAX)
#
# Glossary:
#   G_SHARE_INIT     : starting value; model converges away from it quickly
#   G_SHARE_MIN/MAX  : bounds preventing degenerate fiscal stances
#   G_SHARE_PHI      : adjustment speed (0.05-0.15; higher = faster reversion)
#   M_LG_TARGET_MULT : desired buffer as multiple of one period's transfer bill
G_SHARE_INIT      <- 0.60   # starting value
G_SHARE_MIN       <- 0.4   # floor: always some procurement
G_SHARE_MAX       <- 0.80   # ceiling: bounded by revenue capacity
G_SHARE_PHI       <- 0.15   # adjustment speed
M_LG_TARGET_MULT  <- 1    # target buffer = 2x one period's transfer bill
# ----- Sovereign debt repayment rule ---------------------------------------
# When lg$M exceeds LG_DEBT_FLOOR, the government retires a fraction
# LG_AMORT_RATE of outstanding lg$L each period (smooth amortisation).
# SFC mechanics: lg$M ↓, lg$L ↓ by equal amount → money destroyed,
# bank_R and bank_NW unchanged (symmetric balance sheet contraction).
# This prevents lg$L from accumulating without bound under chronic deficits
# while preserving countercyclical capacity: in downturns lg$M falls below
# LG_DEBT_FLOOR, repayment halts, and the stabiliser operates freely.
# LG_DEBT_FLOOR: minimum cash buffer (one period's transfer bill at ~20%
#   qualifying rate; keeps the government solvent through spending peaks).
# LG_AMORT_RATE: fraction of lg$L retired per period when funded (2-5%).
#   At 0.02, a debt stock of 300k shrinks by ~6k/period when surpluses allow.
#   Interpretable as a fiscal sustainability parameter for Monte Carlo design.
LG_DEBT_FLOOR  <- GOV_TRANSF * N_HH * 0.25   # ~25% of max transfer bill
LG_AMORT_RATE  <- 0.005                         # 2% of outstanding debt/period

STARTUP_CAPITAL   <- 15
SC_WEALTH_FRAC    <- 0.25   # was 0.50 — entry barrier at 25% of mean hH wealth
# Hard floor for STARTUP_CAPITAL: SC >= SC_FLOOR_TRANSF_MULT * GOV_TRANSF.
# Prevents SC from collapsing in downturns and triggering entry explosions.
# At GOV_TRANSF≈5 and mult=8: floor=40. Raise mult to tighten entry.
SC_FLOOR_TRANSF_MULT <- 6 #from 8
INTERCEPT_OPEN    <- -13 #CHANGED FROM -6.25
BETA_M            <-  0.1 #from 0.03
BETA_UNEMP        <-  1.8
BETA_GAIN         <-  0.10
BETA_USPELLS      <-  0.15
BETA_ESPELLS      <- -0.07
USPELLS_CAP       <- 20L    # logit saturation: spells beyond 20 add no further push
ENTRY_RATE_CAP    <- 0.10   # max new firms per period as fraction of N_HH
ENTRY_COOLDOWN    <- 1L
# ----- Firm exit rule ---------------------------------------------------
# 1. Cash illiquidity: f_M - f_L ≤ 0
#    Cash minus loan obligations ≤ 0. Inventory excluded (illiquid, price-distorted).
# 2. Sustained loss: f_Profit_smooth < -LOSS_MULT * f_WB (no f_L > 0 guard).
# 3. Grace period: firms younger than CLOSURE_GRACE are exempt from both
#    conditions. With the production-sales lag, a new firm needs at least
#    2 periods before it has any revenue (period 1: produce into WIP;
#    period 2: WIP matures, first sales). CLOSURE_GRACE = 3 gives one
#    additional buffer period. Without the grace period, every new firm
#    is mechanically insolvent from birth (f_M - f_L ≤ 0 after paying
#    wages with zero revenue) and is killed before WIP ever matures.
#    See Hall & Lerner (2010) on startup runway periods.
CLOSURE_GRACE     <- 5L      # periods before closure rule applies
ALPHA_PROF        <- 0.25   # EWMA weight (~4-period memory)
LOSS_MULT         <- 0.6    # close if smoothed loss > 60% of wage bill (from 20)
BETA_SAT          <- -0.06 #CHANGED FROM -0.05
SIGMA_W           <-  0.10
SAT_RADIUS        <-  10000
OWNER_PROFIT_SHARE <- 0.40
M_DIST_RATE       <- 0.2
M_DIST_PERIOD     <- 12L
getwd()
# ----- TFP: passive learning-by-doing (Arrow 1962) ----------------------
# Each open firm gains LBD_RATE per period (idiosyncratic noise LBD_SIGMA).
# New entrants always start at f_a = 1.0 (frontier technology).
# LBD_ENTRY_FRAC removed: it creates attractor M* = n*LBD/(k*(1-FRAC))
# independent of starting conditions, pulling avg_a well below 1.0.
LBD_RATE        <- 0.0002
LBD_SIGMA       <- 0.003
A_MIN           <- 0.50   # floor raised from 0.1 to prevent extreme price spikes
# Excess cash dividend: firm distributes fM above FM_CAP_WAGE_PERIODS
# weeks of its own wage bill back to its owner every period.
FM_CAP_WAGE_PERIODS <- 12 #from 10

# -�--�- Wage--�productivity linkage (heterodox, non-neoclassical) -�--�--�-

wageperiod <- 1
priceperiod <- 1

# Mechanism 1: Efficiency wage floor (Akerlof-Yellen fair wage)
#   Firms with above-average TFP defend a higher wage floor.
#   f_w >= MIN_WAGE * (f_a / mean_a)^ETA_W
#   ETA_W = 0 -> no linkage; ETA_W = 1 -> full proportional indexation.
#   Empirical target: 0.4-0.6 (Carlsson, Messina & Skans 2016).
#   Set USE_EFF_WAGE <- FALSE to disable without changing parameters.
ETA_W        <- 0.6
USE_EFF_WAGE <- FALSE   # OFF for now; set TRUE to activate


# Mechanism 2: Kaleckian surplus sharing (Bhaduri & Marglin 1990)
#   Workers capture THETA_W share of each period's TFP *gain* as a wage rise.
#   Only positive TFP gains are passed through (downward nominal rigidity).
#   THETA_W = 0 -> no linkage; THETA_W = 1 -> full pass-through of gains.
#   Empirical target: 0.2-0.4 (labour share dynamics literature).
THETA_W      <- 0.5 #from 0.13
USE_SURPLUS  <- TRUE    # ON

QUIT_FRICTION  <- 10.0
QUIT_SIGMA     <- 1.5
SWITCH_SIGMA   <- 0.5
ACCEPT_SIGMA   <- 0.1
COST_PER_M     <- 1e-3 # CHANGED FROM 1E-4
MAX_CONSIDERED <- 30L
MAX_VISITS     <- 100L

# -�--�- Utilities -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
sigmoid   <- function(x) 1 / (1 + exp(-x))
clamp     <- function(x, lo, hi) pmin(pmax(x, lo), hi)
ewma_fn   <- function(p, r, a) (1-a)*p + a*r

softmax_sample <- function(scores, temp = 1.0) {
  z <- scores / temp - max(scores / temp)
  w <- exp(z)
  sample.int(length(w), 1L, prob = w / sum(w))
}

gini_fn <- function(arr) {
  a <- sort(pmax(arr, 0)); n <- length(a)
  if (n == 0 || sum(a) == 0) return(0)
  (2 * sum(seq_len(n) * a) - (n + 1) * sum(a)) / (n * sum(a))
}

spatial_search <- function(d_row, cand, lam0) {
  if (length(cand) == 0) return(integer(0))
  d   <- ifelse(is.finite(d_row[cand]), d_row[cand], Inf)
  lam <- lam0
  found <- integer(0)
  while (length(found) == 0 && lam >= 1e-9) {
    hit   <- rbinom(length(cand), 1L, exp(-lam * d)) == 1L
    found <- cand[hit]
    lam   <- lam * 0.5
  }
  found
}

# -�--�- Geography: Ilhéus boundary from geobr CSV -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
LON0    <- -39.219802
LAT_C   <- -14.761865
SCALE_X <- 111320 * cos(LAT_C * pi / 180)
SCALE_Y <- 110540

ll_to_utm <- function(lon, lat) {
  list(x = (lon - LON0) * SCALE_X + 50000,
       y = (lat - LAT_C) * SCALE_Y + 50000)
}

# Parse MULTIPOLYGON WKT from geobr CSV
csv_path <- "C:\\Users\\botsi\\OneDrive\\Stock Flow Consistent\\br_geobr_mapas_municipio.csv"
#csv_path <- "C:\\Users\\R3420266\\OneDrive\\Stock Flow Consistent\\br_geobr_mapas_municipio.csv"
# if (!file.exists(csv_path)) {
#   # Try project path
#   csv_path <- "/mnt/project/br_geobr_mapas_municipio.csv"
# }
#stopifnot(file.exists(csv_path))

raw <- readLines(csv_path, warn = FALSE)
row <- raw[grep("2913606", raw)[1]]
wkt_start <- regexpr("MULTIPOLYGON", row)
wkt <- substr(row, wkt_start, nchar(row))
wkt <- gsub('"+$', '', wkt)

# Extract all ring coordinate strings
ring_matches <- gregexpr(
  "\\((-?[0-9]+\\.[0-9]+ -?[0-9]+\\.[0-9]+(?:, -?[0-9]+\\.[0-9]+ -?[0-9]+\\.[0-9]+)*)\\)",
  wkt, perl = TRUE)
rings_raw <- regmatches(wkt, ring_matches)[[1]]
rings_raw <- gsub("^\\(|\\)$", "", rings_raw)

parse_ring <- function(s) {
  pts <- strsplit(trimws(s), ",")[[1]]
  do.call(rbind, lapply(pts, function(p) as.numeric(strsplit(trimws(p), " ")[[1]])))
}
polys <- lapply(rings_raw, parse_ring)

# Main polygon = ring with most vertices (ring index 6 in 1-based R)
sizes   <- sapply(polys, nrow)
main_poly <- polys[[which.max(sizes)]]
colnames(main_poly) <- c("lon", "lat")

utm_main <- ll_to_utm(main_poly[,"lon"], main_poly[,"lat"])
MUNI_X   <- utm_main$x
MUNI_Y   <- utm_main$y
XMIN <- min(MUNI_X); XMAX <- max(MUNI_X)
YMIN <- min(MUNI_Y); YMAX <- max(MUNI_Y)
cat(sprintf("Boundary loaded: %d vertices, %.0f x %.0f km\n",
            nrow(main_poly),
            (XMAX-XMIN)/1000, (YMAX-YMIN)/1000))

# Point-in-polygon via sf
muni_sf <- st_sf(geometry = st_sfc(
  st_polygon(list(cbind(MUNI_X, MUNI_Y))), crs = NA_crs_))

point_in_muni <- function(x, y) {
  pts <- st_sf(geometry = st_sfc(lapply(seq_along(x),
                                        function(i) st_point(c(x[i], y[i]))), crs = NA_crs_))
  as.logical(st_within(pts, muni_sf, sparse = FALSE)[, 1])
}

# -�--�- Road network definition -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
roads_ll <- list(
  BR415   = rbind(c(-39.49,-14.43),c(-39.44,-14.52),c(-39.38,-14.60),
                  c(-39.33,-14.65),c(-39.30,-14.68),c(-39.22,-14.78),
                  c(-39.15,-14.87),c(-39.18,-14.96),c(-39.22,-15.05)),
  BA001   = rbind(c(-39.05,-14.43),c(-39.02,-14.52),c(-39.00,-14.60),
                  c(-38.99,-14.68),c(-38.98,-14.78),c(-38.99,-14.88),
                  c(-39.01,-14.97),c(-39.04,-15.07)),
  RAD1    = rbind(c(-39.047,-14.789),c(-39.10,-14.79),c(-39.18,-14.79),
                  c(-39.27,-14.79),c(-39.33,-14.79)),
  RAD2    = rbind(c(-39.047,-14.789),c(-39.07,-14.72),c(-39.10,-14.63),
                  c(-39.15,-14.54),c(-39.21,-14.46),c(-39.28,-14.41)),
  RAD3    = rbind(c(-39.047,-14.789),c(-39.06,-14.86),c(-39.10,-14.93),
                  c(-39.17,-15.00),c(-39.22,-15.06)),
  NW1     = rbind(c(-39.49,-14.43),c(-39.49,-14.37),c(-39.44,-14.37),
                  c(-39.38,-14.38),c(-39.30,-14.40)),
  NW2     = rbind(c(-39.44,-14.52),c(-39.44,-14.47),c(-39.44,-14.42)),
  SEC1    = rbind(c(-39.38,-14.60),c(-39.27,-14.60),c(-39.18,-14.60),
                  c(-39.10,-14.61),c(-39.02,-14.61)),
  SEC2    = rbind(c(-39.22,-15.05),c(-39.13,-15.04),c(-39.06,-15.05)),
  RUR1    = rbind(c(-39.44,-14.52),c(-39.44,-14.62),c(-39.44,-14.72)),
  RUR2    = rbind(c(-39.33,-14.65),c(-39.38,-14.72),c(-39.42,-14.79),
                  c(-39.40,-14.87)),
  CAP_NW1 = rbind(c(-39.49,-14.43),c(-39.54,-14.44),c(-39.54,-14.50)),
  CAP_NW2 = rbind(c(-39.44,-14.52),c(-39.50,-14.55),c(-39.54,-14.58)),
  CAP_FW  = rbind(c(-39.40,-14.87),c(-39.44,-14.82),c(-39.48,-14.76),c(-39.50,-14.70)),
  CAP_CN1 = rbind(c(-39.21,-14.46),c(-39.27,-14.46),c(-39.34,-14.46),c(-39.40,-14.46)),
  CAP_CN2 = rbind(c(-39.28,-14.41),c(-39.28,-14.50),c(-39.28,-14.58),c(-39.28,-14.65)),
  CAP_CE1 = rbind(c(-39.18,-14.60),c(-39.14,-14.67),c(-39.10,-14.74),c(-39.07,-14.79)),
  CAP_CE2 = rbind(c(-39.10,-14.61),c(-39.06,-14.68),c(-39.04,-14.72)),
  CAP_NC1 = rbind(c(-39.02,-14.52),c(-39.01,-14.45),c(-39.02,-14.38)),
  CAP_NC2 = rbind(c(-39.00,-14.60),c(-38.99,-14.52),c(-39.00,-14.43)),
  CAP_SW1 = rbind(c(-39.17,-15.00),c(-39.24,-15.00),c(-39.31,-14.98),c(-39.38,-14.96)),
  CAP_SW2 = rbind(c(-39.13,-15.04),c(-39.19,-15.06),c(-39.26,-15.07)),
  CAP_CW  = rbind(c(-39.30,-14.68),c(-39.35,-14.74),c(-39.39,-14.80),c(-39.40,-14.87)),
  CAP_NI1 = rbind(c(-39.33,-14.46),c(-39.33,-14.56),c(-39.33,-14.66),c(-39.33,-14.76)),
  CAP_NI2 = rbind(c(-39.22,-14.50),c(-39.22,-14.59),c(-39.22,-14.68))
)

road_sigma <- c(BR415=1200,BA001=1200,RAD1=800,RAD2=800,RAD3=800,
                NW1=700,NW2=700,SEC1=600,SEC2=600,RUR1=400,RUR2=400,
                CAP_NW1=350,CAP_NW2=350,CAP_FW=350,CAP_CN1=350,CAP_CN2=350,
                CAP_CE1=350,CAP_CE2=350,CAP_NC1=350,CAP_NC2=350,
                CAP_SW1=350,CAP_SW2=350,CAP_CW=350,CAP_NI1=350,CAP_NI2=350)
road_weight<- c(BR415=3.0,BA001=3.5,RAD1=2.5,RAD2=2.5,RAD3=2.0,
                NW1=2.0,NW2=1.8,SEC1=1.5,SEC2=1.2,RUR1=0.8,RUR2=0.8,
                CAP_NW1=0.55,CAP_NW2=0.55,CAP_FW=0.55,CAP_CN1=0.55,CAP_CN2=0.55,
                CAP_CE1=0.55,CAP_CE2=0.55,CAP_NC1=0.55,CAP_NC2=0.55,
                CAP_SW1=0.55,CAP_SW2=0.55,CAP_CW=0.55,CAP_NI1=0.55,CAP_NI2=0.55)

# Build dense point arrays (1 pt / 200m), clipped to municipality
build_road_pts <- function(waypoints_ll, sigma = 200) {
  u <- ll_to_utm(waypoints_ll[,1], waypoints_ll[,2])
  rx <- u$x; ry <- u$y
  pts <- list()
  for (i in seq_len(nrow(waypoints_ll) - 1)) {
    d <- sqrt((rx[i+1]-rx[i])^2 + (ry[i+1]-ry[i])^2)
    n <- max(2L, as.integer(d / sigma))
    pts[[i]] <- cbind(
      seq(rx[i], rx[i+1], length.out = n),
      seq(ry[i], ry[i+1], length.out = n))
  }
  all_pts <- do.call(rbind, pts)
  inside  <- point_in_muni(all_pts[,1], all_pts[,2])
  all_pts[inside, , drop = FALSE]
}

cat("Building road network...\n")
road_pts <- lapply(roads_ll, build_road_pts)

# Density kernels
CX <- ll_to_utm(-39.047, -14.789)$x
CY <- ll_to_utm(-39.047, -14.789)$y
M1X<- ll_to_utm(-39.020, -14.720)$x; M1Y<- ll_to_utm(-39.020,-14.720)$y
M2X<- ll_to_utm(-39.000, -14.860)$x; M2Y<- ll_to_utm(-39.000,-14.860)$y

hdi_density <- function(x, y) {
  alto  <- exp(-((x-CX)^2 + (y-CY)^2) / 2200^2)
  medio <- (exp(-((x-M1X)^2+(y-M1Y)^2)/3500^2) +
              exp(-((x-M2X)^2+(y-M2Y)^2)/3000^2)) * 0.35
  pmin(0.70*alto + 0.20*medio + 0.04, 1.)
}

road_proximity <- function(x, y) {
  val <- 0
  for (rname in names(road_pts)) {
    pts <- road_pts[[rname]]
    if (nrow(pts) == 0) next
    d_min <- min(sqrt((pts[,1]-x)^2 + (pts[,2]-y)^2))
    val   <- val + road_weight[rname] * exp(-d_min^2 / road_sigma[rname]^2)
  }
  min(val, 1.)
}

combined_density <- function(x, y) {
  min(road_proximity(x, y) * 0.60 + hdi_density(x, y) * 0.60, 1.)
}

sample_inside_muni <- function(n, density_fn, max_attempts = 400L) {
  pts <- matrix(NA_real_, nrow = n, ncol = 2)
  count <- 0L; attempts <- 0L
  while (count < n) {
    attempts <- attempts + 1L
    if (attempts > n * max_attempts) stop("Rejection sampling stalled")
    x <- runif(1, XMIN, XMAX)
    y <- runif(1, YMIN, YMAX)
    if (!point_in_muni(x, y)) next
    if (runif(1) < density_fn(x, y)) {
      count <- count + 1L
      pts[count, ] <- c(x, y)
    }
  }
  pts
}

cat("Sampling agent locations...\n")
hh_coords <- sample_inside_muni(N_HH, combined_density)
hh_x <- hh_coords[, 1]
hh_y <- hh_coords[, 2]
cat(sprintf("  %d HH placed.\n", N_HH))

# -�--�- Agent arrays -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
# N_SLOTS grows with births: each new HH gets a Vacant firm slot.
# N_SLOTS_0 records the initial count; used only where a frozen reference
# is genuinely needed (none currently).
N_SLOTS <- N_HH
N_SLOTS_0 <- N_SLOTS

hh_M          <- pmax(rnorm(N_HH, 60, 15), 0)
hh_L          <- numeric(N_HH)
hh_Employed   <- integer(N_HH)     # 0 = unemployed; else 1-based FirmID
hh_Wage       <- numeric(N_HH)
hh_des_C      <- numeric(N_HH)
hh_C_nominal  <- numeric(N_HH)
hh_firm_owner <- integer(N_HH)     # 0 = not owner; else own FirmID
hh_lambda     <- rep(LAMBDA0, N_HH)
hh_unemp_spells <- integer(N_HH)
hh_emp_spells   <- integer(N_HH)

# Firm arrays (same indexing as Python: 0-based internally, 1-based FirmID)
f_x           <- hh_x
f_y           <- hh_y
f_status      <- rep("Vacant", N_SLOTS)
f_age         <- integer(N_SLOTS)
f_M           <- numeric(N_SLOTS)
f_L           <- numeric(N_SLOTS)
f_L_max       <- numeric(N_SLOTS)
f_N           <- integer(N_SLOTS)
f_N_des       <- numeric(N_SLOTS)
f_Y           <- numeric(N_SLOTS)
f_S_hat       <- rep(1.0, N_SLOTS)
f_I_star      <- numeric(N_SLOTS)
f_Inv         <- numeric(N_SLOTS)
f_Inv_wip     <- numeric(N_SLOTS)   # production-sales lag: WIP not yet sellable
f_Sold        <- numeric(N_SLOTS)
f_Rev         <- numeric(N_SLOTS)
f_WB          <- numeric(N_SLOTS)
f_payroll_tax <- numeric(N_SLOTS)
f_Tax         <- numeric(N_SLOTS)
f_Profits_gross <- numeric(N_SLOTS)
f_Profits     <- numeric(N_SLOTS)
f_owner_pay   <- numeric(N_SLOTS)
f_m_dist      <- numeric(N_SLOTS)
f_Int         <- numeric(N_SLOTS)
f_inv_spend   <- numeric(N_SLOTS)
f_Profits_prev<- numeric(N_SLOTS)
f_Profit_smooth <- numeric(N_SLOTS)   # EWMA of profits (closure signal)
f_p           <- numeric(N_SLOTS)
f_mu          <- rep(0.25, N_SLOTS)
f_w           <- rep(MIN_WAGE, N_SLOTS)
f_a           <- rep(1.0, N_SLOTS)
f_a_prev      <- rep(1.0, N_SLOTS)   # TFP lagged one period (for surplus sharing)
f_OpenPos     <- integer(N_SLOTS)
f_entry_period <- integer(N_SLOTS)   # period in which the firm opened; 0 = never opened

reset_firm <- function(fi) {
  # fi is 1-based in R
  f_age[fi]    <<- 0L; f_N[fi]     <<- 0L; f_N_des[fi]  <<- 0
  f_Y[fi]      <<- 0;  f_S_hat[fi] <<- 1;  f_I_star[fi] <<- 0
  f_Inv[fi]    <<- 0;  f_Inv_wip[fi] <<- 0; f_Sold[fi]  <<- 0;  f_Rev[fi]    <<- 0
  f_WB[fi]     <<- 0;  f_payroll_tax[fi] <<- 0; f_Tax[fi] <<- 0
  f_Profits_gross[fi] <<- 0; f_Profits[fi] <<- 0
  f_owner_pay[fi] <<- 0; f_m_dist[fi] <<- 0
  f_Int[fi]    <<- 0;  f_L[fi]     <<- 0
  f_inv_spend[fi] <<- 0; f_Profits_prev[fi] <<- 0; f_Profit_smooth[fi] <<- 0
  f_p[fi]  <<- max(rnorm(1, 5.5, 2.0), P_FLOOR)
  f_mu[fi] <<- 0.35
  f_w[fi]  <<- max(rnorm(1, 10.0, 3.0), MIN_WAGE)
  # TFP at entry: set by caller via entry_a_override argument.
  f_a[fi]      <<- 1.0   # frontier technology at entry
  f_a_prev[fi] <<- 1.0
  f_OpenPos[fi] <<- 0L
  f_entry_period[fi] <<- 0L   # cleared on reset; set to current period at open
}

# -�--�- Distance matrix (pre-allocated to avoid O(N-�) copies each birth period) -�-
# N_MAX is an upper bound on population after MAXT periods of growth.
# All writes during birth events fill in-place; no matrix copy ever occurs.
cat("Building distance matrix...\n")
N_MAX  <- ceiling(N_HH * (1 + POP_GROWTH_RATE)^MAXT * 1.05) + MAXT

d_hf   <- matrix(0.0, N_MAX, N_MAX)
coords <- cbind(hh_x, hh_y)
d_hf[seq_len(N_HH), seq_len(N_HH)] <- as.matrix(dist(coords))
cat("Done.\n")

# -�--�- Period-1 seed -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
N_SEED    <- N_HH %/% 50L   # 2% seed rate; 10% caused catastrophic over-entry

openfirms<-f_status=="Open"
# if(length(openfirms) >0){
#   STARTUP_CAPITAL<-mean(f_w[f_w>0])*1.2
# }else{
STARTUP_CAPITAL <- abs(mean(hh_M)) * 0.3
# }

eligible<-which(hh_M>STARTUP_CAPITAL)

N_SEED<-pmin(length(eligible),N_SEED)

seed_slots <- sample(eligible, N_SEED, replace = FALSE)


for (fi in seed_slots) {
  
  hh_M[fi]         <- hh_M[fi] - STARTUP_CAPITAL
  f_M[fi]          <- STARTUP_CAPITAL  # owner equity only; bank loan added below
  f_L_max[fi]      <- STARTUP_CAPITAL * 3
  reset_firm(fi)
  f_status[fi]     <- "Open"
  hh_firm_owner[fi]<- fi
  f_entry_period[fi] <- 1L
}
cat(sprintf("Seeded %d firms at period 1.\n", N_SEED))

# ----- Bank initialisation (endogenous money) ----------------------------
# bank_NW = bank equity. Initial value sized to support seed-firm leverage.
# Lending headroom = bank_NW / CAR_MIN - total_outstanding_loans.
# With CAR_MIN=0.08, NW=1000 supports up to 12,500 in loans.
bank_NW <- BANK_NW_INIT_MULT * N_SEED * STARTUP_CAPITAL

# Grant leveraged startup loans to seed firms (credit creation)
for (fi in seed_slots) {
  bank_loan <- STARTUP_LEVERAGE * STARTUP_CAPITAL
  f_L[fi]   <- f_L[fi] + bank_loan    # bank asset (loan)
  f_M[fi]   <- f_M[fi] + bank_loan    # bank liability (deposit created)
}
cat(sprintf("Bank seeded: NW=%.0f, total seed loans=%.0f, headroom=%.0f\n",
            bank_NW, sum(f_L), max(0, bank_NW / CAR_MIN - sum(f_L))))

# ----- Government ----------------------------------------------------------
lg <- list(M = 100000, L = 0, tax_rev = 0, transfer = 0,
           deficit = 0, interest = 0) #lg$M
lg_g_spend <- 0
entry_cooldown_t <- 0L   # counts down after a mass-entry episode

# ----- INITIAL BALANCE SHEET CLOSURE (Godley sectoral balances) -----------
# The stock of private financial wealth = cumulative government deficit.
# All money in the economy was historically created by government spending
# financed through bank loans to government. This is the standard Post-
# Keynesian / MMT position (Godley 1999; Lavoie 2014, ch.1).
#
# Balance sheet identity at t=0:
#   bank_NW = sum(f_L) + lg$L - sum(hh_M) - sum(max(f_M,0)) - lg$M
#
# Solving for lg$L:
#   lg$L = bank_NW - sum(f_L) + sum(hh_M) + sum(max(f_M,0)) + lg$M
#
# Since RL_GOV = RD, the bank breaks even on government-backed deposits
# (same rate earned on gov debt as paid on gov-created deposits).
# The bank earns a positive spread ONLY on firm loans (RL_FIRM > RD).
{
  total_deposits_init <- sum(pmax(hh_M, 0)) + sum(pmax(f_M, 0)) + lg$M
  total_firm_loans    <- sum(f_L)
  lg$L <- bank_NW - total_firm_loans + total_deposits_init
  # Verify identity: bank_R should be 0 (no excess reserves at inception)
  bank_R_check <- bank_NW + total_deposits_init - total_firm_loans - lg$L
  cat(sprintf("Balance sheet closed: lg$L=%.0f, bank_R_check=%.2f (should be 0)\n",
              lg$L, bank_R_check))
  cat(sprintf("  Deposits=%.0f, Firm loans=%.0f, Gov debt=%.0f, Bank NW=%.0f\n",
              total_deposits_init, total_firm_loans, lg$L, bank_NW))
  cat(sprintf("  Gov int income (RL_GOV*lg$L): %.0f/period\n", RL_GOV * lg$L))
  cat(sprintf("  Deposit int cost (RD*deps):   %.0f/period\n", RD * total_deposits_init))
  cat(sprintf("  Firm int income (RL*f_L):     %.0f/period\n", RL_FIRM * total_firm_loans))
}
# -�--�- History -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
hist_keys <- c("unemp","wage","price","gdp","hM","fM","rw","gini_hh",
               "total_loans","n_open","n_vacant","n_openings","n_exits",
               "lg_M","lg_L","lg_tax","lg_transfer","lg_deficit","lg_interest","lg_debt_repay",
               "g_share_current","gov_transf_t",
               "owner_pay_total","m_dist_total","n_transfer_recip","transfer_total",
               "mean_M_bot30","mean_M_top30","median_wage_t","income_thresh_t",
               "supply","nom_demand","real_demand","unmet_nom",
               "inventory_total","avg_mu","n_hiring","utilisation",
               "price_index","inflation","inflation_ma5","gdp_real",
               "real_wage_index","inv_total","avg_a","avg_a_weighted",
               "n_hh","n_births","lg_g_spend","tax_rev_smooth",
               "mean_firm_age","med_firm_age","share_young_firms",
               "g_n_firms_served","total_demand_nom","g_share_demand",
               "bank_NW","bank_R","bank_new_loans","bank_dep_int_total",
               "bank_loan_int_total","bank_profit","bank_dividend",
               "bank_bad_debt","bank_headroom","total_M","sfc_check")
hist <- as.data.frame(matrix(0, nrow = MAXT, ncol = length(hist_keys)))
names(hist) <- hist_keys

# -�--�- Price index state -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
px_base_q   <- numeric(N_SLOTS)
px_base_p   <- numeric(N_SLOTS)
px_base_set <- FALSE
px_prev_pi  <- 1.0

# Fiscal stabiliser EWMA state — initialised to 0, warms up during burn-in.
# Updated every period regardless of FISCAL_STABILIZER toggle so the
# smoothed series is always available for diagnostics.
tax_rev_smooth <- 0

# Bolsa Família inflation correction state — price index at last uprating.
# Initialised to 1.0 (base period). Resets to current price_index each time
# a correction fires. Used to compute cumulative inflation since last update.
bf_last_pi <- 1.0

# Endogenous G_SHARE state — updated every period inside the loop.
G_SHARE    <- G_SHARE_INIT
# Target cash buffer: expected transfer bill at ~20% qualifying rate × multiplier.
# Expressed as a fixed monetary target so it doesn't shift with N_HH dynamics.
LG_M_TARGET <- GOV_TRANSF * N_HH * 0.15 * M_LG_TARGET_MULT
cat(sprintf("Fiscal targets: LG_M_TARGET=%.0f, G_SHARE_INIT=%.2f [%.2f, %.2f]\n",
            LG_M_TARGET, G_SHARE_INIT, G_SHARE_MIN, G_SHARE_MAX))

# ----- Banking sector state (endogenous money) ---------------------------
# bank_NW = bank net worth (equity). Primary state variable.
# Balance sheet identity (verified each period):
#   bank_R = bank_NW + total_deposits - total_loans
#   where total_deposits = sum(hh_M) + sum(max(f_M,0)) + lg$M
#   and   total_loans    = sum(f_L) + lg$L
# bank_NW already initialised in seed block above.
# Period-level flow accumulators (reset each period, recorded in hist):
bank_bad_debt_total   <- 0
bank_dep_int_total    <- 0
bank_loan_int_total   <- 0
bank_profit_period    <- 0
bank_dividend_period  <- 0

# -��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��
# MAIN LOOP
# -��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��-��
for (period in seq_len(MAXT)) {
  
  # hh_M <- hh_M + 5  # TEST — disabled; deposit interest now replaces this
  
  # Reset period-level bank accumulators
  bank_bad_debt_total  <- 0
  bank_dep_int_total   <- 0
  bank_loan_int_total  <- 0
  bank_profit_period   <- 0
  bank_dividend_period <- 0
  bank_new_loans       <- 0
  
  if (period %% 1 == 0) {
    n_om  <- sum(f_status == "Open")
    u_now <- 1 - sum(hh_Employed != 0L) / N_HH
    total_deposits <- sum(pmax(hh_M, 0)) + sum(pmax(f_M[f_status == "Open"], 0)) + max(lg$M, 0)
    total_loans_now <- sum(f_L) + lg$L
    bank_R_now <- bank_NW + total_deposits - total_loans_now
    headroom   <- max(0, bank_NW / CAR_MIN - sum(f_L))  # CAR on private loans only
    cat(sprintf("t=%3d | open=%3d | u=%.1f%% | HH_M=%.0f | F_M=%.0f | LG_M=%.0f LG_L=%.0f | B_NW=%.0f B_R=%.0f | headroom=%.0f | SC=%.0f\n",
                period, n_om, u_now*100, sum(hh_M), sum(f_M), lg$M, lg$L,
                bank_NW, bank_R_now, headroom, STARTUP_CAPITAL))
  }
  
  # -�--�- A. CLOSURE -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  n_exits <- 0L
  if (period > 1L) {
    om_mask  <- f_status == "Open"
    # Closure rule — two independent conditions, either sufficient:
    #
    # 1. Cash illiquidity: f_M - f_L ≤ 0
    #    Cash minus loan obligations ≤ 0. Inventory excluded (illiquid).
    #    Test: f_M - f_L ≤ 0 (cash net of loan obligations is non-positive).
    #
    # 2. Sustained operating loss (no f_L guard).
    #    Jovanovic (1982): if expected profit is persistently negative,
    #    exit dominates continuation.
    #    f_WB floor prevents spurious closure of zero-worker firms.
    #
    # 3. Grace period: firms with f_age ≤ CLOSURE_GRACE are exempt.
    #    With the production-sales lag, firms need at least 2 full periods
    #    before first revenue. CLOSURE_GRACE = 3 gives one buffer period
    #    beyond the minimum. Without it, every new firm is mechanically
    #    insolvent from birth (f_M - f_L ≤ 0 after paying wages with
    #    zero revenue) and killed before WIP ever matures.
    loss_thresh <- -LOSS_MULT * pmax(f_WB, MIN_WAGE)
    mature      <- f_age > CLOSURE_GRACE    # exempt young firms
    insolv <- om_mask & mature & (
      (f_M - f_L <= 0) |
      (f_Profit_smooth < loss_thresh)
    )
    ins_idx  <- which(insolv)
    for (fi in ins_idx) {
      # --- Bad debt resolution (endogenous money SFC-consistent) -------
      # 1. Firm's deposit (if positive) nets against its loan.
      # 2. Owner covers residual from personal deposits.
      # 3. Remaining uncovered loan = bad debt → bank_NW reduction.
      # No money vanishes: the bank's asset (f_L) is written off, reducing
      # bank equity. The deposit side shrinks by whatever f_M was positive.
      #
      # EQUITY RETURN: if the firm has positive net equity (f_M > f_L),
      # the excess belongs to the owner. Without this, deposits are
      # destroyed on closure → SFC leak (money vanishes from the system).
      if (f_M[fi] > f_L[fi]) {
        equity_return <- f_M[fi] - f_L[fi]
        hh_M[fi]     <- hh_M[fi] + equity_return   # owner receives net equity
        f_M[fi]      <- f_L[fi]                     # deposit shrinks to match loan
        # Both sides now equal → zeroing below destroys them symmetrically
      }
      net_exposure  <- max(0, f_L[fi] - max(f_M[fi], 0))  # loan minus deposit offset
      owner_covers  <- min(max(hh_M[fi], 0), net_exposure)
      hh_M[fi]      <- hh_M[fi] - owner_covers
      bad_debt      <- net_exposure - owner_covers
      bank_NW       <- bank_NW - bad_debt   # equity absorbs loss (can go negative)
      bank_bad_debt_total <- bank_bad_debt_total + bad_debt
      # Zero out firm financial position (now SFC-symmetric: f_M ≤ f_L)
      f_L[fi]           <- 0
      f_M[fi]           <- 0
      hh_L[fi]          <- 0
      hh_firm_owner[fi] <- 0L
      f_status[fi]      <- "Vacant"
      f_L_max[fi]       <- 0
    }
    n_exits <- length(ins_idx)
  }
  
  # Zero non-open fields
  not_om <- f_status != "Open"
  for (arr_name in c("f_M","f_L","f_L_max","f_N","f_N_des","f_Y","f_S_hat",
                     "f_I_star","f_Inv","f_Inv_wip","f_Sold","f_Rev","f_WB",
                     "f_payroll_tax","f_Tax","f_Profits_gross","f_Profits",
                     "f_owner_pay","f_m_dist","f_Int","f_p","f_w","f_a","f_OpenPos")) {
    arr <- get(arr_name)
    arr[not_om] <- 0
    assign(arr_name, arr)
  }
  f_S_hat[not_om] <- 1
  f_mu[not_om]    <- 0.25
  
  # Lay off workers at closed slots
  not_om_fids <- which(not_om)   # 1-based in R
  hh_Employed[hh_Employed %in% not_om_fids] <- 0L
  
  om      <- f_status == "Open"
  openidx <- which(om)
  f_age[om] <- f_age[om] + 1L
  
  # -�--�- A2. POPULATION GROWTH -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  # Births occur every period after burn-in. New entrants are placed inside
  # the municipality using the same density kernel as the original population.
  # Their wealth is drawn from the lower quartile of the current distribution
  # (lifecycle entry effect). All arrays are extended in place.
  n_births <- 0L
  if (POP_GROWTH_RATE > 0 && period > 1L) {
    n_births  <- max(1L, as.integer(ceiling(N_HH * POP_GROWTH_RATE)))
    N_HH_old  <- N_HH
    
    # Draw locations
    new_coords <- tryCatch(
      sample_inside_muni(n_births, combined_density, max_attempts = 200L),
      error = function(e) {
        idx <- sample.int(N_HH_old, n_births, replace = TRUE)
        cbind(hh_x[idx] + rnorm(n_births, 0, 200),
              hh_y[idx] + rnorm(n_births, 0, 200))
      }
    )
    
    # ----- SFC-consistent birth endowment: parental heritage ------------------
    # Two parents are sampled at random from the existing population.
    # Each transfers HERITAGE_FRAC of their wealth to the newborn.
    # This is debited from the parents, so no wealth is created ex nihilo
    # (De Nardi 2004; Piketty 2014 on inter-generational wealth transmission).
    # If the combined transfer falls below POP_WEALTH_FLOOR, the government
    # tops up the difference — a child benefit funded from lg$M.
    # Random parental selection (no assortative mating) implies the expected
    # endowment is 2 × HERITAGE_FRAC × mean(hh_M), generating wealth
    # persistence: children of wealthy parents start richer.
    new_M <- numeric(n_births)
    gov_topup_total <- 0
    for (bi in seq_len(n_births)) {
      p1 <- sample.int(N_HH_old, 1L)
      p2 <- sample.int(N_HH_old, 1L)
      transfer <- HERITAGE_FRAC * hh_M[p1] + HERITAGE_FRAC * hh_M[p2]
      hh_M[p1] <- hh_M[p1] - HERITAGE_FRAC * hh_M[p1]
      hh_M[p2] <- hh_M[p2] - HERITAGE_FRAC * hh_M[p2]
      if (transfer < POP_WEALTH_FLOOR) {
        topup <- POP_WEALTH_FLOOR - transfer
        lg$M  <- lg$M - topup
        if (lg$M < 0) { lg$L <- lg$L - lg$M; lg$M <- 0 }
        gov_topup_total <- gov_topup_total + topup
        transfer <- POP_WEALTH_FLOOR
      }
      new_M[bi] <- transfer
    }
    
    # -�--�- Distance matrix: fill new rows/cols in-place (no matrix copy) -�--�--�-
    # For each new birth, compute distances to all prior agents and write
    # directly into the pre-allocated d_hf. O(n_births * N_HH_old) total.
    for (bi in seq_len(n_births)) {
      new_idx  <- N_HH_old + bi
      # distances to the original N_HH_old households
      dx_e <- hh_x[seq_len(N_HH_old)] - new_coords[bi, 1]
      dy_e <- hh_y[seq_len(N_HH_old)] - new_coords[bi, 2]
      d_e  <- sqrt(dx_e^2 + dy_e^2)
      d_hf[new_idx, seq_len(N_HH_old)] <- d_e
      d_hf[seq_len(N_HH_old), new_idx] <- d_e
      # distances to earlier births in this same batch
      if (bi > 1L) {
        for (bj in seq_len(bi - 1L)) {
          prev_idx <- N_HH_old + bj
          d_nn <- sqrt((new_coords[bi, 1] - new_coords[bj, 1])^2 +
                         (new_coords[bi, 2] - new_coords[bj, 2])^2)
          d_hf[new_idx, prev_idx] <- d_nn
          d_hf[prev_idx, new_idx] <- d_nn
        }
      }
      # d_hf[new_idx, new_idx] remains 0 (self-distance, pre-filled)
    }
    
    # -�--�- HH arrays -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
    hh_x            <- c(hh_x,            new_coords[, 1])
    hh_y            <- c(hh_y,            new_coords[, 2])
    hh_M            <- c(hh_M,            new_M)
    hh_L            <- c(hh_L,            numeric(n_births))
    hh_Employed     <- c(hh_Employed,     integer(n_births))
    hh_Wage         <- c(hh_Wage,         numeric(n_births))
    hh_des_C        <- c(hh_des_C,        numeric(n_births))
    hh_C_nominal    <- c(hh_C_nominal,    numeric(n_births))
    hh_firm_owner   <- c(hh_firm_owner,   integer(n_births))
    hh_lambda       <- c(hh_lambda,       rep(LAMBDA0, n_births))
    hh_unemp_spells <- c(hh_unemp_spells, integer(n_births))
    hh_emp_spells   <- c(hh_emp_spells,   integer(n_births))
    
    # -�--�- Firm arrays (Vacant slots for new citizens) -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
    f_x             <- c(f_x,             new_coords[, 1])
    f_y             <- c(f_y,             new_coords[, 2])
    f_status        <- c(f_status,        rep("Vacant", n_births))
    f_age           <- c(f_age,           integer(n_births))
    f_M             <- c(f_M,             numeric(n_births))
    f_L             <- c(f_L,             numeric(n_births))
    f_L_max         <- c(f_L_max,         numeric(n_births))
    f_N             <- c(f_N,             integer(n_births))
    f_N_des         <- c(f_N_des,         numeric(n_births))
    f_Y             <- c(f_Y,             numeric(n_births))
    f_S_hat         <- c(f_S_hat,         rep(1.0, n_births))
    f_I_star        <- c(f_I_star,        numeric(n_births))
    f_Inv           <- c(f_Inv,           numeric(n_births))
    f_Inv_wip       <- c(f_Inv_wip,       numeric(n_births))
    f_Sold          <- c(f_Sold,          numeric(n_births))
    f_Rev           <- c(f_Rev,           numeric(n_births))
    f_WB            <- c(f_WB,            numeric(n_births))
    f_payroll_tax   <- c(f_payroll_tax,   numeric(n_births))
    f_Tax           <- c(f_Tax,           numeric(n_births))
    f_Profits_gross <- c(f_Profits_gross, numeric(n_births))
    f_Profits       <- c(f_Profits,       numeric(n_births))
    f_owner_pay     <- c(f_owner_pay,     numeric(n_births))
    f_m_dist        <- c(f_m_dist,        numeric(n_births))
    f_Int           <- c(f_Int,           numeric(n_births))
    f_inv_spend     <- c(f_inv_spend,     numeric(n_births))
    f_Profits_prev  <- c(f_Profits_prev,  numeric(n_births))
    f_Profit_smooth <- c(f_Profit_smooth, numeric(n_births))
    f_p             <- c(f_p,             numeric(n_births))
    f_mu            <- c(f_mu,            rep(0.25, n_births))
    f_w             <- c(f_w,             rep(MIN_WAGE, n_births))
    f_a             <- c(f_a,             rep(1.0, n_births))
    f_a_prev        <- c(f_a_prev,        rep(1.0, n_births))
    f_OpenPos       <- c(f_OpenPos,       integer(n_births))
    f_entry_period  <- c(f_entry_period,  integer(n_births))
    
    # -�--�- Laspeyres price index state vectors -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
    # px_base_q and px_base_p must stay co-indexed with firm arrays.
    # New slots have no base-period observation yet (zero = not yet set).
    px_base_q <- c(px_base_q, numeric(n_births))
    px_base_p <- c(px_base_p, numeric(n_births))
    
    # -�--�- Dimension counters -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
    N_HH   <- N_HH_old + n_births
    N_SLOTS <- N_HH
  }
  
  # -�--�- B. LOGISTIC ENTRY -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  # All HH indices up to N_SLOTS (grows with births) are eligible.
  # New citizens have a Vacant firm slot allocated at birth.
  
  
  # Dynamic startup capital — EWMA-smoothed to prevent wild oscillation.
  # Target = mean wage bill per firm; falls back to SC_WEALTH_FRAC * mean_M.
  # HARD FLOOR: SC never drops below SC_FLOOR_TRANSF_MULT * GOV_TRANSF.
  # Without this floor, both the wage-bill target and the mean-wealth fallback
  # collapse together in downturns, making the entire population eligible and
  # triggering an explosive entry boom that destroys demand for all incumbents.
  # The floor is theoretically grounded: entrepreneurship requires a working-
  # capital buffer sufficient to cover at least a few wage payments; tying it
  # to GOV_TRANSF (the reservation wage) ensures it scales with the economy.
  SC_FLOOR <- SC_FLOOR_TRANSF_MULT * GOV_TRANSF
  {
    open_wages  <- f_w[f_status == "Open" & f_w > 0]
    open_N      <- f_N[f_status == "Open" & f_N > 0]
    mean_hh_M   <- max(mean(hh_M), MIN_WAGE, na.rm = TRUE)
    SC_fallback <- max(mean_hh_M * SC_WEALTH_FRAC, SC_FLOOR)
    SC_target   <- if (length(open_wages) > 0 && length(open_N) > 0)
      mean(open_wages) * pmax(mean(open_N)*0.6, 1)
    else
      SC_fallback
    SC_target <- max(SC_target, SC_fallback, SC_FLOOR, na.rm = TRUE)
    if (!is.finite(SC_target)) SC_target <- SC_fallback
    STARTUP_CAPITAL <- 0.5 * STARTUP_CAPITAL + 0.5 * SC_target
  }
  
  # ----- GOV_TRANSF: fixed nominal; WEALTH_FLOOR tracks median wage ----------
  # GOV_TRANSF is NOT updated here — nominal value is set by the periodic
  # inflation correction rule (after price_index is computed, block K).
  # WEALTH_FLOOR still tracks median wages to preserve the eligibility logic:
  # households below BF_WEALTH_FRAC × median_wage qualify for transfers.
  {
    emp_wages_pre   <- hh_Wage[hh_Employed != 0L]
    median_wage_pre <- if (length(emp_wages_pre) > 0) median(emp_wages_pre) else MIN_WAGE
    WEALTH_FLOOR    <- max(BF_WEALTH_FRAC * median_wage_pre, GOV_TRANSF * 6)
  }
  
  
  
  eligible <- which(f_status == "Vacant" &
                      hh_firm_owner == 0L &
                      hh_M >= STARTUP_CAPITAL)
  n_openings <- 0L
  
  if (period > 1L && length(eligible) > 0) {
    t_wealth <- BETA_M * hh_M[eligible]
    t_unemp  <- BETA_UNEMP * as.numeric(hh_Employed[eligible] == 0L)
    
    emp_wages  <- hh_Wage[hh_Employed != 0L]
    true_med_w <- if (length(emp_wages) > 0) median(emp_wages) else MIN_WAGE
    noise      <- rlnorm(length(eligible), 0, SIGMA_W)
    w_perceived<- true_med_w * noise
    w_own      <- ifelse(hh_Employed[eligible] != 0L, hh_Wage[eligible], 0)
    t_gain     <- BETA_GAIN * (w_perceived - w_own)
    
    t_uspells  <- BETA_USPELLS * pmin(hh_unemp_spells[eligible], USPELLS_CAP)
    t_espells  <- BETA_ESPELLS * hh_emp_spells[eligible]
    
    om_now_idx <- openidx
    if (length(om_now_idx) > 0) {
      d_to_open <- d_hf[eligible, om_now_idx, drop = FALSE]
      n_local   <- rowSums(d_to_open <= SAT_RADIUS)
    } else {
      n_local <- rep(0, length(eligible))
    }
    t_sat <- BETA_SAT * n_local
    
    logit_score <- INTERCEPT_OPEN + t_wealth + t_unemp + t_gain +
      t_uspells + t_espells + t_sat
    p_open  <- sigmoid(logit_score)
    do_open <- eligible[rbinom(length(eligible), 1L, p_open) == 1L]
    
    # Hard cap with cooldown: after a mass-entry episode, halve the cap
    # for ENTRY_COOLDOWN periods to smooth the trap-to-recovery transition.
    if (entry_cooldown_t > 0L) entry_cooldown_t <- entry_cooldown_t - 1L
    max_new <- max(1L, as.integer(floor(N_HH * ENTRY_RATE_CAP)))
    if (entry_cooldown_t > 0L) max_new <- max(1L, max_new %/% 2L)
    if (length(do_open) > max_new)
      do_open <- sample(do_open, max_new)
    if (length(do_open) >= max_new) entry_cooldown_t <- ENTRY_COOLDOWN
    
    
    
    for (fi in do_open) {
      # --- Leveraged entry (Kalecki 1937, principle of increasing risk) ---
      # Owner contributes equity from personal wealth.
      # Bank matches with a loan = STARTUP_LEVERAGE * equity, subject to CAR.
      # Credit creation: f_M ↑ (new deposit), f_L ↑ (new loan), bank_NW unchanged.
      owner_equity  <- STARTUP_CAPITAL
      bank_loan_req <- STARTUP_LEVERAGE * owner_equity
      lending_headroom_now <- max(0, bank_NW / CAR_MIN - sum(f_L))  # lg$L excluded: sovereign risk-weight = 0 (Basel II/III)
      bank_loan     <- min(bank_loan_req, lending_headroom_now)
      
      hh_M[fi]          <- hh_M[fi] - owner_equity
      reset_firm(fi)             # zeros f_L, f_M, sets f_a=1.0, draws f_p, f_w
      f_M[fi]           <- owner_equity + bank_loan   # deposit = equity + credit
      f_L[fi]           <- bank_loan                  # loan obligation
      f_L_max[fi]       <- (owner_equity + bank_loan) * 3
      f_status[fi]      <- "Open"
      hh_firm_owner[fi] <- fi
      f_entry_period[fi] <- period
      hh_unemp_spells[fi] <- 0L
      hh_emp_spells[fi]   <- 0L
      n_openings <- n_openings + 1L
      bank_new_loans <- bank_new_loans + bank_loan
    }
  }
  
  om      <- f_status == "Open"
  openidx <- which(om)
  
  
  # -�--�- ENDOGENOUS ANIMAL SPIRITS -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  # Floor on sales expectations scales with mean household income,
  # so firm optimism is anchored to aggregate purchasing power.
  # AS_SCALE controls the sensitivity; 0.5-1.5 is a reasonable range.
  # Keynes (1936) ch.12; Dosi et al. (2010) for ABM precedent.
  # Guard: falls back to MIN_WAGE when no one is employed yet (period 1)
  # or when the wage vector is transiently empty after a crash.
  {
    pos_wages <- hh_Wage[hh_Wage > 0]
    animal_spirits_t <- if (length(pos_wages) > 0)
      AS_SCALE * mean(pos_wages)
    else
      AS_SCALE * ANIMAL_SPIRITS_MIN
    animal_spirits_t <- max(animal_spirits_t, ANIMAL_SPIRITS_MIN, na.rm = TRUE)
  }
  
  
  # ----- WIP MATURATION (before planning) ----------------------------------
  # Last period's work-in-progress becomes sellable inventory.
  # Placed here so that the planning block (C) sees the firm's full
  # sellable stock when computing desired production and markup.
  if (any(om)) {
    f_Inv[om] <- f_Inv[om] + f_Inv_wip[om]
    f_Inv_wip[om] <- 0
  }

  # -�--�- C. PLANNING -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  if (any(om)) {
    s <- ewma_fn(f_S_hat[om], f_Sold[om], ALPHA_S)
    #s <- pmax(s, ANIMAL_SPIRITS)
    s <- pmax(s, animal_spirits_t)  
    f_S_hat[om] <- s
    Is <- THETA * s; f_I_star[om] <- Is
    Yd <- pmax(s + (Is - f_Inv[om]), 0)
    f_N_des[om] <- ceiling((1 - LAMBDA_N) * f_N[om] +
                             LAMBDA_N * pmax(Yd / (f_a[om] + 1e-9), 0))
    f_mu[om] <- clamp(f_mu[om] + MU_PHI * (Is - f_Inv[om]), MU_MIN, MU_MAX)
    f_N_des[!om] <- 0
    
    f_OpenPos <- pmax(0L, as.integer(round(f_N_des - f_N)))
    f_OpenPos[!om] <- 0L
    f_OpenPos[f_M <= 0] <- 0L
    
    if (period > 1L) {
      t_tight <- f_N[om] / (f_N_des[om] + 1e-5)
      t_tight[abs(t_tight - 1) <= 0.1] <- 1.0
      
      
      
      if(period %% wageperiod == 0){
        f_w[om] <- clamp(
          f_w[om] * clamp(1 + PHI_W * (1 - t_tight), 0.75, 1.25),
          MIN_WAGE, MAX_WAGE)
      }
    }
  }
  
  # -�--�- D. LABOUR MARKET -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  for (i in sample.int(N_HH)) {
    hr  <- which(f_OpenPos > 0L & f_status == "Open")
    cf  <- hh_Employed[i]
    d_row <- d_hf[i, ]
    lam_i <- hh_lambda[i]
    
    if (cf != 0L) {
      ci <- cf
      if (f_status[ci] != "Open") {
        hh_Employed[i] <- 0L; hh_Wage[i] <- 0; next
      }
      alt <- spatial_search(d_row, setdiff(hr, ci), lam_i)
      A   <- f_w[ci] - d_row[ci] * COST_PER_M - GOV_TRANSF
      if (runif(1) < sigmoid((-A - QUIT_FRICTION) / QUIT_SIGMA)) {
        hh_Employed[i] <- 0L; hh_Wage[i] <- 0
        f_OpenPos[ci]  <- f_OpenPos[ci] + 1L; next
      }
      if (length(alt) == 0) next
      if (length(alt) > MAX_CONSIDERED)
        alt <- alt[order(d_row[alt])][seq_len(MAX_CONSIDERED)]
      Ac     <- c(0, (f_w[alt] - f_w[ci]) - (d_row[alt] - d_row[ci]) * COST_PER_M)
      ch     <- c(ci, alt)
      chosen <- ch[softmax_sample(Ac, SWITCH_SIGMA)]
      if (chosen != ci && f_OpenPos[chosen] > 0L) {
        hh_Employed[i] <- chosen; hh_Wage[i] <- f_w[chosen]
        f_OpenPos[ci]  <- f_OpenPos[ci] + 1L
        f_OpenPos[chosen] <- f_OpenPos[chosen] - 1L
      }
    } else {
      offers <- spatial_search(d_row, hr, lam_i)
      if (length(offers) == 0) next
      if (length(offers) > MAX_CONSIDERED)
        offers <- offers[order(d_row[offers])][seq_len(MAX_CONSIDERED)]
      A_alt  <- f_w[offers] - d_row[offers] * COST_PER_M - GOV_TRANSF
      if (runif(1) >= sigmoid(max(A_alt) / ACCEPT_SIGMA)) next
      chosen <- offers[softmax_sample(A_alt, SWITCH_SIGMA)]
      if (f_OpenPos[chosen] > 0L) {
        hh_Employed[i] <- chosen; hh_Wage[i] <- f_w[chosen]
        f_OpenPos[chosen] <- f_OpenPos[chosen] - 1L
      }
    }
    f_OpenPos <- pmax(f_OpenPos, 0L)
  }
  
  # Recount employment
  emp_vec <- hh_Employed[hh_Employed > 0L]
  f_N <- tabulate(emp_vec, nbins = N_SLOTS)
  
  # -�--�- E. LAYOFFS -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  if (period > 1L) {
    excess_idx <- which(om & (f_N - f_N_des) > 0)
    for (fi in excess_idx) {
      fid  <- fi
      wks  <- which(hh_Employed == fid)
      n_fire <- min(length(wks), as.integer(f_N[fi] - f_N_des[fi]))
      if (n_fire > 0) {
        fired <- if (length(wks) == 1) wks else sample(wks, n_fire)
        hh_Employed[fired] <- 0L
      }
    }
    emp_vec <- hh_Employed[hh_Employed > 0L]
    f_N <- tabulate(emp_vec, nbins = N_SLOTS)
  }
  
  # Spell counters
  is_emp  <- hh_Employed != 0L
  hh_emp_spells[is_emp]    <- hh_emp_spells[is_emp] + 1L
  hh_unemp_spells[is_emp]  <- 0L
  non_owner_unemp <- !is_emp & hh_firm_owner == 0L
  hh_unemp_spells[non_owner_unemp] <- hh_unemp_spells[non_owner_unemp] + 1L
  hh_emp_spells[non_owner_unemp]   <- 0L
  
  # ----- F. PRODUCTION & WAGES --------------------------------------------
  # PRODUCTION-SALES LAG (Graziani 2003; Godley & Lavoie 2007, ch.7-11)
  # Output produced this period goes into f_Inv_wip (work-in-progress).
  # WIP is NOT available for sale until next period, when it transfers
  # to f_Inv (sellable). This creates a structural working capital demand:
  # firms must pay wages THIS period but earn revenue NEXT period.
  # The resulting cash-flow gap is the "finance motive" for credit
  # (Keynes 1937, "Alternative Theories of the Rate of Interest").
  #
  # Timing each period:
  #   0. WIP maturation (done before planning, block C)
  #   1. Produce f_Y -> f_Inv_wip (not yet sellable)
  #   2. Pay wages from f_M (-> cash-flow gap -> bank lending)
  #   3. Sell from f_Inv (last period's output)
  #   4. Receive revenue -> f_M recovers

  f_Y[] <- 0; f_WB[] <- 0; f_Sold[] <- 0
  # ----- TFP: learning-by-doing (Arrow 1962) ------------------------------
  f_inv_spend[] <- 0; total_inv <- 0
  if (any(om)) {
    firm_lbd <- rnorm(sum(om), mean = LBD_RATE, sd = LBD_SIGMA)
    f_a[om]  <- pmax(f_a[om] + firm_lbd, A_MIN)
  }

  # Step 1: Production -> WIP (not yet sellable)
  Y_cap   <- f_a * f_N
  Y_plan  <- pmax(f_S_hat + (f_I_star - f_Inv), 0)
  f_Y[om] <- pmin(Y_plan[om], Y_cap[om])
  f_Inv_wip[om] <- f_Y[om]   # goes to WIP, available for sale NEXT period

  # Step 2: Pay wages (creates the cash-flow gap that drives credit demand)
  for (fi in openidx) {
    wks <- which(hh_Employed == fi)
    wb  <- length(wks) * f_w[fi]
    f_WB[fi]  <- wb
    f_M[fi]   <- f_M[fi] - wb
    if (length(wks) > 0) {
      hh_M[wks]   <- hh_M[wks] + f_w[fi]
      hh_Wage[wks]<- f_w[fi]
    }
  }
  
  # -�--�- G. PRICING -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  pr <- which(om & f_Y > 0)
  sr <- which(om & f_Y == 0 & f_Inv > 0)
  
  
  if(period %% priceperiod == 0){
    
    if (length(pr) > 0) f_p[pr] <- (1 + f_mu[pr]) * (f_w[pr] / f_a[pr])
    if (length(sr) > 0) {
      sig_sr <- pmin(f_Inv[sr] / (f_S_hat[sr] + 1e-5), 5)
      f_p[sr] <- pmax(f_p[sr] * (1 - GAMMA_MD * sig_sr), P_FLOOR)
    }
  }
  
  # ----- H. BANKING SECTOR (endogenous money, SFC-closed) -----------------
  # All financial flows concentrated here. Sequence:
  #   1. Deposit interest paid to all depositors (HH, firms, gov)
  #   2. Loan interest computed (for profit chain in J) + bank_NW credited
  #   3. Government loan interest
  #   4. New lending to cash-negative firms (credit creation, CAR-constrained)
  #
  # Amortisation and dividends are handled post-accounting (after block J)
  # so that repayment can be funded from period profits.

  # Step 1: Deposit interest — bank pays depositors
  # Depositors receive RD * max(balance, 0). Money is created on the bank's
  # balance sheet: deposit liability ↑, bank_NW ↓ by total interest paid.
  # This is the primary channel replacing the +5 helicopter test.
  dep_int_hh <- RD * pmax(hh_M, 0)
  dep_int_f  <- numeric(N_SLOTS); dep_int_f[om] <- RD * pmax(f_M[om], 0)
  dep_int_lg <- RD * max(lg$M, 0)
  total_dep_int <- sum(dep_int_hh) + sum(dep_int_f) + dep_int_lg

  hh_M    <- hh_M + dep_int_hh        # HH deposits grow
  f_M[om] <- f_M[om] + dep_int_f[om]  # firm deposits grow
  lg$M    <- lg$M + dep_int_lg         # gov deposit grows
  bank_NW <- bank_NW - total_dep_int   # bank equity shrinks
  bank_dep_int_total <- total_dep_int

  # Step 2: Loan interest — firms pay bank (computed for accounting block J)
  # f_Int enters the profit chain: f_Profits_gross = f_Rev - f_WB - f_Int.
  # The actual f_M deduction happens in J via the cash-flow update.
  # bank_NW is credited here; the corresponding f_M deduction in J balances it.
  f_Int[] <- 0
  f_Int[om] <- f_L[om] * RL_FIRM
  bank_NW   <- bank_NW + sum(f_Int[om])
  bank_loan_int_total <- sum(f_Int[om])

  # Step 3: Government loan interest — debited from lg$M, credited to bank_NW
  # Replaces the old lg_int → spray-to-households route.
  # Government interest flows to bank equity, not directly to HH.
  # HH receive their share via deposit interest (step 1) and bank dividends.
  lg_int   <- lg$L * RL_GOV
  lg$M     <- lg$M - lg_int
  if (lg$M < 0) { lg$L <- lg$L - lg$M; lg$M <- 0 }
  bank_NW  <- bank_NW + lg_int
  bank_loan_int_total <- bank_loan_int_total + lg_int

  # Step 4: New lending — credit creation for cash-negative firms
  # After wages (block F), some firms have f_M < 0. The bank creates
  # a deposit (f_M ↑) backed by a loan (f_L ↑). Bank_NW unchanged.
  # Lending is CAR-constrained: total loans ≤ bank_NW / CAR_MIN.
  lending_headroom <- max(0, bank_NW / CAR_MIN - sum(f_L))  # lg$L excluded: sovereign risk-weight = 0 (Basel II/III)
  neg_idx <- which(om & f_M < 0)
  for (fi in neg_idx) {
    need <- -f_M[fi]
    lend <- min(need, lending_headroom)
    if (lend > 0) {
      f_L[fi]          <- f_L[fi] + lend   # bank asset (loan)
      f_M[fi]          <- f_M[fi] + lend   # bank liability (deposit created)
      lending_headroom <- lending_headroom - lend
      bank_new_loans   <- bank_new_loans + lend
    }
    # If headroom exhausted, residual negative f_M triggers insolvency next period.
  }

  # ----- I. GOODS MARKET --------------------------------------------------
  f_Inv[!is.finite(f_Inv)] <- 0
  f_p[!is.finite(f_p) | f_p < 0] <- 0
  hh_des_C   <- pmax(0, A1 * hh_M + A2 * hh_Wage)
  hh_C_nominal[] <- 0

  for (hi in sample.int(N_HH)) {
    budget <- min(max(0, hh_des_C[hi]), max(0, hh_M[hi]))
    if (budget <= 0) next
    d_om <- d_hf[hi, openidx]
    d_om[!is.finite(d_om)] <- Inf
    lam  <- hh_lambda[hi]
    vis  <- integer(0)
    while (length(vis) == 0 && lam >= 1e-9) {
      hit <- rbinom(length(openidx), 1L, exp(-lam * d_om)) == 1L
      vis <- openidx[hit]; lam <- lam * 0.3
    }
    if (length(vis) == 0) next
    vis <- vis[order(d_hf[hi, vis])][seq_len(min(length(vis), MAX_VISITS))]
    spent <- 0
    for (j in vis) {
      if (budget <= 0.01) break
      sq <- max(0, f_Inv[j]); p_j <- f_p[j]
      if (sq <= 0 || !is.finite(p_j) || p_j <= 0.001) next
      buy_q <- min(sq, budget / p_j); pay <- p_j * buy_q
      if (pay > 0) {
        hh_M[hi] <- hh_M[hi] - pay; spent <- spent + pay; budget <- budget - pay
        f_Inv[j] <- max(0, f_Inv[j] - buy_q); f_Sold[j] <- f_Sold[j] + buy_q
      }
    }
    hh_C_nominal[hi] <- spent
  }

  # ----- J. ACCOUNTING ----------------------------------------------------
  f_Rev[om]           <- f_Sold[om] * f_p[om]
  f_Profits_gross[om] <- f_Rev[om] - f_WB[om] - f_Int[om]
  f_payroll_tax[om]   <- f_WB[om] * PAYROLL_TAX
  taxable_om          <- pmax(f_Profits_gross[om] - f_payroll_tax[om], 0)
  f_Tax[om]           <- taxable_om * PROFIT_TAX
  f_Profits[om]       <- f_Profits_gross[om] - f_payroll_tax[om] - f_Tax[om]

  lg_tax_now <- sum(f_payroll_tax) + sum(f_Tax)
  lg$M       <- lg$M + lg_tax_now
  lg$tax_rev <- lg_tax_now

  # Owner profit share
  f_owner_pay[] <- 0; total_owner_pay <- 0
  for (fi in openidx) {
    fp <- f_Profits[fi]
    if (fp <= 0) next
    pay <- fp * OWNER_PROFIT_SHARE
    f_owner_pay[fi] <- pay
    hh_M[fi] <- hh_M[fi] + pay
    total_owner_pay <- total_owner_pay + pay
  }

  # Excess-cash dividend
  for (fi in openidx) {
    wb_fi  <- f_WB[fi]
    cap_fi <- FM_CAP_WAGE_PERIODS * max(wb_fi, MIN_WAGE)
    excess <- f_M[fi] - cap_fi
    if (excess <= 0) next
    f_M[fi]  <- f_M[fi] - excess
    hh_M[fi] <- hh_M[fi] + excess
    total_owner_pay <- total_owner_pay + excess
  }

  f_Profits_prev[openidx] <- f_Profits[openidx]
  f_Profits_prev[!om]     <- 0
  # Update EWMA profitability signal for open firms
  f_Profit_smooth[om]  <- (1 - ALPHA_PROF) * f_Profit_smooth[om] +
                           ALPHA_PROF       * f_Profits[om]
  f_Profit_smooth[!om] <- 0

  # Wage-productivity linkage
  if (any(om)) {
    if (USE_EFF_WAGE) {
      mean_a_now <- mean(f_a[om])
      prod_floor <- MIN_WAGE * (f_a[om] / (mean_a_now + 1e-9)) ^ ETA_W
      f_w[om]    <- pmax(f_w[om], prod_floor)
    }
    if (USE_SURPLUS) {
      delta_a  <- pmax(f_a[om] - f_a_prev[om], 0)
      f_w[om]  <- f_w[om] + THETA_W * delta_a * f_w[om]
      f_w[om]  <- clamp(f_w[om], MIN_WAGE, MAX_WAGE)
    }
    f_a_prev[om]  <- f_a[om]
    f_a_prev[!om] <- 0
  }

  # Firm cash-flow (single application — f_Int already in f_Profits_gross)
  for (fi in openidx)
    f_M[fi] <- f_M[fi] + (f_Profits[fi] - f_owner_pay[fi]) + f_WB[fi]

  # Loan amortisation — scheduled principal repayment (Minsky 1986)
  # Each period, the firm owes f_L / LOAN_MATURITY in principal.
  # If cash-short, the shortfall is capitalised (rolled into f_L),
  # creating the debt-trap dynamics that drive endogenous fragility.
  # Both sides of the bank balance sheet shrink: money is DESTROYED.
  amort_idx <- which(om & f_L > 0)
  for (fi in amort_idx) {
    instalment   <- f_L[fi] / LOAN_MATURITY
    actual_repay <- min(instalment, max(f_M[fi], 0))
    f_M[fi]      <- f_M[fi] - actual_repay    # deposit shrinks (money destroyed)
    f_L[fi]      <- f_L[fi] - actual_repay    # loan shrinks
    shortfall    <- instalment - actual_repay
    if (shortfall > 0) f_L[fi] <- f_L[fi] + shortfall  # capitalised interest
  }

  # Bank dividends — distribute net interest income to HH every period.
  # bank_profit = loan interest received - deposit interest paid.
  # (bad debt already hit bank_NW in closure block A)
  # Dividend suspension: if bank_NW < BANK_NW_FLOOR, full profit is retained
  # to rebuild equity (Basel III Pillar 2 supervisory logic; see BCBS 2011).
  # When bank_NW >= BANK_NW_FLOOR, BANK_DIV_RATE fraction is paid out.
  bank_profit_period   <- bank_loan_int_total - bank_dep_int_total
  bank_dividend_period <- if (bank_NW >= BANK_NW_FLOOR) {
    max(0, bank_profit_period) * BANK_DIV_RATE
  } else 0
  if (bank_dividend_period > 0 && N_HH > 0) {
    hh_M    <- hh_M + bank_dividend_period / N_HH   # deposit liability ↑
    bank_NW <- bank_NW - bank_dividend_period        # equity ↓
  }

  # M distribution
  f_m_dist[] <- 0; total_m_dist <- 0
  if (period %% M_DIST_PERIOD == 0) {
    for (fi in openidx) {
      if (f_M[fi] <= 0) next
      wks <- which(hh_Employed == fi)
      if (length(wks) == 0) next
      pool <- f_M[fi] * M_DIST_RATE
      f_m_dist[fi] <- pool; f_M[fi] <- f_M[fi] - pool
      hh_M[wks] <- hh_M[wks] + pool / length(wks)
      total_m_dist <- total_m_dist + pool
    }
  }

  # ----- K. BOLSA FAMILIA --------------------------------------------------
  emp_wages_now  <- hh_Wage[hh_Employed != 0L]
  median_wage_t  <- if (length(emp_wages_now) > 0) median(emp_wages_now) else MIN_WAGE
  income_thresh  <- INCOME_FRAC * median_wage_t
  income_ok      <- (hh_Wage < income_thresh) | (hh_Employed == 0L)
  wealth_ok      <- hh_M < WEALTH_FLOOR
  qualifying     <- which(wealth_ok & income_ok)
  n_qual         <- length(qualifying)
  transfer_bill  <- n_qual * GOV_TRANSF
  lg$M           <- lg$M - transfer_bill
  hh_M[qualifying] <- hh_M[qualifying] + GOV_TRANSF
  unemp_qual     <- qualifying[hh_Employed[qualifying] == 0L]
  hh_Wage[unemp_qual] <- GOV_TRANSF

  if (lg$M < 0) { lg$L <- lg$L - lg$M; lg$M <- 0 }
  # lg_int computed and debited from lg$M in banking block H (step 3).
  # Government interest flows to bank_NW, not directly to HH.
  # HH receive their share via deposit interest and bank dividends.
  lg$interest <- lg_int
  lg$transfer <- transfer_bill

  lg_g_spend      <- 0
  g_n_firms_served <- 0L

  if (length(openidx) > 0) {
    # ----- Endogenous G_SHARE update (stock-flow norm) ----------------------
    # Adjust G_SHARE toward the cash buffer target before computing the budget.
    # Division by (tax_rev_smooth + 1) keeps the adjustment scale-invariant:
    # a given cash surplus triggers a proportional spending response regardless
    # of the absolute revenue level. Guard against tax_rev_smooth = 0 at t=1.
    LG_M_TARGET <- GOV_TRANSF * N_HH * 0.15 * M_LG_TARGET_MULT
    g_share_adj <- G_SHARE_PHI * (lg$M - LG_M_TARGET) / (tax_rev_smooth + 1)
    G_SHARE     <- clamp(G_SHARE + g_share_adj, G_SHARE_MIN, G_SHARE_MAX)

    gov_budget <- if (FISCAL_STABILIZER) {
      G_SHARE * tax_rev_smooth
    } else {
      net_fiscal <- lg_tax_now - transfer_bill
      G_SHARE * max(net_fiscal, 0)
    }

    if (gov_budget > 0.01) {
      # --- Proportional procurement -----------------------------------------
      # Budget is allocated across ALL firms in proportion to their available
      # inventory value. This prevents the first few firms in a random walk
      # from absorbing the entire budget, spreads demand support broadly, and
      # avoids the fM concentration that was suppressing the Keynesian multiplier.
      #
      # For each firm j: alloc_j = gov_budget * (f_Inv[j]*f_p[j]) / total_inv_val
      # We then buy up to alloc_j / p_j units from firm j.
      # Any unspent allocation (e.g. firm has less inventory than allocated)
      # is NOT redistributed — it is simply unspent, keeping the rule simple
      # and avoiding a second-pass loop.
      #
      # NOTE: f_M[j] is NOT updated here. Revenue enters f_M via the accounting
      # step below (f_Rev = f_Sold * f_p), exactly as HH purchases do.
      # The previous code updated f_M[j] directly AND via f_Sold, double-crediting
      # gov procurement revenue. That was the primary cause of fM bloat.
      inv_val <- pmax(f_Inv[openidx] * f_p[openidx], 0)
      total_inv_val <- sum(inv_val)

      if (total_inv_val > 0) {
        alloc <- gov_budget * inv_val / total_inv_val   # budget share per firm
        for (k in seq_along(openidx)) {
          j     <- openidx[k]
          sq    <- max(0, f_Inv[j])
          p_j   <- f_p[j]
          if (sq <= 0 || !is.finite(p_j) || p_j <= 0.001 || alloc[k] <= 0.01) next
          buy_q <- min(sq, alloc[k] / p_j)
          pay   <- p_j * buy_q
          if (pay <= 0) next
          lg$M        <- lg$M - pay
          lg_g_spend  <- lg_g_spend + pay
          f_Inv[j]    <- max(0, f_Inv[j] - buy_q)
          f_Sold[j]   <- f_Sold[j] + buy_q   # revenue flows via accounting below
          g_n_firms_served <- g_n_firms_served + 1L
        }
        if (lg$M < 0) { lg$L <- lg$L - lg$M; lg$M <- 0 }
      }
    }
  }
  lg$deficit <- transfer_bill + lg_int + lg_g_spend - lg_tax_now

  # ----- Sovereign debt repayment (SFC-consistent) -------------------------
  # Triggered only when lg$M exceeds the cash floor, so the stabiliser
  # retains full capacity in downturns (repayment halts when lg$M is tight).
  # Repayment destroys money symmetrically: lg$M ↓, lg$L ↓ by equal amount.
  # bank_R = bank_NW + deposits - loans is unchanged (both sides shrink).
  lg_debt_repay <- 0
  if (lg$L > 0 && lg$M > LG_DEBT_FLOOR) {
    lg_debt_repay <- min(LG_AMORT_RATE * lg$L, lg$M - LG_DEBT_FLOOR)
    lg$M <- lg$M - lg_debt_repay
    lg$L <- lg$L - lg_debt_repay
  }

  # Fiscal EWMA update — smooths tax revenue for fiscal stabiliser.
  # Updated every period regardless of FISCAL_STABILIZER toggle.
  # alpha = 1/FISCAL_EWMA_TAU; warm-up from 0 during burn-in.
  {
    fiscal_alpha   <- 1 / FISCAL_EWMA_TAU
    tax_rev_smooth <- (1 - fiscal_alpha) * tax_rev_smooth + fiscal_alpha * lg_tax_now
  }
  
  
  hh_M <- pmax(hh_M, 0)
  
  # -�--�- L. RECORD -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
  om_now  <- f_status == "Open"; n_open <- sum(om_now)
  unemp_r <- 1 - sum(hh_Employed != 0L) / N_HH
  aw <- if (n_open > 0) mean(f_w[om_now]) else 0
  ap <- if (n_open > 0) mean(f_p[om_now]) else 0
  Ms_sorted <- sort(hh_M); n30 <- as.integer(N_HH * 0.3)
  
  hist[period, "unemp"]          <- unemp_r
  hist[period, "wage"]           <- aw
  hist[period, "price"]          <- ap
  hist[period, "gdp"]            <- sum(hh_C_nominal)
  hist[period, "hM"]             <- sum(hh_M)
  hist[period, "fM"]             <- sum(f_M[om_now])
  hist[period, "rw"]             <- aw / max(ap, 0.001)
  hist[period, "gini_hh"]        <- gini_fn(hh_M)
  hist[period, "total_loans"]    <- sum(f_L[om_now])
  hist[period, "n_open"]         <- n_open
  hist[period, "n_vacant"]       <- sum(f_status == "Vacant")
  hist[period, "n_openings"]     <- n_openings
  hist[period, "n_exits"]        <- n_exits
  hist[period, "lg_M"]           <- lg$M
  hist[period, "lg_L"]           <- lg$L
  hist[period, "lg_tax"]         <- lg$tax_rev
  hist[period, "lg_transfer"]    <- lg$transfer
  hist[period, "lg_deficit"]     <- lg$deficit
  hist[period, "lg_interest"]    <- lg$interest
  hist[period, "lg_debt_repay"]  <- lg_debt_repay
  hist[period, "g_share_current"]<- G_SHARE
  hist[period, "gov_transf_t"]   <- GOV_TRANSF
  hist[period, "owner_pay_total"]<- total_owner_pay
  hist[period, "m_dist_total"]   <- total_m_dist
  hist[period, "n_transfer_recip"]<- n_qual
  hist[period, "transfer_total"] <- transfer_bill
  hist[period, "mean_M_bot30"]   <- mean(Ms_sorted[seq_len(n30)])
  hist[period, "mean_M_top30"]   <- mean(Ms_sorted[(N_HH - n30 + 1L):N_HH])
  hist[period, "median_wage_t"]  <- median_wage_t
  hist[period, "income_thresh_t"]<- income_thresh
  hist[period, "supply"]         <- sum(f_Y[om_now])
  hist[period, "nom_demand"]     <- sum(hh_des_C)
  hist[period, "real_demand"]    <- sum(hh_des_C) / max(ap, 0.001)
  hist[period, "unmet_nom"]      <- sum(hh_des_C) - sum(hh_C_nominal)
  hist[period, "inventory_total"]<- sum(f_Inv[om_now]) + sum(f_Inv_wip[om_now])
  hist[period, "avg_mu"]         <- if (n_open > 0) mean(f_mu[om_now]) else 0
  hist[period, "n_hiring"]       <- sum(f_OpenPos[om_now] > 0L)
  y_des_tot <- if (n_open > 0) sum(f_N_des[om_now]) * mean(f_a[om_now]) else 1
  hist[period, "utilisation"]    <- sum(f_Sold[om_now]) / max(y_des_tot, 1)
  
  
  # Laspeyres price index
  if (!px_base_set) {
    px_base_q[om_now] <- f_Sold[om_now]
    px_base_p[om_now] <- ifelse(f_p[om_now] > 0, f_p[om_now], 1)
    px_base_set <- TRUE
  } else {
    new_base <- (px_base_q == 0) & om_now & (f_Sold > 0)
    if (any(new_base)) {
      px_base_q[new_base] <- f_Sold[new_base]
      px_base_p[new_base] <- ifelse(f_p[new_base] > 0, f_p[new_base], 1)
    }
  }
  valid_px <- (px_base_q > 0) & om_now
  price_index <- if (sum(valid_px) >= 2) {
    sum(px_base_q[valid_px] * f_p[valid_px]) /
      max(sum(px_base_q[valid_px] * px_base_p[valid_px]), 1e-9)
  } else px_prev_pi
  
  inflation <- (price_index / max(px_prev_pi, 1e-9) - 1) * 100
  px_prev_pi <- price_index

  # ----- GOV_TRANSF inflation correction (periodic, upward only) -----------
  # Fires every BF_INFL_FREQ periods. Computes cumulative inflation since
  # bf_last_pi. If positive, uprates by BF_INFL_PASSTHRU fraction.
  # Nominal transfer never falls below GOV_TRANSF_BASE (downward rigidity).
  # SFC note: higher GOV_TRANSF increases transfer_bill → lg$M ↓ → may
  # trigger lg$L ↑ in deficit periods; this is the correct fiscal response
  # to a real transfer uprating (government borrows to maintain real support).
  if (period %% BF_INFL_FREQ == 0L && period > BF_INFL_FREQ) {
    cumul_infl <- price_index / max(bf_last_pi, 1e-9) - 1
    if (cumul_infl > 0) {
      GOV_TRANSF <- max(GOV_TRANSF * (1 + BF_INFL_PASSTHRU * cumul_infl),
                        GOV_TRANSF_BASE)
      bf_last_pi <- price_index   # reset anchor to current price level
      cat(sprintf("  [BF correction t=%d] cumul_infl=%.2f%%  GOV_TRANSF → %.3f\n",
                  period, cumul_infl * 100, GOV_TRANSF))
    }
  }
  
  hist[period, "price_index"]     <- price_index
  hist[period, "inflation"]       <- inflation
  hist[period, "gdp_real"]        <- sum(hh_C_nominal) / max(price_index, 1e-9)
  hist[period, "real_wage_index"] <- if (price_index > 0 && aw > 0) aw / price_index else 0
  hist[period, "inflation_ma5"]   <- inflation   # smoothed post-loop
  hist[period, "inv_total"]       <- total_inv
  hist[period, "avg_a"]           <- if (n_open > 0) mean(f_a[om_now]) else 0
  hist[period, "avg_a_weighted"]  <- if (n_open > 0)
    sum(f_a[om_now] * f_N[om_now]) / max(sum(f_N[om_now]), 1) else 0
  hist[period, "n_hh"]            <- N_HH
  hist[period, "n_births"]        <- n_births
  hist[period, "lg_g_spend"]      <- lg_g_spend
  hist[period, "tax_rev_smooth"]  <- tax_rev_smooth
  hist[period, "g_n_firms_served"]<- g_n_firms_served
  {
    hh_C  <- sum(hh_C_nominal)
    total_demand_nom <- hh_C + lg_g_spend
    hist[period, "total_demand_nom"] <- total_demand_nom
    hist[period, "g_share_demand"]   <- if (total_demand_nom > 0)
      lg_g_spend / total_demand_nom else 0
  }
  # Bank diagnostics — endogenous money SFC consistency
  total_deposits_now <- sum(pmax(hh_M, 0)) + sum(pmax(f_M[om_now], 0)) + max(lg$M, 0)
  total_loans_now    <- sum(f_L) + lg$L
  bank_R_derived     <- bank_NW + total_deposits_now - total_loans_now
  total_M_now        <- sum(hh_M) + sum(f_M[om_now]) + lg$M
  # SFC check: total M = bank_NW + (total_loans - total_deposits) + total_deposits
  #           = bank_NW + total_loans ≈ bank_R + total_deposits
  # In endogenous money: total_M ≈ total_deposits (since f_M, hh_M, lg$M ARE deposits)
  # The balance sheet identity check: bank_R should be non-negative in normal operation
  sfc_residual <- bank_R_derived  # should be ≥ 0; large negative = SFC violation

  hist[period, "bank_NW"]            <- bank_NW
  hist[period, "bank_R"]             <- bank_R_derived
  hist[period, "bank_new_loans"]     <- bank_new_loans
  hist[period, "bank_dep_int_total"] <- bank_dep_int_total
  hist[period, "bank_loan_int_total"]<- bank_loan_int_total
  hist[period, "bank_profit"]        <- bank_profit_period
  hist[period, "bank_dividend"]      <- bank_dividend_period
  hist[period, "bank_bad_debt"]      <- bank_bad_debt_total
  hist[period, "bank_headroom"]      <- max(0, bank_NW / CAR_MIN - sum(f_L))  # CAR on private loans only
  hist[period, "total_M"]            <- total_M_now
  hist[period, "sfc_check"]          <- sfc_residual
  
  # Firm age diagnostics (open firms only)
  if (n_open > 0) {
    ages <- period - f_entry_period[om_now]
    ages <- ages[ages >= 0]   # guard against uninitialised slots
    hist[period, "mean_firm_age"]    <- mean(ages)
    hist[period, "med_firm_age"]     <- median(ages)
    hist[period, "share_young_firms"]<- mean(ages <= 12)  # fraction < 1 year old
  }
  
} # end main loop

cat("Simulation complete.\n")
BURN_IN <- 90L
RANGE   <- BURN_IN:MAXT
# Smooth inflation_ma5
hist$inflation_ma5 <- stats::filter(hist$inflation,
                                    rep(1/5, 5), sides = 2)
hist$inflation_ma5[is.na(hist$inflation_ma5)] <- hist$inflation[is.na(hist$inflation_ma5)]

# -�--�- TFP diagnostic (post burn-in) -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
hist$period <- seq_len(MAXT)
post <- hist[hist$period > BURN_IN, ]
lm_a      <- lm(avg_a ~ period, data = post)
a_slope   <- coef(lm_a)[["period"]]
a_resid   <- post$avg_a - predict(lm_a)
a_ac1     <- cor(post$avg_a[-1], post$avg_a[-nrow(post)])
cat("\n-�--�- TFP diagnostics (post-burn) -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-\n")
cat(sprintf("  mean a:        %.4f\n",  mean(post$avg_a)))
cat(sprintf("  std  a:        %.4f\n",  sd(post$avg_a)))
cat(sprintf("  trend slope:   %+.6f / period\n", a_slope))
cat(sprintf("  AC(a, lag=1):  %.3f\n",  a_ac1))
cat(sprintf("  detrended std: %.4f\n",  sd(a_resid)))
cat(sprintf("  target: slope -�� 0, detrended std -�� 0.015--�0.030\n"))

# -�--�- Summary -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
cat(sprintf("\n%-42s %10s\n", "Metric", "Value"))
cat(strrep("-", 54), "\n")
metrics <- list(
  "Mean unemployment (%)"         = sprintf("%.1f",  mean(hist$unemp[RANGE])*100),
  "Final real wage (w/p)"         = sprintf("%.3f",  tail(hist$rw, 1)),
  "Final price index"             = sprintf("%.3f",  tail(hist$price_index, 1)),
  "Mean inflation (%/period)"     = sprintf("%.2f",  mean(hist$inflation[RANGE])),
  "Cumulative inflation (%)"      = sprintf("%.1f",  (tail(hist$price_index,1)-1)*100),
  "Final HH Gini"                 = sprintf("%.3f",  tail(hist$gini_hh, 1)),
  "Final open firms"              = sprintf("%d",    as.integer(tail(hist$n_open, 1))),
  "Cumul. firm openings"          = sprintf("%d",    as.integer(sum(hist$n_openings[RANGE]))),
  "Cumul. firm exits"             = sprintf("%d",    as.integer(sum(hist$n_exits[RANGE]))),
  "Final LG debt"                 = sprintf("%.0f",  tail(hist$lg_L, 1)),
  "Final transfer recipients"     = sprintf("%d",    as.integer(tail(hist$n_transfer_recip,1))),
  "--- Banking ---"               = "",
  "Final bank NW"                 = sprintf("%.0f",  tail(hist$bank_NW, 1)),
  "Final bank reserves (derived)" = sprintf("%.0f",  tail(hist$bank_R, 1)),
  "Final total loans"             = sprintf("%.0f",  tail(hist$total_loans, 1)),
  "Final total M"                 = sprintf("%.0f",  tail(hist$total_M, 1)),
  "Mean deposit int / period"     = sprintf("%.1f",  mean(hist$bank_dep_int_total[RANGE])),
  "Mean loan int / period"        = sprintf("%.1f",  mean(hist$bank_loan_int_total[RANGE])),
  "Mean bank profit / period"     = sprintf("%.1f",  mean(hist$bank_profit[RANGE])),
  "Cumul. bad debt"               = sprintf("%.0f",  sum(hist$bank_bad_debt[RANGE])),
  "Final lending headroom"        = sprintf("%.0f",  tail(hist$bank_headroom, 1)),
  "Final SFC check (bank_R)"      = sprintf("%.2f",  tail(hist$sfc_check, 1))
)
for (nm in names(metrics))
  cat(sprintf("%-42s %10s\n", nm, metrics[[nm]]))


id <- format(Sys.time(), "%Y%m%d_%H%M%S")
# -�--�- Save CSV -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
write.csv(hist, paste(id, "abm_v5_tfp_cyclic_macro_pop.csv"), row.names = FALSE)
cat("\nMacro panel saved: abm_v5_tfp_cyclic_macro_pop.csv\n")

# -�--�- HH cross-section -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
hh_out <- data.frame(
  hh_id        = seq_len(N_HH),
  x            = round(hh_x, 2),
  y            = round(hh_y, 2),
  M            = round(hh_M, 4),
  L            = round(hh_L, 4),
  employed     = as.integer(hh_Employed != 0L),
  firm_id      = hh_Employed,
  wage         = round(hh_Wage, 4),
  firm_owner   = hh_firm_owner,
  unemp_spells = hh_unemp_spells,
  emp_spells   = hh_emp_spells
)
write.csv(hh_out, paste(id,"abm_v5_tfp_cyclic_hh_final_pop.csv"), row.names = FALSE)
cat("HH cross-section saved: abm_v5_tfp_cyclic_hh_final_pop.csv\n")

# -�--�- Basic plots -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
periods <- seq_len(MAXT)

png(paste(id,"abm_v5_tfp_cyclic_overview.png"), width = 1400, height = 1800, res = 120)
par(mfrow = c(6, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.6, 0))



plot(periods[RANGE], hist$unemp[RANGE]*100, type="l", col="red", lwd=2,
     main="Unemployment (%)", xlab="Period", ylab="%")
abline(h=mean(hist$unemp[RANGE]*100), lty=2, col='red')

plot(periods[RANGE], hist$wage[RANGE], type="l", col="steelblue", lwd=2,
     main="Wage & Price", xlab="Period", ylab="units",
     ylim=range(c(hist$wage[RANGE], hist$price[RANGE])))
lines(periods[RANGE], hist$price[RANGE], col="darkorange", lwd=2)
legend("topright", c("wage","price"), col=c("steelblue","darkorange"),
       lwd=2, cex=0.7, bty="n")

plot(periods[RANGE], hist$price_index[RANGE], type="l", col="darkorange", lwd=2,
     main="Price Index (Laspeyres)", xlab="Period", ylab="index")
abline(h=1, lty=3)

plot(periods[RANGE], hist$inflation[RANGE], type="h", col=ifelse(hist$inflation[RANGE]>=0,"salmon","lightblue"),
     main="Inflation (%/period)", xlab="Period", ylab="%")
lines(periods[RANGE], hist$inflation_ma5[RANGE], col="darkred", lwd=2)
abline(h=0, col="black", lwd=0.6)

plot(periods[RANGE], hist$gini_hh[RANGE], type="l", col="purple", lwd=2,
     main="HH Wealth Gini", xlab="Period", ylab="Gini", ylim=c(0, 1))

plot(periods[RANGE], hist$n_open[RANGE], type="l", col="dimgrey", lwd=2,
     main="Open Firms", xlab="Period", ylab="# firms")

plot(periods[RANGE], hist$gdp[RANGE], type="l", col="purple", lwd=2,
     main="GDP (nominal vs real)", xlab="Period", ylab="currency",
     ylim=range(c(hist$gdp[RANGE], hist$gdp_real[RANGE])))
lines(periods[RANGE], hist$gdp_real[RANGE], col="mediumblue", lwd=2, lty=2)
legend("topright", c("nominal","real"), col=c("purple","mediumblue"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

plot(periods[RANGE], hist$avg_a[RANGE], type="l", col="red", lwd=2,
     main="Mean TFP (pseudo-capital a)", xlab="Period", ylab="a level")
abline(h=1, lty=3)

plot(periods[RANGE], hist$inv_total[RANGE], type="l", col="steelblue", lwd=2,
     main="Investment (I -�� workers)", xlab="Period", ylab="currency units")

plot(periods[RANGE], hist$n_hh[RANGE], type="l", col="darkgreen", lwd=2,
     main="Population (N_HH)", xlab="Period", ylab="# households")
abline(h=hist$n_hh[BURN_IN], lty=2, col="darkgreen")

plot(periods[RANGE], hist$n_births[RANGE], type="h", col="seagreen",
     main="Births per period", xlab="Period", ylab="# new HH")

plot(periods[RANGE], hist$unemp[RANGE]*hist$n_hh[RANGE], type="l", col="red", lwd=2,
     main="Unemployed (headcount)", xlab="Period", ylab="# persons")

# -�--�- Firm age panels -�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�--�-
plot(periods[RANGE], hist$mean_firm_age[RANGE], type="l", col="chocolate", lwd=2,
     main="Mean Firm Age (periods)", xlab="Period", ylab="periods")
lines(periods[RANGE], hist$med_firm_age[RANGE], col="chocolate4", lwd=2, lty=2)
legend("topleft", c("mean","median"), col=c("chocolate","chocolate4"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

plot(periods[RANGE], hist$share_young_firms[RANGE]*100, type="l", col="coral", lwd=2,
     main="Young Firms (age -�� 12, %)", xlab="Period", ylab="%",
     ylim=c(0, 100))
abline(h=50, lty=3, col="grey60")

# Gov. Procurement vs BF transfers
plot(periods[RANGE], hist$lg_g_spend[RANGE], type="l", col="navy", lwd=2,
     main="Gov. Procurement (G) vs Transfers", xlab="Period", ylab="currency units",
     ylim=range(c(hist$lg_g_spend[RANGE], hist$transfer_total[RANGE])))
lines(periods[RANGE], hist$transfer_total[RANGE], col="steelblue", lwd=2, lty=2)
legend("topleft", c("G procurement","BF transfers"), col=c("navy","steelblue"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

# G vs HH consumption decomposition
ylim_d <- range(c(hist$gdp[RANGE], hist$lg_g_spend[RANGE], hist$total_demand_nom[RANGE]))
plot(periods[RANGE], hist$total_demand_nom[RANGE], type="l", col="black", lwd=2,
     main="Demand decomposition (G vs HH-C)", xlab="Period", ylab="currency units",
     ylim=ylim_d)
lines(periods[RANGE], hist$gdp[RANGE],       col="purple",   lwd=2, lty=1)
lines(periods[RANGE], hist$lg_g_spend[RANGE],col="navy",     lwd=2, lty=2)
legend("topleft", c("Total demand","HH consumption (C)","Gov procurement (G)"),
       col=c("black","purple","navy"), lwd=2, lty=c(1,1,2), cex=0.65, bty="n")

# G share of total demand + firms served
par(mar=c(3,3,2,4))
plot(periods[RANGE], hist$g_share_demand[RANGE]*100, type="l", col="navy", lwd=2,
     main="G share of demand & firms served", xlab="Period", ylab="G / (G+C) %")
abline(h=mean(hist$g_share_demand[RANGE]*100), lty=2, col="navy")
par(new=TRUE)
plot(periods[RANGE], hist$g_n_firms_served[RANGE], type="l", col="coral", lwd=1.5,
     axes=FALSE, xlab="", ylab="")
axis(4, col.axis="coral"); mtext("# firms served", side=4, line=2.5, col="coral", cex=0.7)
legend("topright", c("G/demand %","Firms served"), col=c("navy","coral"),
       lwd=c(2,1.5), cex=0.65, bty="n")
par(mar=c(3,3,2,1))

# Firm cash stock
plot(periods[RANGE], hist$fM[RANGE], type="l", col="darkred", lwd=2,
     main="Firm aggregate cash (fM)", xlab="Period", ylab="currency units")
abline(h=0, lty=3)

dev.off()
cat("Overview plot saved: abm_v5_stable2_overview_pop.png\n")

# ----- Banking diagnostics plot ------------------------------------------
png(paste(id,"abm_v5_banking_diagnostics.png"), width = 1400, height = 1200, res = 120)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.6, 0))

# 1. Bank NW (equity)
plot(periods[RANGE], hist$bank_NW[RANGE], type="l", col="darkgreen", lwd=2,
     main="Bank Net Worth (equity)", xlab="Period", ylab="currency units")
abline(h=0, lty=3, col="red")

# 2. Bank reserves (derived)
plot(periods[RANGE], hist$bank_R[RANGE], type="l", col="steelblue", lwd=2,
     main="Bank Reserves (derived)", xlab="Period", ylab="currency units")
abline(h=0, lty=3, col="red")

# 3. Total loans outstanding
plot(periods[RANGE], hist$total_loans[RANGE], type="l", col="darkred", lwd=2,
     main="Total Loans Outstanding", xlab="Period", ylab="currency units")

# 4. Lending headroom (CAR constraint)
plot(periods[RANGE], hist$bank_headroom[RANGE], type="l", col="darkorange", lwd=2,
     main="Lending Headroom (CAR)", xlab="Period", ylab="currency units")
abline(h=0, lty=3, col="red")

# 5. Interest flows: deposit vs loan
ylim_int <- range(c(hist$bank_dep_int_total[RANGE], hist$bank_loan_int_total[RANGE]))
plot(periods[RANGE], hist$bank_loan_int_total[RANGE], type="l", col="darkred", lwd=2,
     main="Interest Flows", xlab="Period", ylab="currency units", ylim=ylim_int)
lines(periods[RANGE], hist$bank_dep_int_total[RANGE], col="steelblue", lwd=2, lty=2)
legend("topright", c("Loan int (received)","Deposit int (paid)"),
       col=c("darkred","steelblue"), lwd=2, lty=c(1,2), cex=0.7, bty="n")

# 6. Bank profit & dividends
ylim_pd <- range(c(hist$bank_profit[RANGE], hist$bank_dividend[RANGE]))
plot(periods[RANGE], hist$bank_profit[RANGE], type="l", col="darkgreen", lwd=2,
     main="Bank Profit & Dividends", xlab="Period", ylab="currency units", ylim=ylim_pd)
lines(periods[RANGE], hist$bank_dividend[RANGE], col="gold3", lwd=2, lty=2)
abline(h=0, lty=3)
legend("topright", c("Net profit","Dividends"), col=c("darkgreen","gold3"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

# 7. Bad debt per period
plot(periods[RANGE], hist$bank_bad_debt[RANGE], type="h", col="red",
     main="Bad Debt (per period)", xlab="Period", ylab="currency units")

# 8. Total M in the economy
plot(periods[RANGE], hist$total_M[RANGE], type="l", col="purple", lwd=2,
     main="Total Money Supply (M)", xlab="Period", ylab="currency units")

# 9. SFC consistency check (bank_R)
plot(periods[RANGE], hist$sfc_check[RANGE], type="l", col="black", lwd=2,
     main="SFC Check (bank_R, should be >= 0)", xlab="Period", ylab="currency units")
abline(h=0, lty=3, col="red")

dev.off()
cat("Banking diagnostics plot saved.\n")


