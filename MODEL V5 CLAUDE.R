# =============================================================================
# ABM v5 stable2 — R port
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

set.seed(42)

# ── Parameters ────────────────────────────────────────────────
N_HH          <- 2000
MAXT          <- 300

# ── Population growth ─────────────────────────────────────────
# Demographic engine: each period, ceil(N_HH * POP_GROWTH_RATE) new households
# are born inside the municipality boundary, inheriting a small wealth endowment
# drawn from the lower tail of the current HH wealth distribution (mimicking
# entry with modest initial assets, consistent with lifecycle models;
# see Huggett 1996, Aiyagari 1994).
# New entrants are unemployed and have zero firm ownership.
# N_HH is updated dynamically; all HH arrays grow accordingly.
# Set POP_GROWTH_RATE <- 0 to disable without changing anything else.
GOV_TRANSF    <- 5
POP_GROWTH_RATE  <- 0.002        # ~0.2% per period; adjust to target annual rate
POP_WEALTH_FRAC  <- 0.25         # new entrant wealth = POP_WEALTH_FRAC * median(hh_M)
POP_WEALTH_FLOOR <- GOV_TRANSF         # absolute minimum endowment for new entrant
MIN_WAGE      <- 2
AS_SCALE       <- 1.0   # animal spirits scale: floor = AS_SCALE * mean(wage)
ANIMAL_SPIRITS_MIN <- MIN_WAGE  # hard floor of last resort if wages collapse
THETA         <- 0.6
ALPHA_S       <- 0.35   # was 0.6 — slower sales expectation update damps demand-signal amplification
LAMBDA_N      <- 0.50   # was 0.95 — gradual hiring spreads entry waves over ~4 periods
# instead of snapping to target immediately; key fix for oscillations
MU_PHI        <- 0.02
MU_MIN        <- 0.01
MU_MAX        <- 1000
PHI_W         <- 0.15   # was 0.25 — slower wage adjustment reduces wage-price feedback

MAX_WAGE      <- 1000
LAMBDA0       <- 1e-3
P_FLOOR       <- 0.1
GAMMA_MD      <- 0.03     # was 0.05 — slower markdown, reduces deflationary spikes
A1            <- 0.15
A2            <- 0.85
RL_FIRM       <- 0.015
TFP_GROWTH    <- 1.00
SIGMA_TFP     <- 0.03   # used only if macro TFP shock block is uncommented

PAYROLL_TAX   <- 0.20
PROFIT_TAX    <- 0.27
WEALTH_FLOOR  <- 10        # overridden dynamically each period
INCOME_FRAC   <- 0.8
RL_GOV        <- 0.005
BF_MEDIAN_FRAC  <- 0.45    # transfer = 45% of median employed wage
BF_WEALTH_FRAC  <- 6.00    # wealth floor = 6x median wage
# ── Smooth fiscal rule ────────────────────────────────────────
# Government spends G_SHARE of net-of-transfer tax revenue on
# procurement every period, unconditionally. No threshold, no cap.
# Funded from lg$M first; deficit-financed if deposits run short.
# LG_M is a pure state variable — its trajectory is an outcome,
# not a control parameter. Remove LG_M_CAP, LG_G_FRAC, G_TAX_FRAC.
G_SHARE         <- 0.40    # G = G_SHARE * (tax_rev - transfer_bill)
# 0.40 gives ~30% of GDP in procurement at
# current tax/GDP ratios; adjust to calibrate

STARTUP_CAPITAL   <- 75
SC_WEALTH_FRAC    <- 0.50   # fallback SC = fraction of mean(hh_M); keeps entry
# expensive when BF transfers inflate deposits
INTERCEPT_OPEN    <- -5
BETA_M            <-  0.03
BETA_UNEMP        <-  1.8
BETA_GAIN         <-  0.10
BETA_USPELLS      <-  0.15
BETA_ESPELLS      <- -0.04
USPELLS_CAP       <- 20L    # logit saturation: spells beyond 20 add no further push
ENTRY_RATE_CAP    <- 0.01   # max new firms per period as fraction of N_HH (~1%)
ENTRY_COOLDOWN    <- 8L     # periods of reduced entry after a mass-entry episode;
# after n_openings >= ENTRY_RATE_CAP*N_HH, cap halves
# for ENTRY_COOLDOWN periods to dampen the step-function
# transition from trap to full-employment
BETA_SAT          <- -0.05
SIGMA_W           <-  0.10
SAT_RADIUS        <-  2000.0
OWNER_PROFIT_SHARE <- 0.15
M_DIST_RATE       <- 0.2
M_DIST_PERIOD     <- 12L

INV_RATE       <- 0.02
A_DEPRECIATION <- 0.006
A_MIN          <- 0.1

# ── Wage–productivity linkage (heterodox, non-neoclassical) ───

wageperiod <- 1
priceperiod <- 1

# Mechanism 1: Efficiency wage floor (Akerlof-Yellen fair wage)
#   Firms with above-average TFP defend a higher wage floor.
#   f_w >= MIN_WAGE * (f_a / mean_a)^ETA_W
#   ETA_W = 0 -> no linkage; ETA_W = 1 -> full proportional indexation.
#   Empirical target: 0.4-0.6 (Carlsson, Messina & Skans 2016).
#   Set USE_EFF_WAGE <- FALSE to disable without changing parameters.
ETA_W        <- 0.6
USE_EFF_WAGE <- TRUE   # OFF for now; set TRUE to activate


# Mechanism 2: Kaleckian surplus sharing (Bhaduri & Marglin 1990)
#   Workers capture THETA_W share of each period's TFP *gain* as a wage rise.
#   Only positive TFP gains are passed through (downward nominal rigidity).
#   THETA_W = 0 -> no linkage; THETA_W = 1 -> full pass-through of gains.
#   Empirical target: 0.2-0.4 (labour share dynamics literature).
THETA_W      <- 0.4
USE_SURPLUS  <- TRUE    # ON

QUIT_FRICTION  <- 10.0
QUIT_SIGMA     <- 1.5
SWITCH_SIGMA   <- 0.5
ACCEPT_SIGMA   <- 0.1
COST_PER_M     <- 1e-9
MAX_CONSIDERED <- 30L
MAX_VISITS     <- 100L

# ── Utilities ─────────────────────────────────────────────────
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

# ── Geography: Ilhéus boundary from geobr CSV ─────────────────
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

# ── Road network definition ────────────────────────────────────
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

# ── Agent arrays ──────────────────────────────────────────────
# N_SLOTS grows with births: each new HH gets a Vacant firm slot.
# N_SLOTS_0 records the initial count; used only where a frozen reference
# is genuinely needed (none currently).
N_SLOTS <- N_HH
N_SLOTS_0 <- N_SLOTS

hh_M          <- pmax(rnorm(N_HH, 150, 50), 0)
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
  f_Inv[fi]    <<- 0;  f_Sold[fi]  <<- 0;  f_Rev[fi]    <<- 0
  f_WB[fi]     <<- 0;  f_payroll_tax[fi] <<- 0; f_Tax[fi] <<- 0
  f_Profits_gross[fi] <<- 0; f_Profits[fi] <<- 0
  f_owner_pay[fi] <<- 0; f_m_dist[fi] <<- 0
  f_Int[fi]    <<- 0;  f_L[fi]     <<- 0
  f_inv_spend[fi] <<- 0; f_Profits_prev[fi] <<- 0
  f_p[fi]  <<- max(rnorm(1, 5.5, 2.0), P_FLOOR)
  f_mu[fi] <<- 0.35
  f_w[fi]  <<- max(rnorm(1, 10.0, 3.0), MIN_WAGE)
  f_a[fi]  <<- 1.0
  f_a_prev[fi] <<- 1.0
  f_OpenPos[fi] <<- 0L
  f_entry_period[fi] <<- 0L   # cleared on reset; set to current period at open
}

# ── Distance matrix (pre-allocated to avoid O(N²) copies each birth period) ─
# N_MAX is an upper bound on population after MAXT periods of growth.
# All writes during birth events fill in-place; no matrix copy ever occurs.
cat("Building distance matrix...\n")
N_MAX  <- ceiling(N_HH * (1 + POP_GROWTH_RATE)^MAXT * 1.05)
d_hf   <- matrix(0.0, N_MAX, N_MAX)
coords <- cbind(hh_x, hh_y)
d_hf[seq_len(N_HH), seq_len(N_HH)] <- as.matrix(dist(coords))
cat("Done.\n")

# ── Period-1 seed ─────────────────────────────────────────────
N_SEED    <- N_HH %/% 50L   # 2% seed rate; 10% caused catastrophic over-entry

openfirms<-f_status=="Open"
# if(length(openfirms) >0){
#   STARTUP_CAPITAL<-mean(f_w[f_w>0])*1.2
# }else{
STARTUP_CAPITAL<-abs(mean(hh_M))*1.3
# }

eligible<-which(hh_M>STARTUP_CAPITAL)

N_SEED<-pmin(length(eligible),N_SEED)

seed_slots <- sample(eligible, N_SEED, replace = FALSE)


for (fi in seed_slots) {
  
  hh_M[fi]         <- hh_M[fi] - STARTUP_CAPITAL
  f_M[fi]          <- STARTUP_CAPITAL
  f_L_max[fi]      <- STARTUP_CAPITAL * 3
  reset_firm(fi)
  f_status[fi]     <- "Open"
  hh_firm_owner[fi]<- fi
  f_entry_period[fi] <- 1L
}
cat(sprintf("Seeded %d firms at period 1.\n", N_SEED))

# ── Government ────────────────────────────────────────────────
lg <- list(M = 500, L = 0, tax_rev = 0, transfer = 0,
           deficit = 0, interest = 0)
lg_g_spend <- 0
entry_cooldown_t <- 0L   # counts down after a mass-entry episode

# ── History ───────────────────────────────────────────────────
hist_keys <- c("unemp","wage","price","gdp","hM","fM","rw","gini_hh",
               "total_loans","n_open","n_vacant","n_openings","n_exits",
               "lg_M","lg_L","lg_tax","lg_transfer","lg_deficit","lg_interest",
               "owner_pay_total","m_dist_total","n_transfer_recip","transfer_total",
               "mean_M_bot30","mean_M_top30","median_wage_t","income_thresh_t",
               "supply","nom_demand","real_demand","unmet_nom",
               "inventory_total","avg_mu","n_hiring","utilisation",
               "price_index","inflation","inflation_ma5","gdp_real",
               "real_wage_index","inv_total","avg_a","avg_a_weighted",
               "n_hh","n_births","lg_g_spend",
               "mean_firm_age","med_firm_age","share_young_firms")
hist <- as.data.frame(matrix(0, nrow = MAXT, ncol = length(hist_keys)))
names(hist) <- hist_keys

# ── Price index state ─────────────────────────────────────────
px_base_q   <- numeric(N_SLOTS)
px_base_p   <- numeric(N_SLOTS)
px_base_set <- FALSE
px_prev_pi  <- 1.0

# ═════════════════════════════════════════════════════════════
# MAIN LOOP
# ═════════════════════════════════════════════════════════════
for (period in seq_len(MAXT)) {
  
  if (period %% 1 == 0) {
    n_om  <- sum(f_status == "Open")
    u_now <- 1 - sum(hh_Employed != 0L) / N_HH
    cat(sprintf("t=%3d | open=%3d | u=%.1f%% | HH_M=%.0f | LG_M=%.0f LG_L=%.0f | STARTUP CAPITAL= %.0f | Mean M = %0.f | Eligible Entrep. = %0.f\n",
                period, n_om, u_now*100, sum(hh_M), lg$M, lg$L,STARTUP_CAPITAL,mean(hh_M),sum(hh_M>STARTUP_CAPITAL)))
    
    
    #print(STARTUP_CAPITAL)
    #print(mean(hh_M))
    #print(sum(hh_M>STARTUP_CAPITAL))
  }
  
  # ── A. CLOSURE ──────────────────────────────────────────────
  n_exits <- 0L
  if (period > 1L) {
    om_mask  <- f_status == "Open"
    insolv   <- om_mask & ((f_M + f_Inv * f_p <= 0) | (f_L > f_L_max))
    ins_idx  <- which(insolv)
    for (fi in ins_idx) {
      deficit <- max(0, -(f_M[fi] + f_Inv[fi] * f_p[fi]))
      cover   <- min(hh_M[fi], deficit)
      hh_M[fi] <- hh_M[fi] - cover
      residual <- deficit - cover
      if (residual > 0) hh_L[fi] <- hh_L[fi] + residual
      hh_firm_owner[fi] <- 0L
      f_status[fi]      <- "Vacant"
      f_L_max[fi]       <- 0
    }
    n_exits <- length(ins_idx)
  }
  
  # Zero non-open fields
  not_om <- f_status != "Open"
  for (arr_name in c("f_M","f_L","f_L_max","f_N","f_N_des","f_Y","f_S_hat",
                     "f_I_star","f_Inv","f_Sold","f_Rev","f_WB",
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
  
  # ── A2. POPULATION GROWTH ───────────────────────────────────
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
    
    # Wealth endowment
    med_M  <- max(median(hh_M), 1e-3)
    new_M  <- pmax(rnorm(n_births, POP_WEALTH_FRAC * med_M,
                         0.1 * POP_WEALTH_FRAC * med_M), POP_WEALTH_FLOOR)
    
    # ── Distance matrix: fill new rows/cols in-place (no matrix copy) ───
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
    
    # ── HH arrays ────────────────────────────────────────────────
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
    
    # ── Firm arrays (Vacant slots for new citizens) ───────────────
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
    f_p             <- c(f_p,             numeric(n_births))
    f_mu            <- c(f_mu,            rep(0.25, n_births))
    f_w             <- c(f_w,             rep(MIN_WAGE, n_births))
    f_a             <- c(f_a,             rep(1.0, n_births))
    f_a_prev        <- c(f_a_prev,        rep(1.0, n_births))
    f_OpenPos       <- c(f_OpenPos,       integer(n_births))
    f_entry_period  <- c(f_entry_period,  integer(n_births))
    
    # ── Laspeyres price index state vectors ───────────────────────
    # px_base_q and px_base_p must stay co-indexed with firm arrays.
    # New slots have no base-period observation yet (zero = not yet set).
    px_base_q <- c(px_base_q, numeric(n_births))
    px_base_p <- c(px_base_p, numeric(n_births))
    
    # ── Dimension counters ────────────────────────────────────────
    N_HH   <- N_HH_old + n_births
    N_SLOTS <- N_HH
  }
  
  # ── B. LOGISTIC ENTRY ───────────────────────────────────────
  # All HH indices up to N_SLOTS (grows with births) are eligible.
  # New citizens have a Vacant firm slot allocated at birth.
  
  
  # Dynamic startup capital — EWMA-smoothed to prevent wild oscillation.
  # Target = mean wage bill per firm; falls back to 3x MIN_WAGE when
  # no employment exists yet. Alpha=0.1 means ~10-period adjustment lag.
  {
    open_wages  <- f_w[f_status == "Open" & f_w > 0]
    open_N      <- f_N[f_status == "Open" & f_N > 0]
    mean_hh_M   <- max(mean(hh_M), MIN_WAGE, na.rm = TRUE)
    SC_fallback <- max(mean_hh_M * SC_WEALTH_FRAC, MIN_WAGE * 2)
    SC_target   <- if (length(open_wages) > 0 && length(open_N) > 0)
      mean(open_wages) * mean(open_N)
    else
      SC_fallback
    SC_target <- max(SC_target, SC_fallback, na.rm = TRUE)
    if (!is.finite(SC_target)) SC_target <- SC_fallback
    STARTUP_CAPITAL <- 0.8 * STARTUP_CAPITAL + 0.2 * SC_target
  }
  
  # ── GOV_TRANSF update (before labour market so quit decisions use
  #    current-period transfer as reservation wage) ──────────────
  {
    emp_wages_pre  <- hh_Wage[hh_Employed != 0L]
    median_wage_pre <- if (length(emp_wages_pre) > 0) median(emp_wages_pre) else MIN_WAGE
    GOV_TRANSF  <- max(BF_MEDIAN_FRAC * median_wage_pre, MIN_WAGE)
    WEALTH_FLOOR <- max(BF_WEALTH_FRAC * median_wage_pre, GOV_TRANSF * 6)
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
      hh_M[fi]          <- hh_M[fi] - STARTUP_CAPITAL
      f_M[fi]           <- STARTUP_CAPITAL
      f_L_max[fi]       <- STARTUP_CAPITAL * 3
      reset_firm(fi)
      f_status[fi]      <- "Open"
      hh_firm_owner[fi] <- fi
      f_entry_period[fi] <- period
      hh_unemp_spells[fi] <- 0L
      hh_emp_spells[fi]   <- 0L
      n_openings <- n_openings + 1L
    }
  }
  
  om      <- f_status == "Open"
  openidx <- which(om)
  
  
  # ── ENDOGENOUS ANIMAL SPIRITS ────────────────────────────────
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
  
  
  # ── C. PLANNING ─────────────────────────────────────────────
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
  
  # ── D. LABOUR MARKET ────────────────────────────────────────
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
  
  # ── E. LAYOFFS ───────────────────────────────────────────────
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
  
  # ── F. PRODUCTION & WAGES ───────────────────────────────────
  f_Y[] <- 0; f_WB[] <- 0; f_Sold[] <- 0
  # ── TFP cyclical growth ──────────────────────────────────────
  # Common macro shock (70% of variance) creates correlated cycles across
  # all open firms; firm-level idiosyncratic shock (30%) adds dispersion.
  # With TFP_GROWTH=1.000 the multiplicative term is neutral and the
  # pseudo-capital investment block below drives any slow upward drift.
  # Raise TFP_GROWTH to 1.002–1.005 for an explicit positive trend.
  # if (any(om)) {
  #   macro_eps <- rnorm(1L,          mean = 0, sd = SIGMA_TFP * 0.7)
  #   firm_eps  <- rnorm(sum(om), mean = 0, sd = SIGMA_TFP * 0.3)
  #   f_a[om]   <- pmax(f_a[om] * TFP_GROWTH * exp(macro_eps + firm_eps), A_MIN)
  # }
  
  
  # Pseudo-capital investment — diminishing returns to TFP accumulation
  # Gain = I / (denom * f_a[fi]) so that high-productivity firms extract
  # progressively less TFP per unit of investment (Romer 1990).
  # Steady-state: a* = sqrt(I / (denom * A_DEPRECIATION))
  f_inv_spend[] <- 0; total_inv <- 0
  for (fi in openidx) {
    f_a[fi] <- max(f_a[fi] * (1 - A_DEPRECIATION), A_MIN)
    prev_profit <- f_Profits_prev[fi]
    if (prev_profit <= 0) next
    I <- INV_RATE * prev_profit
    I <- min(I, max(0, f_M[fi]))
    if (I <= 0) next
    f_M[fi]        <- f_M[fi] - I
    f_inv_spend[fi] <- I
    total_inv      <- total_inv + I
    wks <- which(hh_Employed == fi)
    if (length(wks) > 0) hh_M[wks] <- hh_M[wks] + I / length(wks)
    denom   <- f_p[fi] * max(f_N[fi], 1) + 1e-9
    f_a[fi] <- f_a[fi] + I / (denom * f_a[fi])   # <-- only change
  }
  
  Y_cap   <- f_a * f_N
  Y_plan  <- pmax(f_S_hat + (f_I_star - f_Inv), 0)
  f_Y[om] <- pmin(Y_plan[om], Y_cap[om])
  f_Inv[om] <- f_Inv[om] + f_Y[om]
  
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
  
  # ── G. PRICING ───────────────────────────────────────────────
  pr <- which(om & f_Y > 0)
  sr <- which(om & f_Y == 0 & f_Inv > 0)
  
  
  if(period %% priceperiod == 0){
    
    if (length(pr) > 0) f_p[pr] <- (1 + f_mu[pr]) * (f_w[pr] / f_a[pr])
    if (length(sr) > 0) {
      sig_sr <- pmin(f_Inv[sr] / (f_S_hat[sr] + 1e-5), 5)
      f_p[sr] <- pmax(f_p[sr] * (1 - GAMMA_MD * sig_sr), P_FLOOR)
    }
  }
  
  # ── H. FIRM CREDIT ───────────────────────────────────────────
  neg_idx <- which(om & f_M < 0)
  for (fi in neg_idx) {
    f_L[fi] <- f_L[fi] - f_M[fi]; f_M[fi] <- 0
  }
  f_Int[om] <- f_L[om] * RL_FIRM
  f_M[om]   <- f_M[om] - f_Int[om]
  
  # ── I. GOODS MARKET ──────────────────────────────────────────
  # Defensive scrub: any NA/NaN in inventory or price would crash the
  # inner loop. These should not occur after the animal_spirits fix,
  # but guard anyway for robustness.
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
  
  # ── J. ACCOUNTING ────────────────────────────────────────────
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
  
  f_Profits_prev[openidx] <- f_Profits[openidx]
  f_Profits_prev[!om]     <- 0
  
  # ── WAGE–PRODUCTIVITY LINKAGE ────────────────────────────────
  # Runs after pseudo-capital so wages see the updated f_a values.
  # Both mechanisms are independently switchable via USE_EFF_WAGE / USE_SURPLUS.
  
  if (any(om)) {
    
    # -- Mechanism 1: Efficiency wage floor (Akerlof-Yellen fair wage) --------
    # OFF by default (USE_EFF_WAGE <- FALSE in parameters).
    # Firms with above-average TFP defend a higher wage floor:
    #   w_floor_i = MIN_WAGE * (a_i / mean_a)^ETA_W
    # Workers at more productive firms resist wage cuts below this reference.
    # Activate by setting USE_EFF_WAGE <- TRUE.
    if (USE_EFF_WAGE) {
      mean_a_now   <- mean(f_a[om])
      prod_floor   <- MIN_WAGE * (f_a[om] / (mean_a_now + 1e-9)) ^ ETA_W
      f_w[om]      <- pmax(f_w[om], prod_floor)
    }
    
    # -- Mechanism 2: Kaleckian surplus sharing (Bhaduri & Marglin 1990) ------
    # ON by default (USE_SURPLUS <- TRUE in parameters).
    # Workers capture THETA_W of each firm's proportional TFP *gain* this period.
    # Only positive gains pass through (downward nominal wage rigidity):
    #   delta_a_i = max(a_i_t - a_i_{t-1}, 0)
    #   w_i <- w_i + THETA_W * delta_a_i * w_i
    # When TFP falls, wages are sticky: no symmetric cut is applied.
    if (USE_SURPLUS) {
      delta_a      <- pmax(f_a[om] - f_a_prev[om], 0)
      f_w[om]      <- f_w[om] + THETA_W * delta_a * f_w[om]
      f_w[om]      <- clamp(f_w[om], MIN_WAGE, MAX_WAGE)
    }
    
    # Store this period's TFP for next period's surplus calculation
    f_a_prev[om]  <- f_a[om]
    f_a_prev[!om] <- 0
  }
  # ─────────────────────────────────────────────────────────────
  
  # Firm cash-flow
  for (fi in openidx)
    f_M[fi] <- f_M[fi] + (f_Profits[fi] - f_owner_pay[fi]) + f_WB[fi]
  
  # Loan repayment
  repay_idx <- which(om & f_Profits > 0 & f_L > 0)
  for (fi in repay_idx) {
    repay   <- min(f_Profits[fi] * 0.10, f_L[fi])
    f_L[fi] <- f_L[fi] - repay; f_M[fi] <- f_M[fi] - repay
  }
  
  # Bank interest → HH
  firm_int_total <- sum(f_Int[om])
  if (firm_int_total > 0) hh_M <- hh_M + firm_int_total / N_HH
  
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
  
  
  
  # ── K. BOLSA FAMÍLIA ─────────────────────────────────────────
  # Transfer and eligibility both indexed to median employed wage,
  # so real value is maintained as the economy grows or contracts.
  # GOV_TRANSF was already updated before the labour market above.
  # WEALTH_FLOOR = BF_WEALTH_FRAC * median_wage ensures broad coverage
  # that doesn't collapse as the lower tail of the wealth distribution
  # deteriorates. The GOV_TRANSF * 6 floor prevents the perverse case
  # where a single transfer payment briefly lifts someone above the
  # threshold and disqualifies them the next period.
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
  lg_int <- lg$L * RL_GOV
  if (lg$M >= lg_int) lg$M <- lg$M - lg_int else { lg$L <- lg$L + lg_int - lg$M; lg$M <- 0 }
  lg$interest <- lg_int
  lg$transfer <- transfer_bill
  # lg$deficit computed after procurement so lg_g_spend is included
  
  
  # ── Government procurement (smooth fiscal rule) ──────────────
  # G = G_SHARE * net fiscal revenue, every period, no threshold.
  # Net revenue = tax_rev - transfer_bill (what remains after BF).
  # Positive net revenue: government is in surplus → spend G_SHARE of it.
  # Zero or negative net revenue: no procurement (avoids pro-cyclical
  # borrowing to spend when the government is already running a deficit
  # purely on transfers). Interest payments already deducted above.
  # Funded from lg$M; deficit-financed only if lg$M runs short.
  # This produces a smooth, proportional automatic stabiliser:
  # procurement rises in booms (high taxes, low transfers) and falls
  # in recessions, but never switches off abruptly at a threshold.
  lg_g_spend <- 0
  net_fiscal  <- lg_tax_now - transfer_bill   # net of BF obligations
  if (net_fiscal > 0 && length(openidx) > 0) {
    gov_budget  <- net_fiscal * G_SHARE
    visit_order <- sample(openidx)
    for (j in visit_order) {
      if (gov_budget <= 0.01) break
      sq  <- max(0, f_Inv[j])
      p_j <- f_p[j]
      if (sq <= 0 || !is.finite(p_j) || p_j <= 0.001) next
      buy_q <- min(sq, gov_budget / p_j)
      pay   <- p_j * buy_q
      if (pay <= 0) next
      lg$M       <- lg$M - pay
      gov_budget  <- gov_budget - pay
      lg_g_spend  <- lg_g_spend + pay
      f_Inv[j]   <- max(0, f_Inv[j] - buy_q)
      f_Sold[j]  <- f_Sold[j] + buy_q
      f_M[j]     <- f_M[j] + pay
    }
    if (lg$M < 0) { lg$L <- lg$L - lg$M; lg$M <- 0 }
  }
  lg$deficit  <- transfer_bill + lg_int + lg_g_spend - lg_tax_now
  
  
  hh_M <- pmax(hh_M, 0)
  
  # ── L. RECORD ────────────────────────────────────────────────
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
  hist[period, "inventory_total"]<- sum(f_Inv[om_now])
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

# ── TFP diagnostic (post burn-in) ────────────────────────────
hist$period <- seq_len(MAXT)
post <- hist[hist$period > BURN_IN, ]
lm_a      <- lm(avg_a ~ period, data = post)
a_slope   <- coef(lm_a)[["period"]]
a_resid   <- post$avg_a - predict(lm_a)
a_ac1     <- cor(post$avg_a[-1], post$avg_a[-nrow(post)])
cat("\n── TFP diagnostics (post-burn) ──────────────────────────\n")
cat(sprintf("  mean a:        %.4f\n",  mean(post$avg_a)))
cat(sprintf("  std  a:        %.4f\n",  sd(post$avg_a)))
cat(sprintf("  trend slope:   %+.6f / period\n", a_slope))
cat(sprintf("  AC(a, lag=1):  %.3f\n",  a_ac1))
cat(sprintf("  detrended std: %.4f\n",  sd(a_resid)))
cat(sprintf("  target: slope ≈ 0, detrended std ≈ 0.015–0.030\n"))

# ── Summary ───────────────────────────────────────────────────
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
  "Final transfer recipients"     = sprintf("%d",    as.integer(tail(hist$n_transfer_recip,1)))
)
for (nm in names(metrics))
  cat(sprintf("%-42s %10s\n", nm, metrics[[nm]]))


id <- format(Sys.time(), "%Y%m%d_%H%M%S")
# ── Save CSV ──────────────────────────────────────────────────
write.csv(hist, paste(id, "abm_v5_tfp_cyclic_macro_pop.csv"), row.names = FALSE)
cat("\nMacro panel saved: abm_v5_tfp_cyclic_macro_pop.csv\n")

# ── HH cross-section ─────────────────────────────────────────
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

# ── Basic plots ───────────────────────────────────────────────
periods <- seq_len(MAXT)

png(paste(id,"abm_v5_tfp_cyclic_overview.png"), width = 1400, height = 1500, res = 120)
par(mfrow = c(5, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.6, 0))



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
     main="Investment (I → workers)", xlab="Period", ylab="currency units")

plot(periods[RANGE], hist$n_hh[RANGE], type="l", col="darkgreen", lwd=2,
     main="Population (N_HH)", xlab="Period", ylab="# households")
abline(h=hist$n_hh[BURN_IN], lty=2, col="darkgreen")

plot(periods[RANGE], hist$n_births[RANGE], type="h", col="seagreen",
     main="Births per period", xlab="Period", ylab="# new HH")

plot(periods[RANGE], hist$unemp[RANGE]*hist$n_hh[RANGE], type="l", col="red", lwd=2,
     main="Unemployed (headcount)", xlab="Period", ylab="# persons")

# ── Firm age panels ───────────────────────────────────────────
plot(periods[RANGE], hist$mean_firm_age[RANGE], type="l", col="chocolate", lwd=2,
     main="Mean Firm Age (periods)", xlab="Period", ylab="periods")
lines(periods[RANGE], hist$med_firm_age[RANGE], col="chocolate4", lwd=2, lty=2)
legend("topleft", c("mean","median"), col=c("chocolate","chocolate4"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

plot(periods[RANGE], hist$share_young_firms[RANGE]*100, type="l", col="coral", lwd=2,
     main="Young Firms (age ≤ 12, %)", xlab="Period", ylab="%",
     ylim=c(0, 100))
abline(h=50, lty=3, col="grey60")

plot(periods[RANGE], hist$lg_g_spend[RANGE], type="l", col="navy", lwd=2,
     main="Gov. Procurement (G)", xlab="Period", ylab="currency units")
lines(periods[RANGE], hist$transfer_total[RANGE], col="steelblue", lwd=2, lty=2)
legend("topleft", c("G procurement","BF transfers"), col=c("navy","steelblue"),
       lwd=2, lty=c(1,2), cex=0.7, bty="n")

dev.off()
cat("Overview plot saved: abm_v5_stable2_overview_pop.png\n")


