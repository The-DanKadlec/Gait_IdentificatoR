# ============================================================
# AUTOMATED GPS-ACCELEROMETER GAIT-EVENT IDENTIFICATION
# ============================================================
# Detects gait impact landmarks from a single trunk-mounted
# accelerometer/GPS (resultant signal) and exports step and/or
# stride interval time-series.
#
# Input:  a folder of CSV files, each with a time column and a
#         resultant-acceleration column.
# Output: a new sub-folder containing one CSV per input file,
#         plus a run-summary CSV. Missing intervals are kept as
#         explicit NA (never interpolated) and reported.
# ============================================================

# %% ---- packages -------------------------------------------
req <- c("dplyr", "tibble", "readr", "signal")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))
library(dplyr); library(tibble); library(readr)

# %% ---- USER SETTINGS ( edit these ) -----------------------

# Folders
input_dir  <- "path/to/csv/folder"                       # your input CSVs
output_dir <- file.path(input_dir, "intervals_output")   # created automatically

# Column names in your CSVs
time_col <- "time_s"        # time column, in seconds
acc_col  <- "resultant"     # resultant-acceleration column

# What to export: "step", "stride", or "both"
export_type <- "both"

# Sampling frequency (Hz) - MUST match your recordings
fs <- 100

# Running-bout gate: keeps the longest steady-state running segment,
# excluding standing, walking, and turns.
gate_win_s    <- 1.0        # activity-envelope smoothing window (s)
gate_frac     <- 0.35       # threshold as fraction of median envelope
min_bout_s    <- 10         # minimum acceptable bout length (s)
edge_buffer_s <- 5          # seconds trimmed from each end (accel/decel)

# Detection band (Hz): isolates the impact transient. Lower edge sits
# above the step-frequency harmonics; upper edge below sensor noise.
band_lo <- 15
band_hi <- 40
filt_order <- 4             # Butterworth order, applied zero-phase (filtfilt)

# Physiological interval limits (s): intervals outside this range are
# treated as errors and set to NA, never interpolated. Widen only for
# very slow or fast gait. Keep active.
dt_lo <- 0.20
dt_hi <- 0.60

# Viterbi sequence selection: balances peak prominence against
# deviation from the locally expected step period. Higher lambda =
# stronger adherence to the expected period.
sd_prior   <- 0.05          # expected step-period tolerance (s)
lambda     <- 3             # prior strength
start_cost <- 60            # sequence entry cost

# Sub-sample refinement (Lanczos-windowed sinc): refines each landmark
# below the sample grid. Larger values = finer but slower.
ref_half_range <- 3         # search range about each landmark (samples)
ref_step       <- 0.02      # refinement grid step (samples)
ref_taps       <- 20        # sinc kernel half-width (samples)

# Gap recovery: reinstates genuine landmarks the Viterbi sequence
# skipped (intervals near 2x the expected period). Reinstates only
# REAL detected peaks; never fabricates values. Set to 0 to disable.
gap_tol    <- 0.20
max_insert <- 3

# %% ---- helper functions ( no editing needed ) -------------

bpf <- function(x, lo, hi, fsx, order = filt_order) {
  bf <- signal::butter(order, c(lo, hi) / (fsx / 2), type = "pass")
  as.numeric(signal::filtfilt(bf, x))
}

roll_mean <- function(x, n) {
  k <- floor(n / 2); cs <- cumsum(c(0, x)); i <- seq_along(x)
  a <- pmax(i - k, 1); b <- pmin(i + k, length(x)); (cs[b + 1] - cs[a]) / (b - a + 1)
}

local_peaks <- function(x, half_win) {
  cand <- which(diff(sign(diff(x))) == -2) + 1
  cand <- cand[cand > half_win & cand < (length(x) - half_win)]
  prom <- vapply(cand, function(k)
    x[k] - max(min(x[(k - half_win):k]), min(x[k:(k + half_win)])), numeric(1))
  list(idx = cand, prom = prom)
}

local_period <- function(x, fsx, win_s = 4, hop_s = 1, lo = dt_lo, hi = dt_hi) {
  xl <- bpf(x, 1, 10, fsx); w <- round(win_s * fsx); h <- round(hop_s * fsx)
  st <- seq(1, max(1, length(xl) - w), by = h); ct <- pr <- numeric(length(st))
  for (i in seq_along(st)) {
    s <- xl[st[i]:(st[i] + w - 1)]; s <- s - mean(s)
    a <- as.numeric(stats::acf(s, lag.max = round(hi * fsx), plot = FALSE)$acf)
    rng <- (round(lo * fsx) + 1):(round(hi * fsx) + 1)
    pr[i] <- (rng[which.max(a[rng])] - 1) / fsx; ct[i] <- (st[i] + w / 2) / fsx
  }
  approxfun(ct, pr, rule = 2)
}

viterbi_chain <- function(ct, score, pf, dt_lo, dt_hi, sd_prior, lambda, start_cost) {
  n <- length(ct); best <- score - start_cost; prev <- rep(-1L, n)
  for (j in seq_len(n)) {
    lo <- findInterval(ct[j] - dt_hi, ct) + 1; hi <- findInterval(ct[j] - dt_lo, ct)
    if (hi < lo) next; Tj <- pf(ct[j])
    for (i in lo:hi) {
      pen <- lambda * ((ct[j] - ct[i] - Tj) / sd_prior)^2; v <- best[i] + score[j] - pen
      if (v > best[j]) { best[j] <- v; prev[j] <- i }
    }
  }
  k <- which.max(best); path <- integer(0)
  while (k > 0) { path <- c(k, path); k <- prev[k] }; path
}

sinc <- function(d) ifelse(abs(d) < 1e-12, 1, sin(pi * d) / (pi * d))
sinc_eval <- function(x, tq, taps) {
  lo <- max(1, floor(min(tq)) - taps); hi <- min(length(x), ceiling(max(tq)) + taps)
  kk <- lo:hi; D <- outer(tq, kk, "-")
  W <- sinc(D) * ifelse(abs(D) < taps, sinc(D / taps), 0); as.numeric(W %*% x[kk])
}
refine_sinc <- function(x, idx0, half_range, step, taps) {
  offs <- seq(-half_range, half_range, by = step)
  vapply(idx0, function(i0) { tq <- i0 + offs
  if (min(tq) - taps < 1 || max(tq) + taps > length(x)) return(NA_real_)
  tq[which.max(sinc_eval(x, tq, taps))] }, numeric(1))
}

# %% ---- core steps -----------------------------------------

steady_segment <- function(tt, acc) {
  env <- roll_mean(abs(bpf(acc, 1, 20, fs)), round(gate_win_s * fs))
  thr <- gate_frac * median(env); m <- env > thr
  rr <- rle(m); ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1
  ok <- which(rr$values & rr$lengths >= min_bout_s * fs)
  if (!length(ok)) stop("no running bout found")
  j <- ok[which.max(rr$lengths[ok])]
  t0 <- tt[starts[j]] + edge_buffer_s; t1 <- tt[ends[j]] - edge_buffer_s
  if (t1 - t0 < min_bout_s) stop("bout too short after edge exclusion")
  list(seg = which(tt >= t0 & tt <= t1), t0 = t0, t1 = t1)
}

detect_events <- function(a, t0, pf) {
  xf <- bpf(a, band_lo, band_hi, fs)
  pk <- local_peaks(xf, half_win = round(0.08 * fs))
  sc <- pk$prom / median(pk$prom); ct <- pk$idx / fs
  p  <- viterbi_chain(ct, sc, pf, dt_lo, dt_hi, sd_prior, lambda, start_cost)
  ic <- refine_sinc(xf, pk$idx[p], ref_half_range, ref_step, ref_taps) / fs + t0
  list(ic = ic[is.finite(ic)], xf = xf, cand_idx = pk$idx)
}

recover_gaps <- function(ic_t, cand_idx, xf, pf, t0) {
  if (max_insert < 1) return(ic_t)
  cand_t <- cand_idx / fs + t0; inserted <- numeric(0); i <- 1
  while (i < length(ic_t)) {
    gap <- ic_t[i + 1] - ic_t[i]; Tloc <- pf(ic_t[i]); k <- round(gap / Tloc)
    if (k >= 2 && k <= (max_insert + 1) && abs(gap - k * Tloc) < gap_tol * k * Tloc) {
      for (mm in seq_len(k - 1)) {
        te <- ic_t[i] + mm * gap / k; win <- which(abs(cand_t - te) < 0.5 * Tloc)
        win <- win[!(cand_t[win] %in% c(ic_t, inserted))]
        if (length(win)) inserted <- c(inserted, cand_t[win[which.max(xf[cand_idx[win]])]])
      }
    }
    i <- i + 1
  }
  if (length(inserted)) {
    ii <- refine_sinc(xf, round((inserted - t0) * fs), ref_half_range, ref_step, ref_taps) / fs + t0
    sort(c(ic_t, ii[is.finite(ii)]))
  } else ic_t
}

build_intervals <- function(ic) {
  n <- length(ic)
  step_time   <- c(NA, diff(ic))
  stride_time <- c(NA, NA, ic[3:n] - ic[1:(n - 2)])
  bad_step <- !is.na(step_time) & (step_time < dt_lo | step_time > dt_hi)
  step_time[bad_step] <- NA
  bad_stride <- c(FALSE, bad_step[-n]) | bad_step
  stride_time[bad_stride] <- NA
  tibble(landmark = seq_len(n), landmark_time_s = round(ic, 4),
         step_time_s = round(step_time, 4), stride_time_s = round(stride_time, 4))
}

# %% ---- per-file processing --------------------------------

process_file <- function(path) {
  d <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c(time_col, acc_col) %in% names(d)))
    stop(sprintf("missing '%s' or '%s'", time_col, acc_col))
  tt <- d[[time_col]]; acc <- d[[acc_col]]
  sg  <- steady_segment(tt, acc)
  seg <- acc[sg$seg]; t0 <- sg$t0
  pf  <- local_period(seg, fs)
  det <- detect_events(seg, t0, pf)
  n_pre <- length(det$ic)
  ic  <- recover_gaps(det$ic, det$cand_idx, det$xf, pf, t0)
  list(intervals = build_intervals(ic),
       n_landmarks = length(ic), n_recovered = length(ic) - n_pre)
}

# %% ---- batch run ------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
files <- files[normalizePath(dirname(files)) != normalizePath(output_dir)]
if (!length(files)) stop("no CSV files found in input_dir")

summary_rows <- list()

for (f in files) {
  base <- tools::file_path_sans_ext(basename(f))
  res <- tryCatch(process_file(f), error = function(e) {
    message(sprintf("SKIP %s: %s", base, conditionMessage(e))); NULL
  })
  if (is.null(res)) {
    summary_rows[[base]] <- tibble(file = base, status = "failed",
                                   landmarks = NA, recovered = NA, step_gaps = NA, stride_gaps = NA)
    next
  }
  
  iv <- res$intervals
  step_gaps   <- sum(is.na(iv$step_time_s[-1]))
  stride_gaps <- sum(is.na(iv$stride_time_s[-(1:2)]))
  
  if (step_gaps > 0) {
    pos <- iv$landmark[which(is.na(iv$step_time_s))][-1]
    warning(sprintf("%s: %d step gap(s) at landmark(s) %s - inspect raw signal.",
                    base, step_gaps, paste(head(pos, 20), collapse = ", ")), call. = FALSE)
  }
  message(sprintf("%s: %d landmarks | %d recovered | %d step gaps | %d stride gaps",
                  base, res$n_landmarks, res$n_recovered, step_gaps, stride_gaps))
  
  # export with NAs preserved (position-preserving; no gap dropped)
  if (export_type %in% c("step", "both"))
    readr::write_csv(iv[, c("landmark", "landmark_time_s", "step_time_s")],
                     file.path(output_dir, paste0(base, "_step.csv")))
  if (export_type %in% c("stride", "both"))
    readr::write_csv(iv[, c("landmark", "landmark_time_s", "stride_time_s")],
                     file.path(output_dir, paste0(base, "_stride.csv")))
  
  summary_rows[[base]] <- tibble(file = base, status = "ok",
                                 landmarks = res$n_landmarks, recovered = res$n_recovered,
                                 step_gaps = step_gaps, stride_gaps = stride_gaps)
}

readr::write_csv(bind_rows(summary_rows), file.path(output_dir, "_run_summary.csv"))
message("done -> ", output_dir)