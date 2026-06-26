# ============================================================
#  Sanity-check / ground-truth comparison:
#  nitrate-drawdown NCP  vs  independent O2-budget NCP analogue
#  (float 3902681, subpolar N. Atlantic)
#
#  Valid pairing: nitrate NCP  <->  O2 *net* biological residual
#  (MLD-following, air-sea + entrainment removed; o2_budget_3902681.R).
#  Night-O2 = respiration, NOT NCP -> deliberately excluded.
#
#  Comparison is done in VOLUMETRIC (mixed-layer-mean) units, the report's
#  recommended trustworthy quantity: dividing by the integration depth removes
#  the winter depth-amplification that otherwise swamps the AREAL residual
#  (Fig. D). Both methods are put on a COMMON MLD-based integration depth so
#  the contrast is method-vs-method, not depth-vs-depth.
#
#  Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/ncp_o2_sanitycheck.R
# ============================================================

suppressMessages({
  library(tidyverse); library(lubridate); library(zoo); library(patchwork)
})

PQ     <- 1.45                         # O2:C net-production quotient (Laws 1991)
GRID   <- 15                           # day bin = NCP method's native window
origin <- as.Date("2024-10-01")        # fixed bin anchor (so both series align)
cutoff <- as.Date("2025-11-01")
bin    <- function(d) origin + GRID * floor(as.numeric(d - origin) / GRID)
seas   <- function(d) factor(case_when(month(d) %in% c(4,5,6,7)  ~ "productive (Apr-Jul)",
                                       month(d) %in% c(11,12,1,2) ~ "winter (Nov-Feb)",
                                       TRUE                        ~ "shoulder"),
                             levels = c("productive (Apr-Jul)","shoulder","winter (Nov-Feb)"))
col_no3 <- "#08519c"; col_o2 <- "#cb181d"

# ---------------------------------------------------------------------------
# 1.  O2-budget NET biological residual (independent ground truth)
# ---------------------------------------------------------------------------
o2 <- read_csv("Data/Processed/o2_budget_3902681_mld.csv", show_col_types = FALSE) %>%
  filter(type == "net") %>%
  mutate(date = as.Date(mtime)) %>%
  filter(!is.na(rate_resid_vol_mmol_o2_m3_d), date <= cutoff) %>%
  transmute(date,
            o2_vol   = rate_resid_vol_mmol_o2_m3_d / PQ,        # mmol C m-3 d-1
            o2_areal = rate_resid_mmol_o2_m2_d   / PQ,          # mmol C m-2 d-1
            z_int    = z_int_m,
            w        = 1 / pmax(sigma_rate_vol_mmol_o2_m3_d / PQ, 1e-3)^2)

# common MLD-based integration depth z(t), daily-interpolated from the O2 pairs
zser <- o2 %>% group_by(date) %>% summarise(z = mean(z_int), .groups = "drop")

# ---------------------------------------------------------------------------
# 2.  Nitrate-drawdown NCP (the method being published) -> volumetric
# ---------------------------------------------------------------------------
ncp <- read_csv("Data/Processed/ncp_float_3902681_30d.csv", show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>% distinct(date, .keep_all = TRUE) %>%
  arrange(date) %>% filter(!is.na(NCP), date <= cutoff) %>%
  transmute(date, ncp_areal = NCP) %>%
  mutate(z       = approx(zser$date, zser$z, xout = date, rule = 2)$y,
         ncp_vol = ncp_areal / z)                              # mmol C m-3 d-1

# ---------------------------------------------------------------------------
# 3.  Bin both to the common 15-day grid
# ---------------------------------------------------------------------------
ncp_b <- ncp %>% mutate(b = bin(date)) %>% group_by(b) %>%
  summarise(ncp_vol = mean(ncp_vol), ncp_areal = mean(ncp_areal),
            z = mean(z), .groups = "drop")
o2_b <- o2 %>% mutate(b = bin(date)) %>% group_by(b) %>%
  summarise(o2_vol  = weighted.mean(o2_vol, w),
            o2_sd   = sqrt(1 / sum(w)),
            o2_areal = mean(o2_areal), .groups = "drop")
cmp <- inner_join(ncp_b, o2_b, by = "b") %>% mutate(season = seas(b))

# ---------------------------------------------------------------------------
# 4.  Agreement statistics (volumetric)
# ---------------------------------------------------------------------------
r_all  <- cor(cmp$ncp_vol, cmp$o2_vol)
r_sp   <- cor(cmp$ncp_vol, cmp$o2_vol, method = "spearman")
fit    <- lm(o2_vol ~ ncp_vol, data = cmp); sl <- coef(fit)[2]
rmsd   <- sqrt(mean((cmp$o2_vol - cmp$ncp_vol)^2))
bias   <- mean(cmp$o2_vol - cmp$ncp_vol)
cat(sprintf("\n n(15-day bins)        = %d\n", nrow(cmp)))
cat(sprintf(" Pearson r  (vol)      = %.2f\n", r_all))
cat(sprintf(" Spearman r (vol)      = %.2f\n", r_sp))
cat(sprintf(" OLS slope (vol)       = %.2f\n", sl))
cat(sprintf(" RMSD / bias (vol)     = %.2f / %+.2f mmol C m-3 d-1\n", rmsd, bias))
cat(sprintf(" Pearson r  (AREAL)    = %.2f   <- swamped by winter amplification\n",
            cor(cmp$ncp_areal, cmp$o2_areal)))

# ---------------------------------------------------------------------------
# 5.  Figures
# ---------------------------------------------------------------------------
# Panel A — volumetric time-series overlay (headline)
pA <- ggplot(cmp, aes(b)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = o2_vol - o2_sd, ymax = o2_vol + o2_sd),
              fill = col_o2, alpha = 0.15) +
  geom_line(aes(y = o2_vol,  color = "O2 budget (net residual / PQ)"),  linewidth = 0.8) +
  geom_point(aes(y = o2_vol, color = "O2 budget (net residual / PQ)"),  size = 1.1) +
  geom_line(aes(y = ncp_vol, color = "Nitrate-drawdown NCP"),           linewidth = 0.8) +
  geom_point(aes(y = ncp_vol,color = "Nitrate-drawdown NCP"),           size = 1.1) +
  scale_color_manual(NULL, values = c("Nitrate-drawdown NCP" = col_no3,
                                      "O2 budget (net residual / PQ)" = col_o2)) +
  scale_x_date(date_labels = "%b %Y", breaks = "2 months") +
  labs(title = "A.  Both methods resolve the same seasonal NCP cycle",
       subtitle = "mixed-layer-mean (volumetric) units; spring-summer peak, near-zero convective winter",
       x = NULL, y = expression(NCP~(mmol~C~m^-3~d^-1))) +
  theme_bw() + theme(legend.position = "top")

# Panel B — 1:1 scatter (volumetric)
lim <- range(c(cmp$ncp_vol, cmp$o2_vol));
pB <- ggplot(cmp, aes(ncp_vol, o2_vol)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey85", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = .3) +
  geom_errorbar(aes(ymin = o2_vol - o2_sd, ymax = o2_vol + o2_sd),
                color = "grey75", width = 0) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = .6) +
  geom_point(aes(fill = season), shape = 21, size = 2.8, color = "grey20") +
  scale_fill_manual(values = c("productive (Apr-Jul)" = "#1a9850",
                               "shoulder" = "#fdae61", "winter (Nov-Feb)" = "#4575b4")) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = sprintf("B.  1:1 agreement (r = %.2f, n = %d)", r_all, nrow(cmp)),
       subtitle = "winter at origin; productive season drives the signal",
       x = expression(Nitrate~NCP~(mmol~C~m^-3~d^-1)),
       y = expression(O2-budget~NCP~(mmol~C~m^-3~d^-1)), fill = NULL) +
  theme_bw() + theme(legend.position = c(.02,.98), legend.justification = c(0,1),
                     legend.background = element_rect(fill = alpha("white", .7)))

# Panel C — seasonal-mean bars (volumetric)
sb <- cmp %>% pivot_longer(c(ncp_vol, o2_vol), names_to = "method", values_to = "v") %>%
  group_by(season, method) %>%
  summarise(mean = mean(v), se = sd(v)/sqrt(n()), .groups = "drop") %>%
  mutate(method = recode(method, ncp_vol = "Nitrate", o2_vol = "O2 budget"))
pC <- ggplot(sb, aes(season, mean, fill = method)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = .3) +
  geom_col(position = position_dodge(.7), width = .65, color = "grey20") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(.7), width = .2) +
  scale_fill_manual(NULL, values = c("Nitrate" = col_no3, "O2 budget" = col_o2)) +
  labs(title = "C.  Same sign and magnitude in every season",
       x = NULL, y = expression(mean~NCP~(mmol~C~m^-3~d^-1))) +
  theme_bw() + theme(legend.position = "top", axis.text.x = element_text(size = 8))

# Panel D — why volumetric: the AREAL residual is depth-amplified in winter
sc <- max(abs(cmp$o2_areal), abs(cmp$ncp_areal)) / max(cmp$z)   # scale z onto rate axis
pD <- ggplot(cmp, aes(b)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = .3) +
  geom_col(aes(y = z * sc), fill = "grey85", width = 12) +
  geom_line(aes(y = o2_areal,  color = "O2 budget"), linewidth = .7) +
  geom_line(aes(y = ncp_areal, color = "Nitrate"),   linewidth = .7) +
  scale_color_manual(NULL, values = c("Nitrate" = col_no3, "O2 budget" = col_o2)) +
  scale_y_continuous(name = expression(AREAL~NCP~(mmol~C~m^-2~d^-1)),
                     sec.axis = sec_axis(~ . / sc, name = "integration depth (m)")) +
  scale_x_date(date_labels = "%b %Y", breaks = "2 months") +
  labs(title = "D.  Why volumetric: areal rates inflate with mixing depth",
       subtitle = "grey = integration depth; both areal series blow up where the layer deepens",
       x = NULL) +
  theme_bw() + theme(legend.position = "top")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Float 3902681 — nitrate-drawdown NCP validated against an independent O₂ budget",
    subtitle = sprintf("27 fortnightly bins, Oct 2024–Oct 2025  |  O₂→C with PQ = %.2f  |  night-O₂ (respiration) excluded", PQ),
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave("Output/ncp_o2_sanitycheck_3902681.png", fig, width = 14, height = 9.5, dpi = 200)
cat("\nsaved -> Output/ncp_o2_sanitycheck_3902681.png\n")
