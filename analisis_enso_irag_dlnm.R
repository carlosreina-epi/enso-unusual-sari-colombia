# =============================================================================
# ENSO e IRAG inusitada (evento 348, SIVIGILA) en Colombia, 2011-2019
# Modelo de rezagos distribuidos no lineales (DLNM) sobre series mensuales
# Autor: Carlos Alberto Reina Bolaños | Reanálisis noviembre 2026
#
# Datos de entrada (carpeta del proyecto):
#   irag348_semanal_estratos.csv  casos semanales del evento 348 por estrato
#                                 (fuente: datos.gov.co, recurso 4hyg-wa9d,
#                                 validados contra totales departamento-año)
#   raw/oni.csv                   Oceanic Niño Index mensual (NOAA CPC, ERSSTv5)
#   raw/pop_dpto_ano.csv          población DANE por departamento y año
#   raw/dept_strata.csv           región e IPM (terciles) por departamento
#
# La función cross_basis() reproduce la parametrización de dlnm::crossbasis()
# (Gasparrini 2011) con splines naturales de 'splines' (R base), para que el
# script corra sin dependencias externas. Con el paquete dlnm instalado,
# crossbasis(oni, lag = 12, argvar = list(fun = "ns", knots = kv),
# arglag = list(fun = "ns", knots = kl, intercept = TRUE)) da la misma base.
# =============================================================================

suppressPackageStartupMessages({ library(splines); library(ggplot2) })
# ---------------------------------------------------------------------------
# Carpeta del proyecto: si no se ejecuta desde ella, la busca dentro de
# Descargas (carpeta que contiene irag348_semanal_estratos.csv)
# ---------------------------------------------------------------------------
BASE <- "~/Downloads"
if (!file.exists("irag348_semanal_estratos.csv")) {
  hit <- list.files(BASE, pattern = "^irag348_semanal_estratos\\.csv$",
                    recursive = TRUE, full.names = TRUE)
  if (length(hit) == 0) stop("No encontré irag348_semanal_estratos.csv dentro de ", BASE)
  setwd(dirname(hit[1]))
}
cat("Carpeta de trabajo:", getwd(), "\n")
dir.create("resultados", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Semanas epidemiológicas -> meses
#    Semana 1 = semana (domingo a sábado) que contiene el 4 de enero.
#    Cada semana se asigna al mes de su miércoles (punto medio).
# ---------------------------------------------------------------------------
epi_week_mid <- function(year, week) {
  jan4 <- as.Date(sprintf("%d-01-04", year))
  start_w1 <- jan4 - as.integer(format(jan4, "%w"))   # domingo anterior o igual
  start_w1 + 7 * (week - 1) + 3                         # miércoles de la semana
}

wk <- read.csv("irag348_semanal_estratos.csv")
wk$mid <- epi_week_mid(wk$ano, wk$semana)
wk$ym  <- format(wk$mid, "%Y-%m")
estratos <- setdiff(names(wk), c("ano", "semana", "mid", "ym"))
mo <- aggregate(wk[, estratos], by = list(ym = wk$ym), FUN = sum)
mo <- mo[mo$ym >= "2011-01" & mo$ym <= "2019-12", ]
mo$fecha <- as.Date(paste0(mo$ym, "-15"))
mo$year  <- as.integer(substr(mo$ym, 1, 4))
mo$month <- as.integer(substr(mo$ym, 6, 7))
mo$time  <- seq_len(nrow(mo))
stopifnot(nrow(mo) == 108)

# ---------------------------------------------------------------------------
# 2. ONI mensual (temporada de 3 meses centrada en el mes)
# ---------------------------------------------------------------------------
oni_w <- read.csv("raw/oni.csv")
oni <- data.frame(year = rep(oni_w$year, each = 12), month = rep(1:12, nrow(oni_w)),
                  oni = as.vector(t(as.matrix(oni_w[, paste0("m", 1:12)]))))
oni$ym <- sprintf("%d-%02d", oni$year, oni$month)

# ---------------------------------------------------------------------------
# 3. Población (offset): suma departamental por estrato, interpolada a meses
# ---------------------------------------------------------------------------
pop <- read.csv("raw/pop_dpto_ano.csv")
st  <- read.csv("raw/dept_strata.csv")
reg_map <- c("Andi" = "Andina", "Caribe" = "Caribe")
st$reg <- ifelse(grepl("^Andi", st$region), "Andina",
          ifelse(grepl("^Caribe", st$region), "Caribe",
          ifelse(grepl("^Pac", st$region), "Pacifica", "Orinoquia_Amazonia")))
deps <- list(
  Nacional = st$Codigo_Departamento,
  Antioquia = 5,
  Nacional_sin_Antioquia = setdiff(st$Codigo_Departamento, 5),
  Andina = st$Codigo_Departamento[st$reg == "Andina"],
  Caribe = st$Codigo_Departamento[st$reg == "Caribe"],
  Pacifica = st$Codigo_Departamento[st$reg == "Pacifica"],
  Orinoquia_Amazonia = st$Codigo_Departamento[st$reg == "Orinoquia_Amazonia"],
  IPM_T1 = st$Codigo_Departamento[st$ipm_t == "T1"],
  IPM_T2 = st$Codigo_Departamento[st$ipm_t == "T2"],
  IPM_T3 = st$Codigo_Departamento[st$ipm_t == "T3"])
pop_month <- function(dp) {
  py <- aggregate(Total.General ~ AÑO, data = pop[pop$DP %in% dp, ], FUN = sum)
  # población anual ubicada a mitad de año e interpolada linealmente
  approx(x = py$AÑO + 0.5 - 1/24, y = py$Total.General,
         xout = mo$year + (mo$month - 1) / 12, rule = 2)$y
}
names(pop) <- make.names(names(pop))
for (s in estratos) mo[[paste0("pop_", s)]] <- pop_month(deps[[s]])

# ---------------------------------------------------------------------------
# 4. Cross-basis (equivalente a dlnm::crossbasis)
# ---------------------------------------------------------------------------
L_MAX   <- 12
oni_all <- oni[oni$ym >= "2009-01" & oni$ym <= "2019-12", ]
KV      <- quantile(oni_all$oni[oni_all$ym >= "2011-01"], c(1/3, 2/3))  # nudos exposición
BKV     <- range(oni_all$oni)
lag_knots <- function(L, nk) exp(log(L + 1) * seq_len(nk) / (nk + 1)) - 1  # = dlnm::logknots

VAR_FUN <- "lin"   # "lin": efecto lineal por 1 °C de ONI (principal); "ns": no lineal (sensibilidad)
basis_var <- function(x, knots = KV, bk = BKV, fun = VAR_FUN) {
  if (fun == "lin") return(matrix(x, ncol = 1))
  ns(x, knots = knots, Boundary.knots = bk)
}
basis_lag <- function(L, nk = 1) ns(0:L, knots = lag_knots(L, nk), intercept = TRUE,
                                    Boundary.knots = c(0, L))

cross_basis <- function(dates_ym, L = L_MAX, nk_lag = 1) {
  idx <- match(dates_ym, oni_all$ym)
  Q   <- sapply(0:L, function(l) oni_all$oni[idx - l])          # n x (L+1)
  C   <- basis_lag(L, nk_lag)                                    # (L+1) x dl
  dv  <- ncol(basis_var(0)); dl <- ncol(C)
  cb  <- matrix(0, length(idx), dv * dl)
  for (l in 0:L) {
    Bv <- basis_var(Q[, l + 1])                                  # n x dv
    for (i in 1:dv) for (j in 1:dl)
      cb[, (i - 1) * dl + j] <- cb[, (i - 1) * dl + j] + Bv[, i] * C[l + 1, j]
  }
  colnames(cb) <- paste0("cb_v", rep(1:dv, each = dl), "_l", rep(1:dl, dv))
  attr(cb, "C") <- C; attr(cb, "dv") <- dv; attr(cb, "dl") <- dl; attr(cb, "L") <- L
  cb
}

# Predicciones (equivalente a dlnm::crosspred) centradas en ONI = 0
cross_pred <- function(model, cb, at, cen = 0) {
  C <- attr(cb, "C"); dv <- attr(cb, "dv"); dl <- attr(cb, "dl"); L <- attr(cb, "L")
  VAR_FUN <<- attr(cb, "fun")
  nm <- colnames(cb); b <- coef(model)[nm]; V <- vcov(model)[nm, nm]
  Bx <- basis_var(at) - matrix(basis_var(cen), length(at), dv, byrow = TRUE)
  lagres <- list(); cum <- data.frame()
  for (k in seq_along(at)) {
    Xk <- t(sapply(0:L, function(l) as.vector(outer(Bx[k, ], C[l + 1, ]))))  # (L+1) x (dv*dl)
    est <- drop(Xk %*% b); se <- sqrt(pmax(diag(Xk %*% V %*% t(Xk)), 0))
    lagres[[k]] <- data.frame(oni = at[k], lag = 0:L, rr = exp(est),
                              lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
    xc <- colSums(Xk); ec <- sum(xc * b); sc <- sqrt(drop(t(xc) %*% V %*% xc))
    cum <- rbind(cum, data.frame(oni = at[k], rr = exp(ec), lo = exp(ec - 1.96 * sc),
                                 hi = exp(ec + 1.96 * sc)))
  }
  list(lag = do.call(rbind, lagres), cum = cum)
}

# ---------------------------------------------------------------------------
# 5. Modelo principal
#    casos ~ cross-basis ONI (lag 0-12) + tendencia ns(df) + 2 armónicos
#           + offset(log población), cuasi-Poisson
# ---------------------------------------------------------------------------
fit_dlnm <- function(y_name, L = L_MAX, df_trend = 3, nk_lag = 1, season = "fourier",
                     data = mo, ar = FALSE, fun = "lin") {
  VAR_FUN <<- fun
  cb <- cross_basis(data$ym, L, nk_lag)
  attr(cb, "fun") <- fun
  d  <- cbind(data, cb)
  d$y <- d[[y_name]]; d$lpop <- log(d[[paste0("pop_", y_name)]])
  harm <- "sin(2*pi*month/12) + cos(2*pi*month/12) + sin(4*pi*month/12) + cos(4*pi*month/12)"
  seas <- if (season == "fourier") harm else "factor(month)"
  f <- paste("y ~", paste(colnames(cb), collapse = " + "), "+ ns(time, df =", df_trend, ") +",
             seas, "+ offset(lpop)")
  if (ar) { d$lag_y <- log(c(NA, head(d$y, -1)) + 1); f <- paste(f, "+ lag_y") }
  m <- glm(as.formula(f), family = quasipoisson, data = d)
  list(model = m, cb = cb, data = d)
}

QAIC <- function(m) {                   # Peng et al. 2006
  phi <- summary(m)$dispersion
  ll  <- sum(dpois(m$y, m$fitted.values, log = TRUE))
  -2 * ll + 2 * phi * length(coef(m))
}

AT  <- round(seq(-1.5, 2.5, by = 0.1), 1)
KEY <- c(-1.5, -1.0, 1.0, 1.5, 2.0, 2.5)

main <- fit_dlnm("Nacional")
pm   <- cross_pred(main$model, main$cb, AT)
cat("\nDispersión (phi):", round(summary(main$model)$dispersion, 2),
    " | QAIC:", round(QAIC(main$model), 1), "\n")

# ---------------------------------------------------------------------------
# 6. Descriptivos
# ---------------------------------------------------------------------------
desc <- merge(mo[, c("ym", "fecha", "Nacional", "pop_Nacional")], oni[, c("ym", "oni")], by = "ym")
desc$tasa <- desc$Nacional / desc$pop_Nacional * 1e6
write.csv(desc, "resultados/serie_mensual_nacional.csv", row.names = FALSE)

tab1 <- data.frame(
  estrato = estratos,
  casos_2011_2019 = sapply(estratos, function(s) sum(mo[[s]])),
  media_mensual = sapply(estratos, function(s) round(mean(mo[[s]]), 1)),
  meses_con_cero = sapply(estratos, function(s) sum(mo[[s]] == 0)),
  tasa_anual_por_millon = sapply(estratos, function(s)
    round(sum(mo[[s]]) / mean(mo[[paste0("pop_", s)]]) / 9 * 1e6, 1)))
write.csv(tab1, "resultados/tabla1_descriptivos.csv", row.names = FALSE)

sc <- max(desc$tasa) / 3
g1 <- ggplot(desc, aes(fecha)) +
  geom_rect(data = desc, aes(xmin = fecha - 15, xmax = fecha + 15, ymin = -Inf, ymax = Inf,
            fill = ifelse(oni >= 0.5, "El Niño", ifelse(oni <= -0.5, "La Niña", "Neutral"))),
            alpha = 0.18, inherit.aes = FALSE) +
  geom_line(aes(y = tasa), colour = "grey15", linewidth = 0.6) +
  geom_line(aes(y = oni * sc + sc), colour = "#B2182B", linewidth = 0.5, linetype = "22") +
  scale_y_continuous("Unusual SARI, cases per million inhabitants per month",
                     sec.axis = sec_axis(~ (. - sc) / sc, name = "ONI (°C)")) +
  scale_fill_manual(NULL, values = c("El Niño" = "#D6604D", "La Niña" = "#4393C3", "Neutral" = "white")) +
  labs(x = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave("resultados/fig1_serie_irag_oni.png", g1, width = 9, height = 4.5, dpi = 300)

# ---------------------------------------------------------------------------
# 7. Figuras del modelo principal
# ---------------------------------------------------------------------------
g2 <- ggplot(pm$cum, aes(oni, rr)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, linetype = 2) +
  scale_y_log10() + labs(x = "ONI (°C)", y = "Cumulative RR, lags 0-12 months (ref. ONI = 0)") +
  theme_minimal(base_size = 11)
ggsave("resultados/fig2_exposicion_respuesta_acumulada.png", g2, width = 6, height = 4.5, dpi = 300)

lagsel <- pm$lag[pm$lag$oni %in% c(-1.5, 1.5, 2.5), ]
lagsel$panel <- factor(paste("ONI =", lagsel$oni), levels = paste("ONI =", c(-1.5, 1.5, 2.5)))
g3 <- ggplot(lagsel, aes(lag, rr)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") + geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~panel) + scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 12, 3)) +
  labs(x = "Lag (months)", y = "Lag-specific RR (ref. ONI = 0)") + theme_minimal(base_size = 11)
ggsave("resultados/fig3_rezago_respuesta.png", g3, width = 9, height = 3.8, dpi = 300)

# Superficie 3D (exposición x rezago) como mapa de calor
g4 <- ggplot(pm$lag, aes(lag, oni, fill = log(rr))) + geom_tile() +
  scale_fill_gradient2("log RR", low = "#2166AC", mid = "white", high = "#B2182B") +
  scale_x_continuous(breaks = seq(0, 12, 3)) + labs(x = "Lag (months)", y = "ONI (°C)") +
  theme_minimal(base_size = 11)
ggsave("resultados/fig4_superficie_exposicion_rezago.png", g4, width = 6, height = 4.5, dpi = 300)

# ---------------------------------------------------------------------------
# 8. Estratos (región, IPM) y sensibilidad
# ---------------------------------------------------------------------------
fmt <- function(r) sprintf("%.2f (%.2f-%.2f)", r$rr, r$lo, r$hi)
summ <- function(fit, label) {
  p <- cross_pred(fit$model, fit$cb, KEY)$cum
  data.frame(analisis = label, oni = p$oni, rr = round(p$rr, 3), lo = round(p$lo, 3),
             hi = round(p$hi, 3), rr_ic = fmt(p), phi = round(summary(fit$model)$dispersion, 2),
             qaic = round(QAIC(fit$model), 1))
}
res <- summ(main, "Principal: nacional, ONI lineal, lag 0-12 (ns 3 gl), tendencia df=3, 2 armónicos")
for (s in setdiff(estratos, "Nacional")) res <- rbind(res, summ(fit_dlnm(s), paste("Estrato:", s)))
sens <- list(
  "Sens: lag 0-6"            = fit_dlnm("Nacional", L = 6),
  "Sens: exposición no lineal (ns, 2 nudos)" = fit_dlnm("Nacional", fun = "ns"),
  "Sens: lag 0-18"           = fit_dlnm("Nacional", L = 18),
  "Sens: tendencia df=2"     = fit_dlnm("Nacional", df_trend = 2),
  "Sens: tendencia df=5"     = fit_dlnm("Nacional", df_trend = 5),
  "Sens: tendencia df=9 (1/año)" = fit_dlnm("Nacional", df_trend = 9),
  "Sens: estacionalidad con mes (factor)" = fit_dlnm("Nacional", season = "month"),
  "Sens: término autorregresivo log(y[t-1]+1)" = fit_dlnm("Nacional", ar = TRUE))
for (nm in names(sens)) res <- rbind(res, summ(sens[[nm]], nm))
write.csv(res, "resultados/tabla2_rr_acumulados.csv", row.names = FALSE)

# Diagnóstico: autocorrelación de residuos del modelo principal
rd <- residuals(main$model, type = "deviance")
acf_v <- acf(rd, plot = FALSE, lag.max = 12)$acf[-1]
write.csv(data.frame(lag = 1:12, acf = round(acf_v, 3)), "resultados/acf_residuos_principal.csv",
          row.names = FALSE)

# Forest plot de estratos y sensibilidad para ONI = +2.0 y -1.5
fp <- res[res$oni %in% c(-1.5, 2.0), ]
fp$analisis <- factor(fp$analisis, levels = rev(unique(res$analisis)))
g5 <- ggplot(fp, aes(rr, analisis)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) + geom_point() +
  facet_wrap(~ paste("ONI =", oni)) + scale_x_log10() +
  coord_cartesian(xlim = c(0.1, 10)) +   # el modelo no lineal tiene IC fuera de escala (ver tabla 2)
  labs(x = "Cumulative RR, lags 0-12 months (ref. ONI = 0)", y = NULL) + theme_minimal(base_size = 10)
ggsave("resultados/fig5_forest_estratos_sensibilidad.png", g5, width = 10, height = 6.5, dpi = 300)

saveRDS(list(main = main, pm = pm, res = res), "resultados/modelos.rds")
print(res[res$oni %in% c(-1.5, 1.5, 2.5), c("analisis", "oni", "rr_ic", "phi")], row.names = FALSE)
cat("\nACF residuos (lags 1-6):", round(acf_v[1:6], 2), "\n")
