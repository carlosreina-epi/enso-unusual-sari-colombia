# =============================================================================
# Verificación del modelo principal con el paquete dlnm oficial
# Ejecutar: Rscript verificar_dlnm.R  (encuentra sola la carpeta dentro de /Users/carlos.reina/Downloads)
# Requiere: install.packages(c("ggplot2", "dlnm"))  (una sola vez, en la consola)
# =============================================================================

# Si no se ejecuta desde la carpeta del proyecto, la busca en Descargas
if (!file.exists("analisis_enso_irag_dlnm.R")) {
  BASE <- "~/Downloads"
  hit <- list.files(BASE, pattern = "^analisis_enso_irag_dlnm\\.R$", recursive = TRUE, full.names = TRUE)
  if (length(hit) == 0) stop("Ejecuta este script desde la carpeta del proyecto")
  setwd(dirname(hit[1]))
}
source("analisis_enso_irag_dlnm.R", encoding = "UTF-8")   # corre el análisis completo
library(dlnm)

# Exposiciones rezagadas (ONI de t-0 a t-12) para cada mes 2011-2019
idx <- match(mo$ym, oni_all$ym)
Q   <- sapply(0:12, function(l) oni_all$oni[idx - l])

# Cross-basis oficial: ONI lineal; rezago con spline natural y 1 nudo logarítmico
cb <- crossbasis(Q, lag = 12,
                 argvar = list(fun = "lin"),
                 arglag = list(fun = "ns", knots = logknots(12, nk = 1), intercept = TRUE))

m_dlnm <- glm(Nacional ~ cb + ns(time, df = 3) +
                sin(2*pi*month/12) + cos(2*pi*month/12) +
                sin(4*pi*month/12) + cos(4*pi*month/12) +
                offset(log(pop_Nacional)),
              family = quasipoisson, data = mo)

pred <- crosspred(cb, m_dlnm, at = 1, cen = 0)
cat("\n=== Verificación con dlnm ===\n")
print(round(c(RR = unname(pred$allRRfit), IC_inf = unname(pred$allRRlow),
              IC_sup = unname(pred$allRRhigh)), 2))
cat("Esperado (implementación propia): RR 0.95, IC 0.78 a 1.15\n")

# RR por rezago con dlnm, para comparar con la Figura 3
print(round(cbind(lag = 0:12, RR = pred$matRRfit[1, ], IC_inf = pred$matRRlow[1, ],
                  IC_sup = pred$matRRhigh[1, ]), 3))

png("resultados/verificacion_dlnm_rezago.png", width = 1800, height = 1200, res = 200)
plot(pred, "slices", var = 1, xlab = "Lag (months)", ylab = "RR per 1 °C ONI")
dev.off()
cat("\nGráfico guardado en resultados/verificacion_dlnm_rezago.png\n")
