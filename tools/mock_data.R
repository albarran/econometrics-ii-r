# Stand-in data for TESTING the notebooks while Henry's files are not published.
# Same file names and variable names as the slides; the numbers are NOT Henry's.
# Usage: Rscript tools/mock_data.R <output dir>
# Then run the notebooks with the environment variable ECMTII_DATA=<output dir>.
out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- "mock_data"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

# ---- Unit 3: transparency_firms.csv and transparency_oracle.csv ----
n <- 5000
arm <- sample(rep(0:2, c(2500, 1250, 1250)))
culture <- rnorm(n)
sector <- sample(paste0("S", 1:6), n, replace = TRUE)
province <- sample(c("Alicante", "Castellon", "Valencia"), n, replace = TRUE, prob = c(.35, .15, .5))
size <- pmax(10, round(exp(rnorm(n, log(44), 0.75))))
n_women <- pmin(size - 1, pmax(1, round(size * plogis(rnorm(n, 0, 0.5)))))
n_men <- size - n_women
p <- 0.25 + 0.05 * (arm == 1) + 0.24 * (arm == 2) + 0.17 * culture
transp <- as.integer(runif(n) < pmin(pmax(p, 0.01), 0.99))
parent_abroad <- as.integer(runif(n) < plogis(-2.3 + 0.5 * culture))
cell <- paste(sector, province)
peer_share <- pmin(pmax(0.32 + 0.06 * culture + ave(rnorm(n), cell) + rnorm(n, 0, 0.12), 0), 1)
fair_score <- pmin(pmax(5.9 + 1.1 * culture + rnorm(n, 0, 1.03), 0), 10)
women_mgmt <- pmin(pmax(35 + 8 * culture + rnorm(n, 0, 12.5), 0), 100)
gap2021 <- 0.163 + 0.03 * rnorm(n) + 0.3 * rnorm(n) / sqrt(size)
wage_men_2021 <- round(2150 * exp(rnorm(n, 0, 0.15)), 2)
wage_women_2021 <- round(wage_men_2021 / exp(gap2021), 2)
u <- 0.42 * rnorm(n) / sqrt(size)
g_men <- 0.12 - 0.029 * transp - 0.015 * culture + 0.03 * rnorm(n) + u / 2
g_women <- 0.10 + 0.003 * transp + 0.015 * culture + 0.03 * rnorm(n) - u / 2 +
  (g_men - (0.12 - 0.029 * transp - 0.015 * culture + u / 2))
wage_men_2025 <- round(wage_men_2021 * exp(g_men), 2)
wage_women_2025 <- round(wage_women_2021 * exp(g_women), 2)
firms <- data.frame(id = 1:n, wage_men_2021, wage_women_2021, wage_men_2025, wage_women_2025,
                    transp, arm, size, n_women, n_men, sector, province, parent_abroad,
                    peer_share = round(peer_share, 4), fair_score = round(fair_score, 2),
                    women_mgmt = round(women_mgmt, 1))
write.csv(firms, file.path(out, "transparency_firms.csv"), row.names = FALSE)
write.csv(data.frame(id = 1:n, culture), file.path(out, "transparency_oracle.csv"), row.names = FALSE)

# ---- Unit 4: workers.csv and firms.csv ----
n <- 3000
ability <- rnorm(n)
fam <- rnorm(n)
meduc <- round(pmin(pmax(10 + 2.5 * fam + rnorm(n, 0, 2), 4), 20))
feduc <- round(pmin(pmax(10 + 2.5 * fam + rnorm(n, 0, 2.5), 4), 20))
dist <- round(rexp(n, 1 / 25), 1)
sibs <- pmax(0, round(2.2 - 0.6 * ability + rnorm(n)))
educ <- round(pmin(pmax(9.5 + 0.18 * meduc + 0.14 * feduc - 0.03 * dist + 1.1 * ability + rnorm(n, 0, 1.6), 6), 20))
exper <- pmax(0, round(rnorm(n, 16, 8)))
fsize <- round(rnorm(n, 4, 1.2), 2); exporter <- rbinom(n, 1, 0.3); capital <- round(rnorm(n, 2, 1), 2)
q <- 0.06 * fsize + 0.12 * exporter + 0.05 * capital + rnorm(n, 0, 0.08)
sdw <- 0.14 * exp(0.15 * (educ - 8))
lwage <- 0.9 + 0.08 * educ + 0.06 * exper - 0.0012 * exper^2 + 0.5 * q + 0.07 * ability + sdw * rnorm(n)
workers <- data.frame(id = 1:n, lwage = round(lwage, 4), educ, exper, meduc, feduc, dist, sibs,
                      size = fsize, exporter, capital)
write.csv(workers, file.path(out, "workers.csv"), row.names = FALSE)
m <- 100
fs <- round(rnorm(m, 4, 1.2), 2); fe <- rbinom(m, 1, 0.3); fk <- round(rnorm(m, 2, 1), 2)
fq <- 0.06 * fs + 0.12 * fe + 0.05 * fk + rnorm(m, 0, 0.08) + rnorm(m, 0, 0.15)
write.csv(data.frame(firm = 1:m, size = fs, exporter = fe, capital = fk, q = round(fq, 4)),
          file.path(out, "firms.csv"), row.names = FALSE)

# ---- Unit 2: idealista_alicante_clase.csv and idealista_oraculo_clase.csv ----
n <- 1000
zonas <- c("Benalúa", "Centro", "Carolinas", "San Blas", "Playa de San Juan")
zona <- sample(zonas, n, replace = TRUE, prob = c(.2, .22, .22, .2, .16))
anio <- sample(2013:2026, n, replace = TRUE); mes <- sample(1:12, n, replace = TRUE)
calidad <- rnorm(n)
m2 <- round(pmax(35, exp(rnorm(n, log(86), 0.33) + 0.15 * (zona == "Playa de San Juan"))))
m2_anuncio <- round(m2 + rnorm(n, 4, 11.6))
terraza <- as.integer(runif(n) < plogis(-0.4 + 0.9 * calidad + 0.8 * (zona == "Playa de San Juan")))
m2_terraza <- terraza * round(pmax(3, rnorm(n, 13.4, 6)))
antig <- pmax(0, round(rnorm(n, 38, 20) - 12 * (zona == "Playa de San Juan")))
bajo <- rbinom(n, 1, 0.15)
zeff <- c(0, 0.056, -0.12, -0.175, 0.32)[match(zona, zonas)]
yeff <- cumsum(c(0, rnorm(13, 0.045, 0.03)))[anio - 2012]
sdu <- ifelse(zona == "Playa de San Juan", 0.2, 0.145)
lprecio <- 7.6 + 0.9 * log(m2) + 0.03 * terraza + 0.0075 * m2_terraza - 0.0142 * antig +
  0.00013 * antig^2 - 0.12 * bajo + zeff + yeff + 0.0975 * calidad + sdu * rnorm(n)
cert_energ <- pmin(7, pmax(1, round(4.8 - 0.045 * antig + 0.5 * calidad + rnorm(n))))
num_fotos <- pmax(3, round(18 + 6 * calidad + rnorm(n, 0, 5.5)))
d <- data.frame(id = 1:n, precio = round(exp(lprecio), -2), m2, m2_anuncio, terraza, m2_terraza,
                anio_construccion = anio - antig, antig, bajo, zona,
                fecha = sprintf("%d-%02d", anio, mes), anio, mes, cert_energ, num_fotos)
write.csv(d, file.path(out, "idealista_alicante_clase.csv"), row.names = FALSE)
write.csv(data.frame(id = 1:n, calidad), file.path(out, "idealista_oraculo_clase.csv"), row.names = FALSE)
cat("Mock data written to", out, "\n")
