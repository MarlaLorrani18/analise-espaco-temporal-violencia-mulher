# ==============================================================================
# FIGURAS DO TCC (Figuras 1 a 4 do texto)
# Rode DEPOIS de tcc_viol_mulher_analise.R, na mesma pasta de projeto: este
# script lê as tabelas que ele grava em resultados_R/tabelas.
# Saída: resultados_R/figuras_tcc/fig1.png ... fig4.png (300 dpi)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(sf); library(tidyr)
})
sf::sf_use_s2(FALSE)

DIR_PROJETO <- getwd()
ARQ_MALHA_LOCAL <- NULL   # NULL = malha das UFs do IBGE (2020) pelo geobr
TAB <- file.path(DIR_PROJETO, "resultados_R", "tabelas")
OUT <- file.path(DIR_PROJETO, "resultados_R", "figuras_tcc")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ANOS_PRELIMINARES <- 2025

ler <- function(f) read.csv2(file.path(TAB, f), check.names = FALSE, encoding = "UTF-8")
dec <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")
lab_ano <- function(a) ifelse(a %in% ANOS_PRELIMINARES, paste0(a, "*"), a)

tema <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        axis.text = element_text(colour = "grey25"), axis.title = element_text(colour = "grey15"),
        plot.tag = element_text(face = "bold"), legend.position = "none")

tx <- ler("taxas_unidade_ano.csv")

# --- Figura 1A: taxa no Brasil e nas regiões ------------------------------------
cores <- c(Brasil = "#1f1f1f", Norte = "#2a78d6", Nordeste = "#eb6834", Sudeste = "#1baf7a",
           Sul = "#eda100", `Centro-Oeste` = "#e87ba4")
d1 <- tx |> filter(nivel != "UF") |> mutate(unidade = ifelse(unidade == "BR", "Brasil", unidade))
fim <- d1 |> filter(ano == max(ano)) |> arrange(taxa)
fim$y_lab <- fim$taxa
for (i in 2:nrow(fim)) if (fim$y_lab[i] - fim$y_lab[i - 1] < 22) fim$y_lab[i] <- fim$y_lab[i - 1] + 22
gA <- ggplot(d1, aes(ano, taxa, colour = unidade)) +
  geom_line(aes(linewidth = unidade == "Brasil")) + geom_point(size = 0.9) +
  geom_text(data = fim, aes(x = max(d1$ano) + 0.4, y = y_lab, label = unidade), hjust = 0, size = 2.6, colour = "grey15") +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.6)) + scale_colour_manual(values = cores) +
  scale_x_continuous(breaks = seq(2009, 2025, 2), labels = lab_ano, expand = expansion(mult = c(0.01, 0.2))) +
  scale_y_continuous(labels = function(x) dec(x), limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Ano de notificação", y = "Taxa por 100 mil mulheres", tag = "A") + tema

# --- Figura 1B: municípios notificantes -----------------------------------------
mun <- ler("tab_municipios_notificantes.csv") |> group_by(ano) |> summarise(m = sum(municipios_notificantes))
gB <- ggplot(mun, aes(ano, m)) + geom_col(fill = "#2a78d6", width = 0.7) +
  geom_text(data = mun |> filter(ano %in% c(2009, 2013, 2019, 2025)), aes(label = dec(m)),
            vjust = -0.4, size = 2.5, colour = "grey15") +
  geom_hline(yintercept = 5570, linetype = 2, colour = "grey45", linewidth = 0.4) +
  annotate("text", x = 2009, y = 5570, label = "Total de municípios (5.570)", hjust = 0, vjust = -0.5,
           size = 2.5, colour = "grey30") +
  scale_x_continuous(breaks = seq(2009, 2025, 2), labels = lab_ano) +
  scale_y_continuous(labels = function(x) dec(x), limits = c(0, 6200), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Ano de notificação", y = "Municípios notificantes", tag = "B") + tema

png(file.path(OUT, "fig1.png"), width = 16, height = 13, units = "cm", res = 300)
gridExtra::grid.arrange(gA, gB, ncol = 1, heights = c(1.15, 1))
dev.off()

# --- Figura 2: mapas das taxas em anos selecionados -----------------------------
malha <- if (is.null(ARQ_MALHA_LOCAL)) {
  geobr::read_state(code_state = "all", year = 2020, showProgress = FALSE) |>
    transmute(sigla_uf = abbrev_state)
} else {
  m <- sf::st_read(ARQ_MALHA_LOCAL, quiet = TRUE)
  col <- intersect(c("abbrev_state", "SIGLA_UF", "sigla", "SIGLA"), names(m))[1]
  m |> transmute(sigla_uf = .data[[col]])
}
sel <- c(2009, 2013, 2017, 2021, 2025)
classes <- c("< 50", "50 a < 100", "100 a < 200", "200 a < 300", "≥ 300")
mp <- malha |>
  inner_join(tx |> filter(nivel == "UF", ano %in% sel), by = c("sigla_uf" = "unidade")) |>
  mutate(classe = factor(cut(taxa, c(0, 50, 100, 200, 300, Inf), right = FALSE, labels = classes), levels = classes),
         ano = factor(lab_ano(ano), lab_ano(sel)))
stopifnot(length(unique(mp$sigla_uf)) == 27)
g2 <- ggplot(mp) + geom_sf(aes(fill = classe), colour = "white", linewidth = 0.15) +
  facet_wrap(~ ano, ncol = 3) +
  scale_fill_manual(values = c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"), drop = FALSE,
                    name = "Taxa por 100 mil mulheres") +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  theme_void(base_size = 9) +
  theme(legend.position = c(0.84, 0.25), strip.text = element_text(face = "bold", size = 9))
ggsave(file.path(OUT, "fig2.png"), g2, width = 16, height = 11, units = "cm", dpi = 300, bg = "white")

# --- Figura 3: índice global de Moran -------------------------------------------
mo <- ler("moran_global.csv") |>
  mutate(matriz = recode(matriz, queen = "Contiguidade (queen)", knn3 = "3 vizinhos mais próximos",
                         knn2 = "2 vizinhos mais próximos"),
         matriz = factor(matriz, c("Contiguidade (queen)", "3 vizinhos mais próximos", "2 vizinhos mais próximos")))
g3 <- ggplot(mo, aes(ano, I, colour = matriz, shape = matriz)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c(16, 17, 21)) + scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a")) +
  scale_x_continuous(breaks = seq(2009, 2025, 2), labels = lab_ano) +
  scale_y_continuous(labels = function(x) dec(x, 1)) +
  labs(x = "Ano de notificação", y = "Índice global de Moran (I)", colour = NULL, shape = NULL) +
  tema + theme(legend.position = "top")
ggsave(file.path(OUT, "fig3.png"), g3, width = 16, height = 7.5, units = "cm", dpi = 300, bg = "white")

# --- Figura 4: APC por UF --------------------------------------------------------
regioes <- c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste")
pw <- ler("prais_winsten_completo.csv") |>
  filter(modelo == "linear", !unidade %in% regioes, periodo %in% c("2009-2025", "2009-2019 (pré-pandemia)")) |>
  mutate(periodo = factor(ifelse(periodo == "2009-2025", "2009–2025", "2009–2019"), c("2009–2025", "2009–2019")),
         unidade = recode(unidade, BR = "Brasil"))
ordem <- pw |> filter(periodo == "2009–2025") |> arrange(APC) |> pull(unidade)
pw$unidade <- factor(pw$unidade, ordem)
g4 <- ggplot(pw, aes(APC, unidade, colour = periodo, shape = periodo)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = APC_LI, xmax = APC_LS), width = 0, linewidth = 0.5, orientation = "y",
                position = position_dodge(width = 0.6)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834")) + scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(labels = function(x) dec(x)) +
  labs(x = "Variação percentual anual (APC, %) e IC 95%", y = NULL, colour = "Período", shape = "Período") +
  tema + theme(legend.position = "top", panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
               panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7))
ggsave(file.path(OUT, "fig4.png"), g4, width = 16, height = 15, units = "cm", dpi = 300, bg = "white")

message("Figuras salvas em: ", OUT)
