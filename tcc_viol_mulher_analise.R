# ==============================================================================
# NOTIFICAÇÕES DE VIOLÊNCIA CONTRA MULHERES NO BRASIL (SINAN/VIVA), 2009–2025
# Análise descritiva, taxas, autocorrelação espacial (Moran global, diagrama
# de espalhamento de Moran, LISA) e tendência temporal (Prais-Winsten, APC)
#
# TCC — MBA em Data Science & Analytics (USP/Esalq)
# Autora: Marla Lorrani Santos de Souza
# Referência metodológica: TCC II (2023), com as correções apontadas no
# parecer de pré-avaliação de 06/10/2026.
#
# Como usar: ajuste apenas a seção 1 (caminhos e parâmetros) e rode o script
# inteiro (Ctrl+Shift+Enter no RStudio). Todas as saídas vão para DIR_SAIDA.
# ==============================================================================
#
# ------------------------------------------------------------------------------
# DIAGNÓSTICO: MAPEAMENTO ENTRE O TCC II E AS NOVAS BASES
# ------------------------------------------------------------------------------
# Etapa do TCC II       | No TCC II                       | Nas novas bases                              | Adaptação
# ----------------------|---------------------------------|----------------------------------------------|------------------------------------------
# Fonte do numerador    | SINAN violência (VIVA), 2009–21 | VIOLBR09–VIOLBR25.dbc (DATASUS), 2009–25     | Leitura arquivo a arquivo; auditoria de
#                       |                                 |                                              | estrutura, pois variáveis mudam entre fichas
# Período (tempo)       | Ano da notificação              | NU_ANO (ano da notificação); o arquivo é     | NU_ANO conferido contra o ano do arquivo;
#                       |                                 | organizado por ano de notificação            | 2025 preliminar no DATASUS (marcado)     
# Sexo                  | Filtro sexo = feminino          | CS_SEXO ("F", "M", "I")                      | Mesmo significado; "I" (ignorado) excluído e
#                       |                                 |                                              | quantificado
# Autoprovocada         | Excluída                        | LES_AUTOP (1 = sim, 2 = não, 9 = ignorado)   | Exclui só LES_AUTOP = 1; ignorados mantidos
#                       |                                 |                                              | e quantificados por ano
# Idade                 | Sem restrição (todas as idades) | NU_IDADE_N (código: 1º dígito = unidade,     | Decodificada para anos; faixa definida em
#                       |                                 | 4 = anos; demais = < 1 ano)                  | IDADE_MIN/IDADE_MAX (padrão: todas)
# Unidade espacial      | UF (27)                         | SG_UF (residência), SG_UF_NOT (notificação), | UF de residência por padrão (o denominador
#                       |                                 | SG_UF_OCOR (ocorrência), códigos IBGE        | é população residente); alterável
# Denominador           | População feminina por UF       | Projeções IBGE 2024, idade simples, sexo,    | Soma de "Mulheres" na faixa etária do
#                       |                                 | UF/região/Brasil, 2000–2025 (90 = 90+)       | numerador, por UF e ano
# Taxa                  | casos / pop. fem. × 100 mil     | Idem                                         | Unidade: "por 100 mil mulheres"
# Perfil (descritiva)   | Raça/cor, escolaridade etc.     | CS_RACA, CS_ESCOL_N, SIT_CONJUG, LOCAL_OCOR, | Recodificadas com rótulos; "Ignorado"
#                       |                                 | OUT_VEZES, VIOL_*, REL_*, AUTOR_SEXO,        | separado e completude por ano;
#                       |                                 | AUTOR_ALCO, NUM_ENVOLV                       | pretas + pardas = negras (IBGE)
# Normalização min–max  | Contagens por UF normalizadas   | Mesmas variáveis                             | NÃO reproduzida: min–max de contagem mede
#                       | para mapas por característica   |                                              | tamanho da população, não risco. Mapas por
#                       |                                 |                                              | característica usam proporções na UF
# Moran global          | kNN k = 3 e k = 2               | Malha de UFs (geobr/IBGE)                    | Contiguidade queen (principal) + kNN k = 3
#                       |                                 |                                              | e k = 2 (reprodução do TCC/sensibilidade);
#                       |                                 |                                              | teste analítico e por permutação (999)
# Diagrama de Moran     | Quadrantes AA, BB, AB, BA       | Taxa padronizada × defasagem espacial        | Igual; por ano, com mapa de quadrantes
# LISA                  | Mapas de clusters               | Idem                                         | Permutação condicional; p bruto e FDR
# Prais-Winsten         | log da taxa ~ ano, por UF;      | 17 anos (2009–2025)                          | ln(taxa) ~ ano; períodos 2009–2025,
#                       | 2009–2021 e 2009–2019           |                                              | 2009–2019 e 2009–2024 (sem preliminar);  
#                       |                                 |                                              | sensibilidade com indicador 2020–2021
# APC                   | Texto: −1 + 10^β; tabelas: e^β  | —                                            | APC = (e^β − 1) × 100 (β em ln), IC com
#                       | (INCONSISTENTE)                 |                                              | quantil t de n − 2 g.l.
# Classificação         | Crescente/estável, critério não | —                                            | Critério único: IC 95% da APC exclui 0
#                       | declarado (contradizia p-valor) |                                              | (equivalente a p < 0,05 no mesmo teste t)
# Durbin-Watson         | Não informado se antes/depois   | —                                            | DW original e transformado, + teste DW (OLS)
# ------------------------------------------------------------------------------
#
# AUDITORIA DAS BASES REAIS (VIOLBR09–25, lidas em 06/10/2026)
#  - 17 arquivos, 5.335.262 notificações; as 159 variáveis são as mesmas, com
#    os mesmos nomes, em todos os anos (ficha unificada; campos novos da ficha
#    2015, como CICL_VID, ORIENT_SEX e IDENT_GEN, aparecem vazios antes).
#  - TP_NOT = 2 e ID_AGRAVO = Y09 em todas as linhas (só violência individual).
#  - NU_ANO vazio em todo o arquivo de 2020 (e em 60 linhas de 2019); DT_NOTIFIC
#    coincide com o ano do arquivo em 100% das linhas -> ano vem de DT_NOTIFIC.
#  - Texto livre em latin1 -> convertido para UTF-8 na leitura.
#  - NDUPLIC e TPUNINOT 100% vazios: o SINAN público não traz chave para
#    deduplicação. Duplicatas exatas: 1.418 linhas (0,03%), apenas reportadas.
#  - NU_IDADE_N com unidade 0 (ex.: 0000) = idade ignorada (~0,1%).
#  - Códigos fora do dicionário em quantidade desprezível (AUTOR_SEXO 4/5,
#    NUM_ENVOLV 3, LES_AUTOP 0) -> tratados como ignorado.
#  - LES_AUTOP vazio/ignorado em 257.306 notificações de mulheres (mantidas),
#    com salto em 2025 (dado preliminar).
#  - SG_UF (residência) ausente em 389 notificações incluídas (0,01%);
#    SG_UF_OCOR tem muito mais ausentes (até 50 mil em 2025) -> residência é a
#    UF de referência mais completa.
#  - Base analítica (mulheres, exceto autoprovocada, todas as idades):
#    2.814.808 notificações.
#
# INCONSISTÊNCIAS DO TCC II SINALIZADAS E CORRIGIDAS AQUI
#  (1) A fórmula da APC no texto (−1 + 10^β) não é a usada nas tabelas (e^β).
#      Aqui: regressão em logaritmo natural e APC = (exp(β) − 1) × 100.
#  (2) O IC usava 1,96 (normal) e o p-valor, t com n − 2 g.l.; por isso
#      discordavam. Aqui ambos usam a mesma distribuição t.
#  (3) A classificação "Crescente/Estável" contradizia p-valores. Aqui o
#      critério é único e declarado.
#  (4) O número de vizinhos (k = 3 / k = 2) foi escolhido sem justificativa.
#      Aqui a matriz principal (contiguidade queen) é pré-definida, e kNN
#      entra como sensibilidade.
#  (5) Taxas mostram NOTIFICAÇÕES, não ocorrências: a expansão do VIVA após
#      2009 confunde tendência de notificação com tendência de violência.
#      O script calcula o nº de municípios notificantes por UF e ano para
#      apoiar essa discussão.
# ==============================================================================


# ============================================================
# 1. CONFIGURAÇÃO DO AMBIENTE
# ============================================================

rm(list = ls()); invisible(gc())
options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)

# --- Caminhos ------------------------------------------------------------------
# Pasta do projeto: a que contém as subpastas bases_viol_2009-2025 e
# bases_estimativa_pop_fem. Por padrão, a pasta de trabalho atual (abra o
# projeto do RStudio nessa pasta ou use setwd()).
DIR_PROJETO <- getwd()
DIR_SINAN <- file.path(DIR_PROJETO, "bases_viol_2009-2025")
DIR_POP   <- file.path(DIR_PROJETO, "bases_estimativa_pop_fem")
DIR_SAIDA <- file.path(DIR_PROJETO, "resultados_R")

ARQ_POP <- NULL          # NULL = primeiro .xlsx encontrado em DIR_POP
ARQ_MALHA_LOCAL <- NULL  # NULL = baixa a malha de UFs pelo geobr. Se não houver
                         # internet, informe um shapefile/geojson de UFs do IBGE,
                         # ex.: "C:/.../BR_UF_2022.shp"

# --- Parâmetros do estudo ------------------------------------------------------
ANOS_ESTUDO       <- 2009:2025
ANOS_PRELIMINARES <- 2025  # só 2025 está como preliminar no DATASUS (conferido pela autora)
DATA_DOWNLOAD     <- "2026-09-24"   # data de download dos .dbc (para a Metodologia)

IDADE_MIN <- 0      # idade mínima em anos (0 = sem restrição, como no TCC II)
IDADE_MAX <- 200    # idade máxima em anos (200 = sem restrição)
# Atenção: se restringir a idade (ex.: 18+), o denominador é restringido igual.

EXCLUIR_AUTOPROVOCADA     <- TRUE          # como no TCC II
UF_REFERENCIA             <- "residencia"  # "residencia", "notificacao" ou "ocorrencia"
REMOVER_DUPLICATAS_EXATAS <- FALSE         # FALSE = apenas reporta duplicatas
USAR_CACHE                <- TRUE          # reaproveita .rds já processados

FATOR_TAXA        <- 1e5     # taxa por 100 mil mulheres
ALPHA             <- 0.05
N_PERMUTACOES     <- 999
SEMENTE           <- 20261006
ALTERNATIVA_MORAN <- "two.sided"  # H1: I != 0 ("greater" = só autocorrelação positiva)
MATRIZ_PRINCIPAL  <- "queen"      # "queen", "knn3" ou "knn2"
K_VIZINHOS        <- c(3, 2)      # kNN usados no TCC II
SNAP_CONTIGUIDADE <- 0.01         # tolerância (graus) para fronteiras da malha simplificada

PERIODOS_TENDENCIA <- list(
  "2009-2025"                    = 2009:2025,
  "2009-2019 (pré-pandemia)"     = 2009:2019,  # como no TCC II
  "2009-2024 (sem preliminar)" = 2009:2024
)
PERIODO_PRINCIPAL  <- "2009-2025"
ANOS_PANDEMIA      <- 2020:2021   # indicador da análise de sensibilidade
MIN_ANOS_TENDENCIA <- 7           # mínimo de pontos para ajustar Prais-Winsten

FONTE_GRAF <- "Fonte: SINAN/DATASUS (notificações de violência); IBGE, Projeções da População 2024."

set.seed(SEMENTE)

DIR_CACHE <- file.path(DIR_SAIDA, "cache")
DIR_TAB   <- file.path(DIR_SAIDA, "tabelas")
DIR_FIG   <- file.path(DIR_SAIDA, "figuras")
DIR_AUD   <- file.path(DIR_SAIDA, "auditoria")
for (d in c(DIR_SAIDA, DIR_CACHE, DIR_TAB, DIR_FIG, DIR_AUD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

stopifnot(UF_REFERENCIA %in% c("residencia", "notificacao", "ocorrencia"),
          MATRIZ_PRINCIPAL %in% c("queen", paste0("knn", K_VIZINHOS)),
          PERIODO_PRINCIPAL %in% names(PERIODOS_TENDENCIA),
          IDADE_MIN <= IDADE_MAX)


# ============================================================
# 2. PACOTES
# ============================================================
# read.dbc: leitura dos .dbc do DATASUS. Se não estiver no CRAN para sua
# versão do R, instala do GitHub (danicat/read.dbc; no Windows exige Rtools).

instalar_se_faltar <- function(pkgs) {
  falta <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(falta)) install.packages(falta)
}
instalar_se_faltar(c("foreign", "readxl", "dplyr", "tidyr", "purrr", "stringr",
                     "ggplot2", "scales", "sf", "spdep", "prais", "lmtest",
                     "writexl", "vctrs"))
if (!requireNamespace("read.dbc", quietly = TRUE)) {
  try(install.packages("read.dbc"), silent = TRUE)
  if (!requireNamespace("read.dbc", quietly = TRUE)) {
    instalar_se_faltar("remotes")
    remotes::install_github("danicat/read.dbc")
  }
}
if (is.null(ARQ_MALHA_LOCAL)) instalar_se_faltar("geobr")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(sf); library(spdep)
})

# --- Funções auxiliares gerais -------------------------------------------------
msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
aviso <- function(...) {
  txt <- sprintf(...)
  AVISOS[[length(AVISOS) + 1]] <<- txt
  message("  AVISO: ", txt)
}
AVISOS <- list()

num_br  <- scales::label_number(big.mark = ".", decimal.mark = ",")
num_br1 <- scales::label_number(accuracy = 0.1, big.mark = ".", decimal.mark = ",")
pct_br  <- scales::label_percent(accuracy = 0.1, big.mark = ".", decimal.mark = ",")
fmt2    <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")

tema_tcc <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
tema_mapa <- theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

salvar_fig <- function(p, nome, w = 8, h = 5) {
  ggsave(file.path(DIR_FIG, paste0(nome, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  invisible(p)
}
salvar_csv <- function(x, nome, dir = DIR_TAB) {
  write.csv2(x, file.path(dir, paste0(nome, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  invisible(x)
}

# --- Tabela de referência das UFs (códigos IBGE) -------------------------------
UFS <- data.frame(
  cod_uf = c(11,12,13,14,15,16,17,21,22,23,24,25,26,27,28,29,31,32,33,35,41,42,43,50,51,52,53),
  sigla_uf = c("RO","AC","AM","RR","PA","AP","TO","MA","PI","CE","RN","PB","PE","AL","SE","BA",
               "MG","ES","RJ","SP","PR","SC","RS","MS","MT","GO","DF"),
  nome_uf = c("Rondônia","Acre","Amazonas","Roraima","Pará","Amapá","Tocantins","Maranhão","Piauí",
              "Ceará","Rio Grande do Norte","Paraíba","Pernambuco","Alagoas","Sergipe","Bahia",
              "Minas Gerais","Espírito Santo","Rio de Janeiro","São Paulo","Paraná","Santa Catarina",
              "Rio Grande do Sul","Mato Grosso do Sul","Mato Grosso","Goiás","Distrito Federal")
)
REGIOES <- data.frame(cod_regiao = 1:5, regiao = c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste"))
UFS$cod_regiao <- UFS$cod_uf %/% 10
UFS <- merge(UFS, REGIOES, by = "cod_regiao")
UFS <- UFS[order(UFS$cod_uf), ]


# ============================================================
# 3. IMPORTAÇÃO DAS BASES
# ============================================================
# Cada arquivo é lido, auditado, filtrado e salvo em .rds; só o subconjunto
# analítico fica na memória (os 17 arquivos completos somam vários GB).

# Variáveis de interesse (as ausentes em algum ano são criadas como NA e
# registradas na auditoria de estrutura).
VARS_VIOL <- c("VIOL_FISIC","VIOL_PSICO","VIOL_TORT","VIOL_SEXU","VIOL_TRAF","VIOL_FINAN",
               "VIOL_NEGLI","VIOL_INFAN","VIOL_LEGAL","VIOL_OUTR")
VARS_REL  <- c("REL_PAI","REL_MAE","REL_PAD","REL_MAD","REL_CONJ","REL_EXCON","REL_NAMO",
               "REL_EXNAM","REL_FILHO","REL_DESCO","REL_IRMAO","REL_CONHEC","REL_CUIDA",
               "REL_PATRAO","REL_INST","REL_POL","REL_PROPRI","REL_OUTROS")
VARS_ID   <- c("NU_ANO","DT_NOTIFIC","DT_OCOR","SG_UF_NOT","ID_MUNICIP","SG_UF","ID_MN_RESI",
               "SG_UF_OCOR","ID_MN_OCOR")
VARS_PERFIL <- c("CS_SEXO","NU_IDADE_N","CS_GESTANT","CS_RACA","CS_ESCOL_N","SIT_CONJUG",
                 "DEF_TRANS","LOCAL_OCOR","OUT_VEZES","LES_AUTOP","AUTOR_SEXO","AUTOR_ALCO",
                 "NUM_ENVOLV","CICL_VID")
VARS_ANALISE <- c(VARS_ID, VARS_PERFIL, VARS_VIOL, VARS_REL)
# Variáveis categóricas cujas categorias brutas são tabuladas na auditoria
VARS_CATEG_AUDIT <- c("NU_ANO","SG_UF","SG_UF_NOT","SG_UF_OCOR","CS_SEXO","CS_GESTANT","CS_RACA",
                      "CS_ESCOL_N","SIT_CONJUG","DEF_TRANS","LOCAL_OCOR","OUT_VEZES","LES_AUTOP",
                      "AUTOR_SEXO","AUTOR_ALCO","NUM_ENVOLV","CICL_VID", VARS_VIOL, VARS_REL)

arquivos <- list.files(DIR_SINAN, pattern = "^VIOLBR\\d{2}\\.(dbc|DBC|dbf|DBF)$", full.names = TRUE)
if (!length(arquivos)) stop("Nenhum arquivo VIOLBRxx.dbc encontrado em: ", DIR_SINAN)
ano_dos_arquivos <- 2000 + as.integer(stringr::str_extract(basename(arquivos), "\\d{2}"))
arquivos <- arquivos[order(ano_dos_arquivos)]; ano_dos_arquivos <- sort(ano_dos_arquivos)
msg("%d arquivos encontrados: %s", length(arquivos), paste(basename(arquivos), collapse = ", "))
falta_ano <- setdiff(ANOS_ESTUDO, ano_dos_arquivos)
if (length(falta_ano)) aviso("Sem arquivo para os anos: %s", paste(falta_ano, collapse = ", "))

ler_arquivo_sinan <- function(f) {
  ext <- tolower(tools::file_ext(f))
  df <- if (ext == "dbc") read.dbc::read.dbc(f, as.is = TRUE) else foreign::read.dbf(f, as.is = TRUE)
  names(df) <- toupper(trimws(names(df)))  # padroniza nomes
  # Os .dbc do DATASUS vêm em latin1 (campos de texto livre, como *_ESPEC);
  # converte para UTF-8 para evitar erro de "invalid UTF-8"
  chr <- vapply(df, is.character, logical(1))
  df[chr] <- lapply(df[chr], iconv, from = "latin1", to = "UTF-8")
  df
}

# Converte data em formatos "AAAA-MM-DD", "AAAAMMDD" ou Date
para_data <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- trimws(as.character(x)); x[x == ""] <- NA
  out <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  i <- is.na(out) & !is.na(x)
  out[i] <- suppressWarnings(as.Date(x[i], format = "%Y%m%d"))
  out
}
# Converte código categórico ("01", "1", " 1") para inteiro
para_cod <- function(x) suppressWarnings(as.integer(trimws(as.character(x))))

# Decodifica NU_IDADE_N: 1º dígito = unidade (1 hora, 2 dia, 3 mês, 4 ano)
decodifica_idade <- function(x) {
  x <- para_cod(x)
  unid <- x %/% 1000; val <- x %% 1000
  idade <- ifelse(unid == 4, val, ifelse(unid %in% 1:3, 0L, NA_integer_))
  idade[!is.na(idade) & idade > 120] <- NA_integer_
  as.integer(idade)
}

# --- Auditoria de um arquivo bruto ---------------------------------------------
auditar_arquivo <- function(df, arquivo, ano_arq) {
  n <- nrow(df)
  n_unicos <- length(vctrs::vec_unique_loc(df))
  geral <- data.frame(
    arquivo = arquivo, formato = toupper(tools::file_ext(arquivo)), ano_arquivo = ano_arq,
    n_linhas = n, n_colunas = ncol(df), n_duplicatas_exatas = n - n_unicos,
    n_NU_ANO_ausente = if ("NU_ANO" %in% names(df)) sum(is.na(df$NU_ANO) | trimws(df$NU_ANO) == "") else n,
    pct_NU_ANO_igual_arquivo = if ("NU_ANO" %in% names(df) && any(!is.na(df$NU_ANO))) round(100 * mean(para_cod(df$NU_ANO) == ano_arq, na.rm = TRUE), 2) else NA,
    pct_DT_NOTIFIC_igual_arquivo = if ("DT_NOTIFIC" %in% names(df)) round(100 * mean(as.integer(format(para_data(df$DT_NOTIFIC), "%Y")) == ano_arq, na.rm = TRUE), 2) else NA
  )
  variaveis <- do.call(rbind, lapply(names(df), function(v) {
    x <- df[[v]]
    xc <- if (is.character(x)) trimws(x) else x
    n_vazio <- if (is.character(xc)) sum(xc == "", na.rm = TRUE) else 0L
    u <- unique(xc)
    cats <- if (length(u) <= 30) {
      tb <- sort(table(xc, useNA = "ifany"), decreasing = TRUE)
      paste0(names(tb), " (", tb, ")", collapse = "; ")
    } else paste0("> 30 categorias; ex.: ", paste(head(u[!is.na(u)], 5), collapse = "; "))
    data.frame(arquivo = arquivo, ano_arquivo = ano_arq, variavel = v, tipo = class(x)[1],
               n_ausentes = sum(is.na(x)), n_vazios = n_vazio,
               pct_ausente_ou_vazio = round(100 * (sum(is.na(x)) + n_vazio) / n, 2),
               n_categorias = length(u), categorias = cats)
  }))
  cat_vars <- intersect(VARS_CATEG_AUDIT, names(df))
  categorias <- do.call(rbind, lapply(cat_vars, function(v) {
    tb <- table(trimws(as.character(df[[v]])), useNA = "ifany")
    data.frame(ano_arquivo = ano_arq, variavel = v, valor = names(tb), n = as.integer(tb))
  }))
  list(geral = geral, variaveis = variaveis, categorias = categorias,
       nomes = data.frame(ano_arquivo = ano_arq, variavel = names(df)))
}

# --- Processamento de um arquivo: auditoria + seleção + filtros ----------------
processar_arquivo <- function(f, ano_arq) {
  rds_base <- file.path(DIR_CACHE, sprintf("viol_%d.rds", ano_arq))
  rds_aud  <- file.path(DIR_CACHE, sprintf("auditoria_%d.rds", ano_arq))
  if (USAR_CACHE && file.exists(rds_base) && file.exists(rds_aud)) {
    msg("  %s: usando cache", basename(f))
    return(list(base = readRDS(rds_base), auditoria = readRDS(rds_aud)))
  }
  msg("  Lendo %s ...", basename(f))
  df <- tryCatch(ler_arquivo_sinan(f), error = function(e) {
    aviso("Falha ao ler %s: %s", basename(f), conditionMessage(e)); NULL })
  if (is.null(df)) return(NULL)
  aud <- auditar_arquivo(df, basename(f), ano_arq)

  # Seleção de variáveis (cria NA para as ausentes)
  ausentes <- setdiff(VARS_ANALISE, names(df))
  for (v in ausentes) df[[v]] <- NA
  df <- df[, VARS_ANALISE]

  # Padronização de tipos
  cod_vars <- setdiff(VARS_ANALISE, c("DT_NOTIFIC", "DT_OCOR", "CS_SEXO", "ID_MUNICIP",
                                      "ID_MN_RESI", "ID_MN_OCOR", "NU_IDADE_N"))
  for (v in cod_vars) df[[v]] <- para_cod(df[[v]])
  for (v in c("ID_MUNICIP", "ID_MN_RESI", "ID_MN_OCOR")) df[[v]] <- trimws(as.character(df[[v]]))
  df$DT_NOTIFIC <- para_data(df$DT_NOTIFIC)
  df$DT_OCOR    <- para_data(df$DT_OCOR)
  df$CS_SEXO    <- toupper(trimws(as.character(df$CS_SEXO)))
  df$idade_anos <- decodifica_idade(df$NU_IDADE_N)
  df$NU_IDADE_N <- NULL

  # Ano da notificação: ano de DT_NOTIFIC; se ausente, NU_ANO; se ausente, ano do arquivo.
  # (Nas bases reais, NU_ANO está vazio em todo o arquivo de 2020 e em 60 linhas
  # de 2019, enquanto DT_NOTIFIC coincide com o ano do arquivo em 100% das linhas.)
  df$ano <- dplyr::coalesce(as.integer(format(df$DT_NOTIFIC, "%Y")), df$NU_ANO, as.integer(ano_arq))
  df$ano_arquivo <- as.integer(ano_arq)

  # UF de residência, notificação e ocorrência (com recuperação pelo código do município)
  uf_mun <- function(m) para_cod(substr(m, 1, 2))
  df$uf_residencia  <- dplyr::coalesce(df$SG_UF, uf_mun(df$ID_MN_RESI))
  df$uf_notificacao <- dplyr::coalesce(df$SG_UF_NOT, uf_mun(df$ID_MUNICIP))
  df$uf_ocorrencia  <- dplyr::coalesce(df$SG_UF_OCOR, uf_mun(df$ID_MN_OCOR))
  df$cod_uf <- switch(UF_REFERENCIA, residencia = df$uf_residencia,
                      notificacao = df$uf_notificacao, ocorrencia = df$uf_ocorrencia)
  df$cod_uf[!df$cod_uf %in% UFS$cod_uf] <- NA_integer_

  # Fluxo de seleção (critérios de inclusão/exclusão)
  fluxo <- data.frame(ano_arquivo = integer(), ordem = integer(), etapa = character(), n = integer())
  reg <- function(etapa, n) fluxo <<- rbind(fluxo, data.frame(ano_arquivo = ano_arq, ordem = nrow(fluxo) + 1L,
                                                              etapa = etapa, n = n))
  reg("1. Notificações no arquivo", nrow(df))
  df <- df[df$ano %in% ANOS_ESTUDO, ];                       reg("2. Ano de notificação no período", nrow(df))
  reg("   excl.: sexo masculino", sum(df$CS_SEXO == "M", na.rm = TRUE))
  reg("   excl.: sexo ignorado/em branco", sum(!df$CS_SEXO %in% c("F", "M")))
  df <- df[df$CS_SEXO %in% "F", ];                            reg("3. Sexo feminino", nrow(df))
  reg("   (info) LES_AUTOP ignorado/em branco, mantidas", sum(is.na(df$LES_AUTOP) | df$LES_AUTOP == 9))
  if (EXCLUIR_AUTOPROVOCADA) {
    reg("   excl.: lesão autoprovocada (LES_AUTOP = 1)", sum(df$LES_AUTOP %in% 1))
    df <- df[!df$LES_AUTOP %in% 1, ];                         reg("4. Exceto lesão autoprovocada", nrow(df))
  }
  restringe_idade <- IDADE_MIN > 0 || IDADE_MAX < 120
  if (restringe_idade) {
    reg("   excl.: idade ignorada", sum(is.na(df$idade_anos)))
    reg("   excl.: idade fora da faixa", sum(!is.na(df$idade_anos) & (df$idade_anos < IDADE_MIN | df$idade_anos > IDADE_MAX)))
    df <- df[!is.na(df$idade_anos) & df$idade_anos >= IDADE_MIN & df$idade_anos <= IDADE_MAX, ]
    reg(sprintf("5. Idade entre %d e %d anos", IDADE_MIN, IDADE_MAX), nrow(df))
  } else {
    reg("   (info) idade ignorada, mantidas", sum(is.na(df$idade_anos)))
  }
  reg(sprintf("   (info) UF de %s ausente/inválida (só na taxa Brasil)", UF_REFERENCIA), sum(is.na(df$cod_uf)))
  reg("6. Base analítica", nrow(df))

  aud$fluxo <- fluxo
  aud$ausentes <- if (length(ausentes)) data.frame(ano_arquivo = ano_arq, variavel_ausente = ausentes) else NULL
  saveRDS(df, rds_base); saveRDS(aud, rds_aud)
  rm(list = setdiff(ls(), c("df", "aud"))); invisible(gc())
  list(base = df, auditoria = aud)
}

msg("Importando e auditando os arquivos do SINAN...")
proc <- Map(processar_arquivo, arquivos, ano_dos_arquivos)
proc <- proc[!vapply(proc, is.null, logical(1))]
if (!length(proc)) stop("Nenhum arquivo pôde ser lido.")


# ============================================================
# 4. INSPEÇÃO E AUDITORIA DOS DADOS
# ============================================================

juntar <- function(el) do.call(rbind, lapply(proc, function(p) p$auditoria[[el]]))
aud_geral      <- juntar("geral")
aud_variaveis  <- juntar("variaveis")
aud_categorias <- juntar("categorias")
aud_fluxo      <- juntar("fluxo")
aud_ausentes   <- juntar("ausentes")
rownames(aud_geral) <- NULL

# Estrutura: presença de cada variável por ano (mudanças de ficha ao longo do tempo)
aud_estrutura <- juntar("nomes") |>
  mutate(presente = "x") |>
  pivot_wider(names_from = ano_arquivo, values_from = presente, values_fill = "") |>
  arrange(variavel)
n_anos_arq <- length(proc)
var_instaveis <- aud_estrutura$variavel[rowSums(aud_estrutura[, -1] == "x") < n_anos_arq]
if (length(var_instaveis))
  aviso("%d variáveis não estão presentes em todos os anos (ver auditoria_estrutura.csv).", length(var_instaveis))
if (!is.null(aud_ausentes) && nrow(aud_ausentes)) {
  ausentes_resumo <- aud_ausentes |> group_by(variavel_ausente) |>
    summarise(anos = ifelse(n() == n_anos_arq, "todos os anos", paste(ano_arquivo, collapse = ", ")))
  aviso("Variáveis de análise ausentes em algum ano: %s",
        paste0(ausentes_resumo$variavel_ausente, " [", ausentes_resumo$anos, "]", collapse = "; "))
}

# Consistência entre NU_ANO e o ano do arquivo
inconsist_ano <- aud_geral |> filter(!is.na(pct_DT_NOTIFIC_igual_arquivo), pct_DT_NOTIFIC_igual_arquivo < 99)
if (nrow(inconsist_ano))
  aviso("Ano de DT_NOTIFIC difere do ano do arquivo em > 1%% das linhas em: %s", paste(inconsist_ano$arquivo, collapse = ", "))
nu_ano_vazio <- aud_geral |> filter(n_NU_ANO_ausente > 0)
if (nrow(nu_ano_vazio))
  aviso("NU_ANO vazio em: %s (ano obtido de DT_NOTIFIC).",
        paste0(nu_ano_vazio$arquivo, " [", nu_ano_vazio$n_NU_ANO_ausente, "]", collapse = ", "))
if (any(aud_geral$n_duplicatas_exatas > 0))
  aviso("Há linhas exatamente duplicadas em %d arquivo(s) (total %s). Ver auditoria_arquivos.csv.",
        sum(aud_geral$n_duplicatas_exatas > 0), sum(aud_geral$n_duplicatas_exatas))

# Fluxo de seleção em formato largo (linhas = etapas, colunas = anos)
tab_fluxo <- aud_fluxo |>
  group_by(etapa, ano_arquivo) |> summarise(n = sum(n), ordem = min(ordem)) |>
  group_by(etapa) |> mutate(ordem = min(ordem)) |> ungroup() |>
  pivot_wider(names_from = ano_arquivo, values_from = n, values_fill = 0) |>
  arrange(ordem) |> select(-ordem) |>
  mutate(Total = rowSums(across(where(is.numeric))))

salvar_csv(aud_geral, "auditoria_arquivos", DIR_AUD)
salvar_csv(aud_variaveis, "auditoria_variaveis", DIR_AUD)
salvar_csv(aud_categorias, "auditoria_categorias", DIR_AUD)
salvar_csv(aud_estrutura, "auditoria_estrutura", DIR_AUD)
salvar_csv(tab_fluxo, "fluxo_selecao", DIR_AUD)
msg("Auditoria salva em %s", DIR_AUD)
print(aud_geral[, c("arquivo", "n_linhas", "n_colunas", "n_duplicatas_exatas", "n_NU_ANO_ausente", "pct_DT_NOTIFIC_igual_arquivo")])

# --- Auditoria da projeção populacional ----------------------------------------
if (is.null(ARQ_POP)) {
  cand <- list.files(DIR_POP, pattern = "\\.xlsx$", full.names = TRUE)
  cand <- cand[!grepl("^~\\$", basename(cand))]
  if (!length(cand)) stop("Nenhum .xlsx de população em: ", DIR_POP)
  ARQ_POP <- cand[1]
}
msg("Lendo projeção populacional: %s", basename(ARQ_POP))
pop_bruta <- readxl::read_excel(ARQ_POP)
names(pop_bruta)[1:5] <- c("idade", "sexo", "cod_local", "sigla", "local")
col_anos_pop <- grep("^\\d{4}$", names(pop_bruta), value = TRUE)
aud_pop <- data.frame(
  arquivo = basename(ARQ_POP), n_linhas = nrow(pop_bruta), n_colunas = ncol(pop_bruta),
  anos = paste(range(as.integer(col_anos_pop)), collapse = "–"),
  sexos = paste(unique(pop_bruta$sexo), collapse = "; "),
  idades = paste(range(as.numeric(pop_bruta$idade)), collapse = "–"),
  n_locais = length(unique(pop_bruta$cod_local)),
  n_ausentes = sum(is.na(pop_bruta)),
  n_duplicatas = sum(duplicated(pop_bruta[, c("idade", "sexo", "cod_local")]))
)
salvar_csv(aud_pop, "auditoria_populacao", DIR_AUD)
falta_pop <- setdiff(ANOS_ESTUDO, as.integer(col_anos_pop))
if (length(falta_pop)) stop("Projeção sem os anos: ", paste(falta_pop, collapse = ", "))
if (IDADE_MAX >= 90 && IDADE_MAX < 120)
  aviso("A projeção agrega 90+ na idade 90; IDADE_MAX entre 90 e 119 usa todo o grupo 90+ no denominador.")


# ============================================================
# 5. LIMPEZA E PADRONIZAÇÃO
# ============================================================

base <- bind_rows(lapply(proc, `[[`, "base"))
rm(proc); invisible(gc())
msg("Base combinada: %s notificações", num_br(nrow(base)))

# Duplicatas exatas na base analítica (após seleção de variáveis)
dup_analitica <- nrow(base) - length(vctrs::vec_unique_loc(base))
msg("Duplicatas exatas nas variáveis analíticas: %s", num_br(dup_analitica))
if (REMOVER_DUPLICATAS_EXATAS && dup_analitica > 0) {
  base <- base[vctrs::vec_unique_loc(base), ]
  aviso("%s duplicatas exatas removidas.", num_br(dup_analitica))
}

# --- Rótulos (dicionário do SINAN, ficha de notificação individual de violência)
ROT <- list(
  CS_RACA    = c(`1` = "Branca", `2` = "Preta", `3` = "Amarela", `4` = "Parda", `5` = "Indígena"),
  CS_ESCOL_N = c(`0` = "Analfabeta", `1` = "Fundamental incompleto", `2` = "Fundamental incompleto",
                 `3` = "Fundamental incompleto", `4` = "Fundamental completo",
                 `5` = "Médio incompleto", `6` = "Médio completo", `7` = "Superior incompleto",
                 `8` = "Superior completo", `10` = "Não se aplica"),
  SIT_CONJUG = c(`1` = "Solteira", `2` = "Casada/união consensual", `3` = "Viúva",
                 `4` = "Separada", `8` = "Não se aplica"),
  LOCAL_OCOR = c(`1` = "Residência", `2` = "Habitação coletiva", `3` = "Escola",
                 `4` = "Local de prática esportiva", `5` = "Bar ou similar", `6` = "Via pública",
                 `7` = "Comércio/serviços", `8` = "Indústria/construção", `9` = "Outro"),
  SIM_NAO    = c(`1` = "Sim", `2` = "Não", `8` = "Não se aplica"),
  AUTOR_SEXO = c(`1` = "Masculino", `2` = "Feminino", `3` = "Ambos os sexos"),
  NUM_ENVOLV = c(`1` = "Um", `2` = "Dois ou mais")
)
IGN <- "Ignorado/em branco"
rotular <- function(x, mapa, ordem = unique(mapa)) {
  y <- unname(mapa[as.character(x)])
  y[is.na(y)] <- IGN   # NA, vazio e códigos 9/99 (ignorado) ou não previstos
  factor(y, levels = c(union(ordem, unique(mapa)), IGN))  # nenhuma categoria vira NA
}
# Códigos fora do dicionário (diferentes de 9/99/NA) são reportados, não descartados em silêncio
checar_codigos <- function(x, mapa, nome) {
  fora <- setdiff(unique(x[!is.na(x)]), c(as.integer(names(mapa)), 9L, 99L))
  if (length(fora)) aviso("%s tem códigos fora do dicionário: %s (tratados como ignorado).",
                          nome, paste(fora, collapse = ", "))
}
checar_codigos(base$CS_RACA, ROT$CS_RACA, "CS_RACA")
checar_codigos(base$CS_ESCOL_N, ROT$CS_ESCOL_N, "CS_ESCOL_N")
checar_codigos(base$LOCAL_OCOR, ROT$LOCAL_OCOR, "LOCAL_OCOR")
checar_codigos(base$SIT_CONJUG, ROT$SIT_CONJUG, "SIT_CONJUG")


# ============================================================
# 6. CONSTRUÇÃO DA BASE ANALÍTICA
# ============================================================

sim <- function(x) as.integer(x %in% 1)   # 1 = Sim; 2, 9 e NA = 0
algum_sim <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(rep(NA_integer_, nrow(df)))
  as.integer(rowSums(sapply(df[vars], function(v) v %in% 1)) > 0)
}
# Indicador de vínculo ignorado: nenhuma das REL_* preenchida com 1 ou 2
vinc_preenchido <- rowSums(sapply(base[VARS_REL], function(v) v %in% c(1, 2))) > 0

base <- base |>
  mutate(
    faixa_etaria = cut(idade_anos, c(-Inf, 9, 19, 29, 39, 49, 59, Inf),
                       labels = c("0–9", "10–19", "20–29", "30–39", "40–49", "50–59", "60+")),
    faixa_etaria = factor(ifelse(is.na(faixa_etaria), IGN, as.character(faixa_etaria)),
                          levels = c(levels(faixa_etaria), IGN)),
    raca_cor     = rotular(CS_RACA, ROT$CS_RACA, c("Branca", "Preta", "Parda", "Amarela", "Indígena")),
    raca_negra   = factor(case_when(raca_cor %in% c("Preta", "Parda") ~ "Negra",
                                    raca_cor == IGN ~ IGN, TRUE ~ as.character(raca_cor)),
                          levels = c("Branca", "Negra", "Amarela", "Indígena", IGN)),
    escolaridade = rotular(CS_ESCOL_N, ROT$CS_ESCOL_N),
    situacao_conjugal = rotular(SIT_CONJUG, ROT$SIT_CONJUG),
    local_ocorrencia  = rotular(LOCAL_OCOR, ROT$LOCAL_OCOR),
    violencia_repetida = rotular(OUT_VEZES, ROT$SIM_NAO, c("Sim", "Não")),
    autor_sexo   = rotular(AUTOR_SEXO, ROT$AUTOR_SEXO),
    autor_alcool = rotular(AUTOR_ALCO, ROT$SIM_NAO, c("Sim", "Não")),
    n_envolvidos = rotular(NUM_ENVOLV, ROT$NUM_ENVOLV),
    deficiencia  = rotular(DEF_TRANS, ROT$SIM_NAO, c("Sim", "Não")),
    # Tipos de violência (múltipla resposta)
    tv_fisica = sim(VIOL_FISIC), tv_psicologica = sim(VIOL_PSICO), tv_tortura = sim(VIOL_TORT),
    tv_sexual = sim(VIOL_SEXU), tv_trafico = sim(VIOL_TRAF), tv_financeira = sim(VIOL_FINAN),
    tv_negligencia = sim(VIOL_NEGLI), tv_trab_infantil = sim(VIOL_INFAN),
    tv_intervencao_legal = sim(VIOL_LEGAL), tv_outra = sim(VIOL_OUTR),
    # Vínculo com a vítima (múltipla resposta)
    vin_parceiro_intimo = algum_sim(base, c("REL_CONJ", "REL_EXCON", "REL_NAMO", "REL_EXNAM")),
    vin_familiar        = algum_sim(base, c("REL_PAI", "REL_MAE", "REL_PAD", "REL_MAD", "REL_FILHO", "REL_IRMAO")),
    vin_conhecido       = algum_sim(base, "REL_CONHEC"),
    vin_desconhecido    = algum_sim(base, "REL_DESCO"),
    vin_outros          = algum_sim(base, c("REL_CUIDA", "REL_PATRAO", "REL_INST", "REL_POL", "REL_PROPRI", "REL_OUTROS")),
    vin_preenchido      = as.integer(vinc_preenchido)
  ) |>
  left_join(UFS[, c("cod_uf", "sigla_uf", "nome_uf", "cod_regiao", "regiao")], by = "cod_uf")
rm(vinc_preenchido)

# Variáveis descritivas (nome na base = rótulo nas tabelas)
VARS_DESC <- c(faixa_etaria = "Faixa etária (anos)", raca_cor = "Raça/cor",
               raca_negra = "Raça/cor (pretas + pardas = negras)", escolaridade = "Escolaridade",
               situacao_conjugal = "Situação conjugal", deficiencia = "Deficiência/transtorno",
               local_ocorrencia = "Local de ocorrência", violencia_repetida = "Violência de repetição",
               autor_sexo = "Sexo do(a) provável autor(a)", autor_alcool = "Suspeita de uso de álcool pelo autor",
               n_envolvidos = "Número de envolvidos")
VARS_TIPO <- c(tv_fisica = "Física", tv_psicologica = "Psicológica/moral", tv_tortura = "Tortura",
               tv_sexual = "Sexual", tv_trafico = "Tráfico de seres humanos", tv_financeira = "Financeira/econômica",
               tv_negligencia = "Negligência/abandono", tv_trab_infantil = "Trabalho infantil",
               tv_intervencao_legal = "Intervenção legal", tv_outra = "Outra")
VARS_VINC <- c(vin_parceiro_intimo = "Parceiro(a) ou ex-parceiro(a) íntimo(a)", vin_familiar = "Familiar",
               vin_conhecido = "Amigo(a)/conhecido(a)", vin_desconhecido = "Desconhecido(a)",
               vin_outros = "Outros vínculos")

saveRDS(base, file.path(DIR_SAIDA, "base_analitica.rds"))
msg("Base analítica: %s notificações; %s sem UF de %s válida.",
    num_br(nrow(base)), num_br(sum(is.na(base$cod_uf))), UF_REFERENCIA)

# Indicadores de múltipla resposta ficam NA nos anos em que a variável de
# origem não existe no arquivo (evita contar "não" onde não houve pergunta).
ORIGEM_IND <- list(tv_fisica = "VIOL_FISIC", tv_psicologica = "VIOL_PSICO", tv_tortura = "VIOL_TORT",
                   tv_sexual = "VIOL_SEXU", tv_trafico = "VIOL_TRAF", tv_financeira = "VIOL_FINAN",
                   tv_negligencia = "VIOL_NEGLI", tv_trab_infantil = "VIOL_INFAN",
                   tv_intervencao_legal = "VIOL_LEGAL", tv_outra = "VIOL_OUTR")
if (!is.null(aud_ausentes) && nrow(aud_ausentes)) {
  for (ind in names(ORIGEM_IND)) {
    anos_sem <- aud_ausentes$ano_arquivo[aud_ausentes$variavel_ausente == ORIGEM_IND[[ind]]]
    if (length(anos_sem)) base[[ind]][base$ano_arquivo %in% anos_sem] <- NA_integer_
  }
}


# ============================================================
# 7. CÁLCULO DAS TAXAS
# ============================================================
# Taxa de notificação de violência contra mulheres (UF u, ano t):
#
#   taxa(u,t) = [ notificações em mulheres residentes em u, notificadas em t,
#                 na faixa IDADE_MIN–IDADE_MAX, exceto autoprovocadas ]
#               / [ população feminina de u em t, mesma faixa etária ] × 100.000
#
# Denominador: Projeções da População do IBGE (revisão 2024), sexo "Mulheres",
# soma das idades simples da faixa. Região e Brasil: soma de numeradores e
# denominadores (não média das taxas). No Brasil o numerador inclui também
# as notificações sem UF válida.

pop_fem <- pop_bruta |>
  filter(sexo == "Mulheres", as.numeric(idade) >= IDADE_MIN, as.numeric(idade) <= IDADE_MAX) |>
  pivot_longer(all_of(col_anos_pop), names_to = "ano", values_to = "pop") |>
  mutate(ano = as.integer(ano), cod_local = as.integer(cod_local)) |>
  filter(ano %in% ANOS_ESTUDO) |>
  group_by(cod_local, sigla, local, ano) |> summarise(pop = sum(pop)) |> ungroup()

# Checagem: soma das UFs = Brasil
chk_pop <- pop_fem |> filter(cod_local %in% UFS$cod_uf) |> group_by(ano) |> summarise(soma_uf = sum(pop)) |>
  left_join(pop_fem |> filter(cod_local == 0) |> select(ano, pop_br = pop), by = "ano") |>
  mutate(dif_rel = abs(soma_uf - pop_br) / pop_br)
if (any(chk_pop$dif_rel > 0.001)) aviso("Soma das UFs difere do total Brasil na projeção em algum ano.")
if (!all(UFS$cod_uf %in% pop_fem$cod_local)) stop("UF ausente na projeção populacional.")

casos_uf <- base |> filter(!is.na(cod_uf)) |> count(cod_uf, ano, name = "casos") |>
  complete(cod_uf = UFS$cod_uf, ano = ANOS_ESTUDO, fill = list(casos = 0L))

taxas_uf <- casos_uf |>
  left_join(UFS, by = "cod_uf") |>
  left_join(pop_fem |> select(cod_local, ano, pop), by = c("cod_uf" = "cod_local", "ano")) |>
  mutate(nivel = "UF", unidade = sigla_uf, nome = nome_uf)

taxas_reg <- taxas_uf |> group_by(cod_regiao, regiao, ano) |>
  summarise(casos = sum(casos)) |> ungroup() |>
  left_join(pop_fem |> select(cod_local, ano, pop), by = c("cod_regiao" = "cod_local", "ano")) |>
  mutate(nivel = "Região", unidade = regiao, nome = regiao)

taxas_br <- base |> count(ano, name = "casos") |>
  complete(ano = ANOS_ESTUDO, fill = list(casos = 0L)) |>
  left_join(pop_fem |> filter(cod_local == 0) |> select(ano, pop), by = "ano") |>
  mutate(nivel = "Brasil", unidade = "BR", nome = "Brasil")

taxas <- bind_rows(taxas_br, taxas_reg, taxas_uf) |>
  mutate(taxa = casos / pop * FATOR_TAXA,
         preliminar = ano %in% ANOS_PRELIMINARES,
         nivel = factor(nivel, c("Brasil", "Região", "UF"))) |>
  select(nivel, unidade, nome, cod_uf, sigla_uf, regiao, ano, casos, pop, taxa, preliminar) |>
  arrange(nivel, unidade, ano)

if (any(is.na(taxas$pop))) stop("Há unidade/ano sem população; verifique a junção.")
uf_ano_zero <- taxas |> filter(nivel == "UF", casos == 0)
if (nrow(uf_ano_zero))
  aviso("UF-ano com zero notificações (taxa = 0; impede log na tendência): %s",
        paste0(uf_ano_zero$unidade, "/", uf_ano_zero$ano, collapse = ", "))

tab_taxas_largo <- taxas |> mutate(taxa = round(taxa, 2)) |>
  select(nivel, unidade, ano, taxa) |>
  pivot_wider(names_from = ano, values_from = taxa)
salvar_csv(taxas, "taxas_unidade_ano")
salvar_csv(tab_taxas_largo, "taxas_largo")


# ============================================================
# 8. ANÁLISE DESCRITIVA E EXPLORATÓRIA
# ============================================================
msg("Análise descritiva...")
anos_lab <- function(a) ifelse(a %in% ANOS_PRELIMINARES, paste0(a, "*"), a)
nota_prelim <- if (length(ANOS_PRELIMINARES)) paste0(" * Dados preliminares (", paste(ANOS_PRELIMINARES, collapse = ", "), ").") else ""

# --- 8.1 Frequências: perfil das vítimas, ocorrência e autor --------------------
freq_tab <- function(df, var, rotulo) {
  df |> count(categoria = .data[[var]], name = "n") |>
    mutate(variavel = rotulo, pct = 100 * n / sum(n),
           pct_validos = ifelse(categoria == IGN, NA, 100 * n / sum(n[categoria != IGN]))) |>
    select(variavel, categoria, n, pct, pct_validos)
}
tab_perfil <- bind_rows(lapply(names(VARS_DESC), function(v) freq_tab(base, v, VARS_DESC[[v]])))

# Tipos de violência e vínculo (múltipla resposta: % não soma 100)
mult_tab <- function(df, vars, grupo) {
  bind_rows(lapply(names(vars), function(v) data.frame(
    grupo = grupo, categoria = vars[[v]], n = sum(df[[v]], na.rm = TRUE),
    n_validos = sum(!is.na(df[[v]])),
    pct = 100 * mean(df[[v]], na.rm = TRUE))))
}
tab_tipos <- mult_tab(base, VARS_TIPO, "Tipo de violência")
tab_vinc  <- bind_rows(mult_tab(base, VARS_VINC, "Vínculo com a vítima"),
                       data.frame(grupo = "Vínculo com a vítima", categoria = "Nenhum vínculo preenchido (ignorado)",
                                  n = sum(base$vin_preenchido == 0), n_validos = nrow(base),
                                  pct = 100 * mean(base$vin_preenchido == 0)))

# Evolução anual das características (proporções por ano, excluindo ignorados)
evol_carac <- bind_rows(lapply(names(VARS_DESC), function(v) {
  base |> filter(.data[[v]] != IGN) |> count(ano, categoria = .data[[v]]) |>
    group_by(ano) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |>
    mutate(variavel = VARS_DESC[[v]])
}))
evol_tipos <- base |> group_by(ano) |>
  summarise(across(all_of(names(VARS_TIPO)), ~ 100 * mean(.x, na.rm = TRUE))) |>
  pivot_longer(-ano, names_to = "var", values_to = "pct") |>
  mutate(tipo = VARS_TIPO[var])

# Completude: % ignorado/em branco por variável e ano
completude <- bind_rows(
  lapply(names(VARS_DESC), function(v) base |> group_by(ano) |>
           summarise(pct_ignorado = 100 * mean(.data[[v]] == IGN)) |> mutate(variavel = VARS_DESC[[v]])),
  base |> group_by(ano) |> summarise(pct_ignorado = 100 * mean(vin_preenchido == 0)) |>
    mutate(variavel = "Vínculo com a vítima"),
  base |> group_by(ano) |> summarise(pct_ignorado = 100 * mean(is.na(cod_uf))) |>
    mutate(variavel = paste("UF de", c(residencia = "residência", notificacao = "notificação", ocorrencia = "ocorrência")[[UF_REFERENCIA]]))
)
tab_completude <- completude |> mutate(pct_ignorado = round(pct_ignorado, 1)) |>
  pivot_wider(names_from = ano, values_from = pct_ignorado)

# Cobertura do sistema: municípios notificantes por UF e ano (município de notificação)
cobertura <- base |> filter(!is.na(ID_MUNICIP), ID_MUNICIP != "") |>
  mutate(uf_not = para_cod(substr(ID_MUNICIP, 1, 2))) |>
  group_by(uf_not, ano) |> summarise(municipios_notificantes = n_distinct(ID_MUNICIP)) |> ungroup() |>
  left_join(UFS[, c("cod_uf", "sigla_uf")], by = c("uf_not" = "cod_uf"))
cobertura_br <- base |> filter(!is.na(ID_MUNICIP), ID_MUNICIP != "") |> group_by(ano) |>
  summarise(municipios_notificantes = n_distinct(ID_MUNICIP))

salvar_csv(tab_perfil, "tab_perfil_vitimas"); salvar_csv(tab_tipos, "tab_tipos_violencia")
salvar_csv(tab_vinc, "tab_vinculo_autor"); salvar_csv(tab_completude, "tab_completude")
salvar_csv(cobertura, "tab_municipios_notificantes")

# --- 8.2 Gráficos ---------------------------------------------------------------
g_casos <- ggplot(taxas_br, aes(factor(ano), casos, fill = ano %in% ANOS_PRELIMINARES)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c(`FALSE` = "#2C5F8A", `TRUE` = "#9DB7CF"), guide = "none") +
  scale_x_discrete(labels = anos_lab) +
  scale_y_continuous(labels = num_br, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Notificações de violência contra mulheres por ano, Brasil",
       x = "Ano de notificação", y = "Número de notificações", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_casos, "01_notificacoes_ano")

g_taxa_reg <- taxas |> filter(nivel != "UF") |>
  ggplot(aes(ano, taxa, colour = unidade, linewidth = unidade == "Brasil")) +
  geom_line() + geom_point(size = 1.3) +
  scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.7), guide = "none") +
  scale_x_continuous(breaks = ANOS_ESTUDO, labels = anos_lab) +
  scale_y_continuous(labels = num_br) +
  scale_colour_manual(values = c(Brasil = "black", Norte = "#1b9e77", Nordeste = "#d95f02",
                                 Sudeste = "#7570b3", Sul = "#e7298a", `Centro-Oeste` = "#66a61e")) +
  labs(title = "Taxa de notificação de violência contra mulheres, Brasil e regiões",
       x = "Ano de notificação", y = "Taxa por 100 mil mulheres", colour = NULL,
       caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_taxa_reg, "02_taxa_brasil_regioes")

g_taxa_uf <- taxas |> filter(nivel == "UF") |>
  mutate(unidade = factor(unidade, UFS$sigla_uf)) |>
  ggplot(aes(ano, taxa)) + geom_line(colour = "#2C5F8A") + geom_point(size = 0.6, colour = "#2C5F8A") +
  facet_wrap(~ unidade, ncol = 6) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025)) +
  scale_y_continuous(labels = num_br) +
  labs(title = "Taxa de notificação de violência contra mulheres por UF",
       subtitle = "Mesma escala em todos os painéis",
       x = "Ano de notificação", y = "Taxa por 100 mil mulheres", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text = element_text(size = 7))
salvar_fig(g_taxa_uf, "03_taxa_por_uf", w = 10, h = 8)

g_box <- taxas |> filter(nivel == "UF") |>
  ggplot(aes(factor(ano), taxa)) + geom_boxplot(fill = "#C9D9E8", outlier.colour = "#B23A48") +
  scale_x_discrete(labels = anos_lab) + scale_y_continuous(labels = num_br) +
  labs(title = "Distribuição das taxas de notificação entre as UFs, por ano",
       x = "Ano de notificação", y = "Taxa por 100 mil mulheres", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_box, "04_distribuicao_taxas_uf")

graf_barras <- function(tab, titulo, nome, col_pct = "pct") {
  g <- ggplot(tab, aes(x = reorder(categoria, .data[[col_pct]]), y = .data[[col_pct]])) +
    geom_col(fill = "#2C5F8A", width = 0.7) + coord_flip() +
    geom_text(aes(label = fmt2(.data[[col_pct]], 1)), hjust = -0.15, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = titulo, x = NULL, y = "Percentual (%)", caption = FONTE_GRAF) + tema_tcc
  salvar_fig(g, nome, w = 7, h = 4.5)
}
for (v in c("faixa_etaria", "raca_negra", "escolaridade", "local_ocorrencia", "situacao_conjugal")) {
  graf_barras(tab_perfil |> filter(variavel == VARS_DESC[[v]]),
              paste0(VARS_DESC[[v]], ", notificações de violência contra mulheres, ",
                     min(ANOS_ESTUDO), "–", max(ANOS_ESTUDO)),
              paste0("05_perfil_", v))
}
graf_barras(tab_tipos, "Tipo de violência (múltipla resposta)", "06_tipos_violencia")
graf_barras(tab_vinc, "Vínculo do(a) provável autor(a) com a vítima (múltipla resposta)", "07_vinculo_autor")

g_evol_tipos <- evol_tipos |> filter(var %in% c("tv_fisica", "tv_psicologica", "tv_sexual", "tv_negligencia", "tv_financeira")) |>
  ggplot(aes(ano, pct, colour = tipo)) + geom_line() + geom_point(size = 1) +
  scale_x_continuous(breaks = ANOS_ESTUDO, labels = anos_lab) +
  labs(title = "Proporção das notificações por tipo de violência, por ano",
       subtitle = "Múltipla resposta: uma notificação pode ter mais de um tipo",
       x = "Ano de notificação", y = "Notificações com o tipo (%)", colour = NULL,
       caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_evol_tipos, "08_evolucao_tipos")

g_compl <- ggplot(completude, aes(factor(ano), variavel, fill = pct_ignorado)) +
  geom_tile(colour = "white") + geom_text(aes(label = round(pct_ignorado)), size = 2.5) +
  scale_fill_gradient(low = "#F1F5F9", high = "#B23A48", name = "% ignorado/em branco") +
  scale_x_discrete(labels = anos_lab) +
  labs(title = "Incompletude das variáveis por ano", x = "Ano de notificação", y = NULL,
       caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_compl, "09_completude", w = 9, h = 5)

g_cob <- ggplot(cobertura_br, aes(ano, municipios_notificantes)) + geom_line(colour = "#2C5F8A") + geom_point() +
  scale_x_continuous(breaks = ANOS_ESTUDO, labels = anos_lab) + scale_y_continuous(labels = num_br) +
  labs(title = "Municípios com ao menos uma notificação de violência contra mulher",
       subtitle = "Indicador da expansão do sistema de vigilância (VIVA/SINAN)",
       x = "Ano de notificação", y = "Número de municípios notificantes", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_cob, "10_municipios_notificantes")


# ============================================================
# 9. PREPARAÇÃO DOS DADOS ESPACIAIS
# ============================================================
msg("Preparando malha e matrizes de vizinhança...")
sf::sf_use_s2(FALSE)

obter_malha <- function() {
  m <- if (!is.null(ARQ_MALHA_LOCAL)) {
    sf::st_read(ARQ_MALHA_LOCAL, quiet = TRUE)
  } else {
    tryCatch(geobr::read_state(code_state = "all", year = 2020, showProgress = FALSE),
             error = function(e) stop("Falha ao baixar a malha pelo geobr (", conditionMessage(e),
                                      "). Baixe a malha de UFs do IBGE e informe ARQ_MALHA_LOCAL."))
  }
  nm <- names(m)
  col_cod <- intersect(c("code_state", "CD_UF", "cod_uf", "COD_UF", "codigo_ibg", "CD_GEOCUF"), nm)
  col_sig <- intersect(c("abbrev_state", "SIGLA_UF", "sigla", "SIGLA", "sigla_uf"), nm)
  if (length(col_cod)) {
    m$cod_uf <- as.integer(as.character(m[[col_cod[1]]]))
  } else if (length(col_sig)) {
    m$cod_uf <- UFS$cod_uf[match(toupper(m[[col_sig[1]]]), UFS$sigla_uf)]
  } else stop("A malha não tem coluna de código ou sigla de UF reconhecível.")
  m <- m |> select(cod_uf) |> filter(cod_uf %in% UFS$cod_uf)
  m <- sf::st_make_valid(sf::st_transform(m, 4674))      # SIRGAS 2000
  m <- m |> group_by(cod_uf) |> summarise() |> ungroup()  # 1 feição por UF
  m <- m |> left_join(UFS, by = "cod_uf") |> arrange(cod_uf)
  falta <- setdiff(UFS$sigla_uf, m$sigla_uf)
  if (length(falta)) stop("UFs ausentes na malha: ", paste(falta, collapse = ", "))
  m
}
malha <- obter_malha()

# --- Matrizes de vizinhança -----------------------------------------------------
# queen: UFs que compartilham fronteira ou vértice (matriz principal; critério
#        natural para áreas contíguas e sem ilhas no nível de UF).
# kNN:   k UFs com centroides mais próximos (k = 3 e k = 2, como no TCC II);
#        garante o mesmo nº de vizinhos para todas, mas é assimétrica.
# Todas padronizadas na linha (style = "W"): cada linha soma 1, e a defasagem
# espacial é a média dos vizinhos.
# Ponto interior de cada UF calculado em projeção plana (Policônica, EPSG:5880)
# e devolvido em graus, para distâncias geodésicas no kNN
coords <- sf::st_geometry(malha) |> sf::st_transform(5880) |> sf::st_point_on_surface() |>
  sf::st_transform(4674) |> sf::st_coordinates()
nb_list <- list(queen = spdep::poly2nb(malha, queen = TRUE, snap = SNAP_CONTIGUIDADE))
for (k in K_VIZINHOS) nb_list[[paste0("knn", k)]] <- spdep::knn2nb(spdep::knearneigh(coords, k = k, longlat = TRUE))
for (k in names(nb_list)) attr(nb_list[[k]], "region.id") <- malha$sigla_uf

diag_viz <- bind_rows(lapply(names(nb_list), function(k) {
  nb <- nb_list[[k]]; card <- spdep::card(nb)
  data.frame(matriz = k, n_regioes = length(nb), n_ligacoes = sum(card),
             media_vizinhos = round(mean(card), 2), min_vizinhos = min(card), max_vizinhos = max(card),
             regioes_sem_vizinhos = paste(malha$sigla_uf[card == 0], collapse = ", "),
             simetrica = spdep::is.symmetric.nb(nb, verbose = FALSE, force = TRUE),
             n_componentes = spdep::n.comp.nb(nb)$nc)
}))
print(diag_viz)
if (any(diag_viz$min_vizinhos == 0)) stop("Há UF sem vizinhos; ajuste SNAP_CONTIGUIDADE ou a matriz.")
if (any(diag_viz$n_componentes > 1))
  aviso("Matriz(es) com mais de um componente conexo: %s", paste(diag_viz$matriz[diag_viz$n_componentes > 1], collapse = ", "))
vizinhos_queen <- data.frame(uf = malha$sigla_uf,
                             vizinhos = sapply(nb_list$queen, function(i) paste(malha$sigla_uf[i], collapse = ", ")))
lw_list <- lapply(nb_list, spdep::nb2listw, style = "W", zero.policy = FALSE)

png(file.path(DIR_FIG, "11_matrizes_vizinhanca.png"), width = 3000, height = 1100, res = 250)
par(mfrow = c(1, length(nb_list)), mar = c(0.5, 0.5, 2, 0.5))
for (k in names(nb_list)) {
  plot(sf::st_geometry(malha), border = "grey60", main = k)
  plot(nb_list[[k]], coords, add = TRUE, col = "#B23A48", points = TRUE, cex = 0.5)
}
dev.off()

# Dados espaciais por ano, na MESMA ordem da malha (essencial para Moran/LISA)
taxa_ano <- function(a) {
  x <- taxas |> filter(nivel == "UF", ano == a) |> arrange(cod_uf)
  stopifnot(identical(as.integer(x$cod_uf), as.integer(malha$cod_uf)))
  if (anyNA(x$taxa)) stop("Taxa ausente em ", a)
  x$taxa
}

# Mapas temáticos das taxas (classes por quintis de todos os UF-anos, para comparar anos)
br_quint <- unique(quantile(taxas$taxa[taxas$nivel == "UF"], probs = seq(0, 1, 0.2)))
mapa_taxas <- malha |> select(cod_uf) |>
  inner_join(taxas |> filter(nivel == "UF") |> select(cod_uf, ano, taxa), by = "cod_uf") |>
  mutate(classe = cut(taxa, br_quint, include.lowest = TRUE, dig.lab = 5),
         ano_lab = anos_lab(ano))
g_mapa_taxas <- ggplot(mapa_taxas) + geom_sf(aes(fill = classe), colour = "white", linewidth = 0.1) +
  facet_wrap(~ ano_lab, ncol = 6) +
  scale_fill_brewer(palette = "Reds", name = "Taxa por 100 mil mulheres (quintis do período)") +
  labs(title = "Taxa de notificação de violência contra mulheres por UF e ano", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_mapa
salvar_fig(g_mapa_taxas, "12_mapa_taxas_por_ano", w = 11, h = 7.5)

# Mapas de características: PROPORÇÃO dentro da UF (não contagem normalizada)
prop_uf <- base |> filter(!is.na(cod_uf)) |> group_by(cod_uf) |>
  summarise(`Vítimas negras (% entre raça/cor informada)` = 100 * mean(raca_negra[raca_negra != IGN] == "Negra"),
            `Ocorrência na residência (% entre local informado)` = 100 * mean(local_ocorrencia[local_ocorrencia != IGN] == "Residência"),
            `Autor parceiro(a)/ex íntimo(a) (% do total)` = 100 * mean(vin_parceiro_intimo),
            `Violência sexual (% do total)` = 100 * mean(tv_sexual, na.rm = TRUE)) |>
  pivot_longer(-cod_uf, names_to = "indicador", values_to = "pct")
salvar_csv(prop_uf |> left_join(UFS[, c("cod_uf", "sigla_uf")], by = "cod_uf"), "tab_proporcoes_uf")
mapas_prop <- lapply(unique(prop_uf$indicador), function(ind) {
  ggplot(malha |> left_join(prop_uf |> filter(indicador == ind), by = "cod_uf")) +
    geom_sf(aes(fill = pct), colour = "white", linewidth = 0.1) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "%") +
    labs(title = ind) + tema_mapa + theme(legend.position = "right", plot.title = element_text(size = 9))
})
png(file.path(DIR_FIG, "13_mapas_caracteristicas_uf.png"), width = 3000, height = 2600, res = 300)
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) {
  gridExtra::grid.arrange(grobs = mapas_prop, ncol = 2,
                          bottom = grid::textGrob(paste(FONTE_GRAF, "Período completo."), x = 0.02, hjust = 0, gp = grid::gpar(fontsize = 7)))
} else {
  print(mapas_prop[[1]])
  aviso("Instale 'gridExtra' para o painel de mapas de características (salvo só o 1º mapa).")
}
dev.off()


# ============================================================
# 10. ÍNDICE GLOBAL DE MORAN
# ============================================================
# I = (n / S0) * sum_i sum_j w_ij (x_i - xbar)(x_j - xbar) / sum_i (x_i - xbar)^2
# H0: I = 0 (ausência de autocorrelação espacial; taxas distribuídas ao acaso)
# H1: I != 0 (two.sided) ou I > 0 (greater), conforme ALTERNATIVA_MORAN
# Inferência: (a) aproximação normal sob randomização; (b) permutação (999).
# Decisão pelo p-valor da permutação.
msg("Moran global...")
interp_moran <- function(I, p) {
  if (is.na(p)) return("Não calculado")
  if (p >= ALPHA) return("Sem evidência de autocorrelação espacial: padrão compatível com aleatoriedade")
  if (I > 0) "Autocorrelação positiva: UFs vizinhas têm taxas semelhantes (agrupamento)"
  else "Autocorrelação negativa: UFs vizinhas têm taxas dissimilares"
}
moran_global <- bind_rows(lapply(names(lw_list), function(k) bind_rows(lapply(ANOS_ESTUDO, function(a) {
  x <- taxa_ano(a)
  if (sd(x) == 0) return(data.frame(matriz = k, ano = a, mensagem = "Taxas sem variabilidade"))
  mt <- spdep::moran.test(x, lw_list[[k]], randomisation = TRUE, alternative = ALTERNATIVA_MORAN)
  set.seed(SEMENTE + a)
  mc <- spdep::moran.mc(x, lw_list[[k]], nsim = N_PERMUTACOES, alternative = ALTERNATIVA_MORAN)
  I <- unname(mt$estimate["Moran I statistic"])
  data.frame(matriz = k, ano = a, preliminar = a %in% ANOS_PRELIMINARES,
             I = I, E_I = unname(mt$estimate["Expectation"]), Var_I = unname(mt$estimate["Variance"]),
             z = unname(mt$statistic), p_normal = mt$p.value, p_permutacao = mc$p.value,
             H0 = "I = 0 (sem autocorrelação espacial)",
             H1 = ifelse(ALTERNATIVA_MORAN == "two.sided", "I != 0", ifelse(ALTERNATIVA_MORAN == "greater", "I > 0", "I < 0")),
             decisao = ifelse(mc$p.value < ALPHA, "Rejeita H0", "Não rejeita H0"),
             interpretacao = interp_moran(I, mc$p.value), mensagem = NA_character_)
}))))
salvar_csv(moran_global, "moran_global")

g_moran <- ggplot(moran_global, aes(ano, I, colour = matriz)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line() + geom_point(aes(shape = p_permutacao < ALPHA), size = 2.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = paste("p <", ALPHA), `FALSE` = paste("p ≥", ALPHA)), name = NULL) +
  scale_x_continuous(breaks = ANOS_ESTUDO, labels = anos_lab) +
  labs(title = "Índice global de Moran das taxas de notificação por UF",
       subtitle = sprintf("Teste por permutação (%d); ponto cheio = significativo", N_PERMUTACOES),
       x = "Ano de notificação", y = "I de Moran", colour = "Matriz", caption = paste0(FONTE_GRAF, nota_prelim)) +
  tema_tcc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
salvar_fig(g_moran, "14_moran_global")


# ============================================================
# 11. DIAGRAMA DE ESPALHAMENTO DE MORAN
# ============================================================
# Eixo x: taxa padronizada z_i; eixo y: média padronizada dos vizinhos (Wz)_i.
# Com W padronizada na linha, a inclinação da reta de Wz em z é o I de Moran.
# Q1 Alto-Alto (z > 0, Wz > 0) e Q2 Baixo-Baixo (z < 0, Wz < 0): associação
#   positiva — a UF se parece com seus vizinhos.
# Q3 Alto-Baixo (z > 0, Wz < 0) e Q4 Baixo-Alto (z < 0, Wz > 0): associação
#   negativa — a UF destoa dos vizinhos (possíveis transições ou outliers).
# Os quadrantes são descritivos: não indicam significância (ver LISA).
lw_p <- lw_list[[MATRIZ_PRINCIPAL]]
quadrante <- function(z, wz) factor(case_when(z >= 0 & wz >= 0 ~ "Alto-Alto", z < 0 & wz < 0 ~ "Baixo-Baixo",
                                              z >= 0 & wz < 0 ~ "Alto-Baixo", TRUE ~ "Baixo-Alto"),
                                    levels = c("Alto-Alto", "Baixo-Baixo", "Alto-Baixo", "Baixo-Alto"))
espalhamento <- bind_rows(lapply(ANOS_ESTUDO, function(a) {
  x <- taxa_ano(a); z <- as.numeric(scale(x)); wz <- spdep::lag.listw(lw_p, z)
  data.frame(ano = a, cod_uf = malha$cod_uf, sigla_uf = malha$sigla_uf, taxa = x, z = z, wz = wz,
             quadrante = quadrante(z, wz))
}))
salvar_csv(espalhamento, "diagrama_espalhamento_moran")
cores_quad <- c("Alto-Alto" = "#B2182B", "Baixo-Baixo" = "#2166AC", "Alto-Baixo" = "#F4A582", "Baixo-Alto" = "#92C5DE",
                "Não significativo" = "grey90")

g_esp <- ggplot(espalhamento, aes(z, wz)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
  geom_point(aes(colour = quadrante), size = 1.4) +
  geom_text(aes(label = sigla_uf), size = 1.9, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~ anos_lab(ano), ncol = 6) + scale_colour_manual(values = cores_quad, name = NULL) +
  labs(title = "Diagrama de espalhamento de Moran — taxas de notificação por UF",
       subtitle = paste0("Matriz ", MATRIZ_PRINCIPAL, "; inclinação da reta = I de Moran"),
       x = "Taxa padronizada (z)", y = "Defasagem espacial da taxa padronizada (Wz)",
       caption = paste0(FONTE_GRAF, nota_prelim)) + tema_tcc
salvar_fig(g_esp, "15_diagrama_espalhamento_moran", w = 11, h = 8)

g_quad_mapa <- ggplot(malha |> select(cod_uf) |> inner_join(espalhamento, by = "cod_uf")) +
  geom_sf(aes(fill = quadrante), colour = "white", linewidth = 0.1) +
  facet_wrap(~ anos_lab(ano), ncol = 6) + scale_fill_manual(values = cores_quad, name = NULL, drop = FALSE) +
  labs(title = "Quadrantes do diagrama de espalhamento de Moran (box map)",
       subtitle = paste0("Matriz ", MATRIZ_PRINCIPAL, "; classificação descritiva, sem teste de significância"),
       caption = paste0(FONTE_GRAF, nota_prelim)) + tema_mapa
salvar_fig(g_quad_mapa, "16_mapa_quadrantes_moran", w = 11, h = 7.5)


# ============================================================
# 12. LISA / MORAN LOCAL
# ============================================================
# I_i = z_i * sum_j w_ij z_j. Significância por permutação condicional
# (N_PERMUTACOES): os valores da UF i ficam fixos e os vizinhos são sorteados.
# Cluster só é declarado quando p < ALPHA; o quadrante vem do diagrama.
# Como há 27 testes por ano, informa-se também o p ajustado por FDR
# (Benjamini-Hochberg) como análise de sensibilidade.
msg("LISA...")
lisa <- bind_rows(lapply(names(lw_list), function(k) bind_rows(lapply(ANOS_ESTUDO, function(a) {
  x <- taxa_ano(a)
  lm_ <- spdep::localmoran_perm(x, lw_list[[k]], nsim = N_PERMUTACOES, alternative = "two.sided",
                                iseed = SEMENTE + a)
  col_p <- grep("Sim", colnames(lm_), value = TRUE)
  col_p <- col_p[!grepl("folded", col_p)][1]
  z <- as.numeric(scale(x)); wz <- spdep::lag.listw(lw_list[[k]], z)
  p <- as.numeric(lm_[, col_p])
  data.frame(matriz = k, ano = a, cod_uf = malha$cod_uf, sigla_uf = malha$sigla_uf, taxa = x,
             Ii = as.numeric(lm_[, "Ii"]), E_Ii = as.numeric(lm_[, "E.Ii"]), Z_Ii = as.numeric(lm_[, "Z.Ii"]),
             p_perm = p, p_fdr = p.adjust(p, "fdr"), quadrante = quadrante(z, wz))
}))))
lisa <- lisa |> mutate(
  cluster = factor(ifelse(p_perm < ALPHA, as.character(quadrante), "Não significativo"), levels = names(cores_quad)),
  cluster_fdr = factor(ifelse(p_fdr < ALPHA, as.character(quadrante), "Não significativo"), levels = names(cores_quad)))
salvar_csv(lisa, "lisa_moran_local")

tab_lisa_resumo <- lisa |> filter(cluster != "Não significativo") |>
  group_by(matriz, ano) |>
  summarise(clusters = paste0(sigla_uf, " (", cluster, ")", collapse = "; "),
            n_sig = n(), n_sig_fdr = sum(cluster_fdr != "Não significativo")) |> ungroup()
salvar_csv(tab_lisa_resumo, "lisa_resumo_significativos")

mapa_lisa <- function(k, col = "cluster", sufixo = "") {
  g <- ggplot(malha |> select(cod_uf) |> inner_join(lisa |> filter(matriz == k), by = "cod_uf")) +
    geom_sf(aes(fill = .data[[col]]), colour = "white", linewidth = 0.1) +
    facet_wrap(~ anos_lab(ano), ncol = 6) + scale_fill_manual(values = cores_quad, name = NULL, drop = FALSE) +
    labs(title = paste0("Mapa LISA (Moran local) — taxas de notificação por UF", sufixo),
         subtitle = sprintf("Matriz %s; permutação condicional (%d); clusters com p < %s", k, N_PERMUTACOES, fmt2(ALPHA)),
         caption = paste0(FONTE_GRAF, nota_prelim)) + tema_mapa
  salvar_fig(g, paste0("17_mapa_lisa_", k, ifelse(col == "cluster_fdr", "_fdr", "")), w = 11, h = 7.5)
}
for (k in names(lw_list)) mapa_lisa(k)
mapa_lisa(MATRIZ_PRINCIPAL, "cluster_fdr", " (p ajustado por FDR)")


# ============================================================
# 13. REGRESSÃO DE PRAIS-WINSTEN
# ============================================================
# Modelo: ln(taxa_t) = b0 + b1 * ano_t + e_t,  e_t = rho * e_{t-1} + u_t
# Estimação por mínimos quadrados generalizados factíveis com transformação
# de Prais-Winsten (preserva a 1ª observação), iterando rho até convergência.
# Verificações antes do ajuste: nº mínimo de anos, série completa, taxa > 0.
msg("Prais-Winsten...")
ajustar_pw <- function(d, unidade, periodo, anos, pandemia = FALSE) {
  base_res <- data.frame(unidade = unidade, periodo = periodo, modelo = ifelse(pandemia, "com indicador 2020–2021", "linear"),
                         n_anos = NA_integer_, b1 = NA_real_, ep_b1 = NA_real_, t = NA_real_, gl = NA_integer_,
                         p_valor = NA_real_, rho = NA_real_, convergiu = NA, dw_original = NA_real_,
                         dw_transformado = NA_real_, p_dw_mqo = NA_real_, p_shapiro = NA_real_,
                         b_pandemia = NA_real_, p_pandemia = NA_real_, mensagem = NA_character_)
  d <- d |> filter(ano %in% anos) |> arrange(ano)
  falha <- function(m) { base_res$n_anos <- nrow(d); base_res$mensagem <- m; base_res }
  if (nrow(d) < MIN_ANOS_TENDENCIA) return(falha(sprintf("Série com %d anos (< %d): não ajustada", nrow(d), MIN_ANOS_TENDENCIA)))
  if (length(setdiff(anos, d$ano))) return(falha(paste("Série incompleta; faltam:", paste(setdiff(anos, d$ano), collapse = ", "))))
  if (anyNA(d$taxa)) return(falha("Taxa ausente em algum ano: não ajustada"))
  if (any(d$taxa <= 0)) return(falha(paste("Taxa zero em", paste(d$ano[d$taxa <= 0], collapse = ", "), "(log indefinido): não ajustada")))
  d$log_taxa <- log(d$taxa)
  d$pandemia <- as.integer(d$ano %in% ANOS_PANDEMIA)
  if (pandemia && length(unique(d$pandemia)) < 2) return(falha("Período sem anos de pandemia"))
  form <- if (pandemia) log_taxa ~ ano + pandemia else log_taxa ~ ano
  m <- NULL
  invisible(capture.output(suppressMessages(
    m <- tryCatch(prais::prais_winsten(form, data = d, index = "ano", max_iter = 100, tol = 1e-6),
                  error = function(e) e))))
  if (inherits(m, "error")) return(falha(paste("Erro no ajuste:", conditionMessage(m))))
  rhos <- as.numeric(m$rho)
  conv <- length(rhos) < 101 && (length(rhos) < 2 || abs(diff(tail(rhos, 2))) < 1e-5)
  s <- summary(m); cf <- s$coefficients
  rho <- tail(rhos, 1)
  base_res$n_anos <- nrow(d); base_res$b1 <- cf["ano", 1]; base_res$ep_b1 <- cf["ano", 2]
  base_res$t <- cf["ano", 3]; base_res$p_valor <- cf["ano", 4]; base_res$gl <- m$df.residual
  base_res$rho <- rho; base_res$convergiu <- conv
  base_res$dw_original <- unname(s$dw["original"]); base_res$dw_transformado <- unname(s$dw["transformed"])
  base_res$p_dw_mqo <- tryCatch(lmtest::dwtest(form, data = d, alternative = "two.sided")$p.value, error = function(e) NA)
  res <- residuals(m)
  base_res$p_shapiro <- if (length(res) >= 3 && sd(res) > 0) shapiro.test(res)$p.value else NA
  if (pandemia) { base_res$b_pandemia <- cf["pandemia", 1]; base_res$p_pandemia <- cf["pandemia", 4] }
  msgs <- c(if (!conv) "rho não convergiu", if (abs(rho) >= 0.99) "rho próximo de 1 (série possivelmente não estacionária)")
  base_res$mensagem <- if (length(msgs)) paste(msgs, collapse = "; ") else NA_character_
  base_res
}

series <- split(taxas, taxas$unidade)
ordem_unid <- c("BR", REGIOES$regiao, UFS$sigla_uf)
pw <- bind_rows(
  lapply(names(PERIODOS_TENDENCIA), function(per) bind_rows(lapply(ordem_unid, function(u)
    ajustar_pw(series[[u]], u, per, PERIODOS_TENDENCIA[[per]])))),
  bind_rows(lapply(ordem_unid, function(u)   # sensibilidade: indicador de pandemia
    ajustar_pw(series[[u]], u, PERIODO_PRINCIPAL, PERIODOS_TENDENCIA[[PERIODO_PRINCIPAL]], pandemia = TRUE)))
)


# ============================================================
# 14. APC E INTERVALOS DE CONFIANÇA
# ============================================================
# Com a resposta em logaritmo natural:
#   APC (%)     = (exp(b1) − 1) × 100
#   IC 95% APC  = (exp(b1 ± t_{gl; 0,975} × EP(b1)) − 1) × 100,  gl = n − p
# (o TCC II usava −1 + 10^b1 no texto e 1,96 no IC; ver cabeçalho).
# Classificação (critério único): Crescente se LI > 0; Decrescente se LS < 0;
# Estável se o IC contém 0. Equivale a p < ALPHA no teste t de b1, pois
# ambos usam a mesma distribuição t.
pw <- pw |> mutate(
  t_crit = qt(1 - ALPHA / 2, gl),
  APC = (exp(b1) - 1) * 100,
  APC_LI = (exp(b1 - t_crit * ep_b1) - 1) * 100,
  APC_LS = (exp(b1 + t_crit * ep_b1) - 1) * 100,
  tendencia = case_when(is.na(APC) ~ "Não estimada", APC_LI > 0 ~ "Crescente",
                        APC_LS < 0 ~ "Decrescente", TRUE ~ "Estável"),
  APC_IC = ifelse(is.na(APC), NA, sprintf("%s (%s; %s)", fmt2(APC), fmt2(APC_LI), fmt2(APC_LS))),
  nivel = case_when(unidade == "BR" ~ "Brasil", unidade %in% REGIOES$regiao ~ "Região", TRUE ~ "UF"),
  efeito_pandemia_pct = ifelse(is.na(b_pandemia), NA, (exp(b_pandemia) - 1) * 100)
)
# Coerência entre IC e p-valor (deve ser 100%)
coer <- pw |> filter(!is.na(p_valor)) |> mutate(ok = (p_valor < ALPHA) == (tendencia != "Estável"))
if (!all(coer$ok)) aviso("Classificação pelo IC e pelo p-valor divergem em %d séries (verificar).", sum(!coer$ok))
nao_ajust <- pw |> filter(tendencia == "Não estimada")
if (nrow(nao_ajust)) aviso("Séries sem tendência estimada: %s",
                           paste0(nao_ajust$unidade, " [", nao_ajust$periodo, ": ", nao_ajust$mensagem, "]", collapse = "; "))
salvar_csv(pw, "prais_winsten_completo")


# ============================================================
# 15. DIAGNÓSTICOS
# ============================================================
# (a) Autocorrelação serial: DW do MQO (antes) e do modelo transformado (depois);
#     H0 do teste DW: rho = 0. DW transformado perto de 2 indica correção.
# (b) Normalidade dos resíduos (Shapiro-Wilk; baixa potência com n pequeno).
# (c) Convergência de rho.
# (d) Valores extremos nas taxas: variação anual > +100% ou < −50% e
#     outliers entre UFs no ano (|z robusto| > 3, mediana/MAD).
diag_pw <- pw |> filter(modelo == "linear") |>
  transmute(unidade, periodo, n_anos, rho = round(rho, 3), convergiu,
            dw_original = round(dw_original, 3), dw_transformado = round(dw_transformado, 3),
            p_dw_mqo = signif(p_dw_mqo, 3),
            autocorrelacao_mqo = ifelse(p_dw_mqo < ALPHA, "Sim (justifica Prais-Winsten)", "Não detectada"),
            autocorrelacao_residual = ifelse(is.na(dw_transformado), NA,
                                             ifelse(dw_transformado < 1.5 | dw_transformado > 2.5,
                                                    "Possível (DW transformado fora de 1,5–2,5): interpretar com cautela",
                                                    "Não indicada")),
            p_shapiro = signif(p_shapiro, 3),
            residuos_normais = ifelse(is.na(p_shapiro), NA, ifelse(p_shapiro < ALPHA, "Rejeita normalidade", "Não rejeita")),
            mensagem)
salvar_csv(diag_pw, "diagnosticos_prais_winsten")
dw_ruim <- diag_pw |> filter(grepl("^Possível", autocorrelacao_residual))
if (nrow(dw_ruim))
  aviso("Possível autocorrelação residual após Prais-Winsten (DW transformado fora de 1,5–2,5) em %s. Uma única reta log-linear pode não descrever bem a série (ex.: expansão rápida do sistema em 2009–2012). Detalhes em diagnosticos_prais_winsten.csv.",
        paste(tapply(dw_ruim$unidade, dw_ruim$periodo, function(u) length(u)) |>
                (\(x) paste0(x, " séries [", names(x), "]"))(), collapse = "; "))

extremos <- taxas |> filter(nivel == "UF") |> group_by(unidade) |> arrange(ano, .by_group = TRUE) |>
  mutate(var_anual_pct = 100 * (taxa / lag(taxa) - 1)) |> group_by(ano) |>
  mutate(z_robusto = (taxa - median(taxa)) / mad(taxa)) |> ungroup() |>
  filter((!is.na(var_anual_pct) & is.finite(var_anual_pct) & (var_anual_pct > 100 | var_anual_pct < -50)) | abs(z_robusto) > 3) |>
  select(unidade, ano, casos, taxa, var_anual_pct, z_robusto)
salvar_csv(extremos, "diagnostico_valores_extremos")
if (nrow(extremos)) aviso("%d UF-ano com valores extremos (ver diagnostico_valores_extremos.csv).", nrow(extremos))

salvar_csv(diag_viz, "diagnostico_matrizes_vizinhanca")
salvar_csv(vizinhos_queen, "vizinhos_queen")


# ============================================================
# 16. TABELAS FINAIS
# ============================================================
# Tabela A: taxas em anos selecionados
anos_sel <- unique(c(min(ANOS_ESTUDO), 2013, 2017, 2019, 2021, 2023, max(ANOS_ESTUDO)))
anos_sel <- anos_sel[anos_sel %in% ANOS_ESTUDO]
tab_final_taxas <- taxas |> filter(ano %in% anos_sel) |>
  mutate(taxa = fmt2(taxa), ano = anos_lab(ano)) |>
  select(nivel, unidade, ano, taxa) |> pivot_wider(names_from = ano, values_from = taxa)

# Tabela B: tendência (período principal e pré-pandemia)
tab_final_tend <- pw |> filter(modelo == "linear") |>
  select(nivel, unidade, periodo, APC_IC, p_valor, tendencia) |>
  mutate(p_valor = ifelse(p_valor < 0.001, "< 0,001", fmt2(p_valor, 3))) |>
  pivot_wider(names_from = periodo, values_from = c(APC_IC, p_valor, tendencia), names_vary = "slowest")

# Tabela C: Moran global (matriz principal + sensibilidade)
tab_final_moran <- moran_global |>
  mutate(cel = sprintf("%s (p = %s)", fmt2(I, 3), ifelse(p_permutacao < 0.001, "< 0,001", fmt2(p_permutacao, 3))),
         ano = anos_lab(ano)) |>
  select(ano, matriz, cel) |> pivot_wider(names_from = matriz, values_from = cel)

# Tabela D: resumo das tendências por período
tab_final_resumo_tend <- pw |> filter(modelo == "linear", nivel == "UF") |>
  count(periodo, tendencia) |> pivot_wider(names_from = tendencia, values_from = n, values_fill = 0)

print(tab_final_resumo_tend)


# ============================================================
# 17. GRÁFICOS E MAPAS FINAIS
# ============================================================
cores_tend <- c(Crescente = "#B2182B", Estável = "grey75", Decrescente = "#2166AC", `Não estimada` = "white")

g_forest <- pw |> filter(modelo == "linear", periodo == PERIODO_PRINCIPAL, !is.na(APC)) |>
  mutate(unidade = factor(unidade, rev(ordem_unid))) |>
  ggplot(aes(APC, unidade, colour = tendencia)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = APC_LI, xmax = APC_LS), height = 0.25) + geom_point(size = 1.8) +
  scale_colour_manual(values = cores_tend, name = "Tendência") +
  labs(title = paste0("Variação percentual anual (APC) da taxa de notificação, ", PERIODO_PRINCIPAL),
       subtitle = "Regressão de Prais-Winsten; IC 95% com quantil t",
       x = "APC (%) e IC 95%", y = NULL, caption = paste0(FONTE_GRAF, nota_prelim)) + tema_tcc
salvar_fig(g_forest, "18_apc_forest_plot", w = 7.5, h = 8.5)

mapa_tend <- malha |> select(sigla_uf) |>
  inner_join(pw |> filter(modelo == "linear", nivel == "UF") |> select(unidade, periodo, APC, tendencia),
             by = c("sigla_uf" = "unidade")) |>
  mutate(tendencia = factor(tendencia, names(cores_tend)))
g_mapa_tend <- ggplot(mapa_tend) + geom_sf(aes(fill = tendencia), colour = "grey40", linewidth = 0.1) +
  geom_sf_text(aes(label = ifelse(is.na(APC), "", fmt2(APC, 1))), size = 1.8) +
  facet_wrap(~ periodo) + scale_fill_manual(values = cores_tend, name = "Tendência", drop = FALSE) +
  labs(title = "Tendência das taxas de notificação por UF (Prais-Winsten)",
       subtitle = "Rótulo = APC (%)", caption = paste0(FONTE_GRAF, nota_prelim)) + tema_mapa
suppressWarnings(salvar_fig(g_mapa_tend, "19_mapa_tendencia_apc", w = 11, h = 5))  # aviso de rótulo em lon/lat


# ============================================================
# 18. EXPORTAÇÃO DOS RESULTADOS
# ============================================================
parametros <- data.frame(
  parametro = c("Período", "Anos preliminares", "Data de download", "Faixa etária", "Exclui autoprovocada",
                "UF de referência", "Remove duplicatas exatas", "Fator da taxa", "Alfa", "Permutações",
                "Alternativa Moran", "Matriz principal", "kNN", "Períodos de tendência", "Anos de pandemia",
                "Semente", "R", "spdep", "prais", "sf"),
  valor = c(paste(range(ANOS_ESTUDO), collapse = "–"), paste(ANOS_PRELIMINARES, collapse = ", "), DATA_DOWNLOAD,
            paste(IDADE_MIN, "a", IDADE_MAX, "anos"), EXCLUIR_AUTOPROVOCADA, UF_REFERENCIA,
            REMOVER_DUPLICATAS_EXATAS, FATOR_TAXA, ALPHA, N_PERMUTACOES, ALTERNATIVA_MORAN, MATRIZ_PRINCIPAL,
            paste(K_VIZINHOS, collapse = ", "), paste(names(PERIODOS_TENDENCIA), collapse = " | "),
            paste(ANOS_PANDEMIA, collapse = ", "), SEMENTE, R.version.string,
            as.character(packageVersion("spdep")), as.character(packageVersion("prais")),
            as.character(packageVersion("sf"))))
tab_avisos <- data.frame(aviso = if (length(AVISOS)) unlist(AVISOS) else "Nenhum aviso")

limitar <- function(x) as.data.frame(x)  # writexl exige data.frame simples
writexl::write_xlsx(list(
  parametros = parametros, avisos = tab_avisos, fluxo_selecao = limitar(tab_fluxo),
  auditoria_arquivos = aud_geral, auditoria_estrutura = limitar(aud_estrutura),
  auditoria_populacao = aud_pop, completude = limitar(tab_completude),
  perfil = tab_perfil, tipos_violencia = tab_tipos, vinculo = tab_vinc,
  municipios_notificantes = limitar(cobertura), taxas = limitar(taxas),
  taxas_anos_selec = limitar(tab_final_taxas), vizinhanca = diag_viz, vizinhos_queen = vizinhos_queen,
  moran_global = moran_global, moran_global_resumo = limitar(tab_final_moran),
  espalhamento = espalhamento, lisa = lisa, lisa_resumo = limitar(tab_lisa_resumo),
  prais_winsten = limitar(pw), tendencia_resumo = limitar(tab_final_tend),
  tendencia_contagem = limitar(tab_final_resumo_tend), diagnosticos_pw = diag_pw,
  valores_extremos = limitar(extremos)
), file.path(DIR_SAIDA, "resultados_tcc_violencia_mulher.xlsx"))

writeLines(c(capture.output(sessionInfo()), "", "AVISOS:", unlist(AVISOS)),
           file.path(DIR_SAIDA, "sessao_e_avisos.txt"))

msg("Concluído. Resultados em: %s", DIR_SAIDA)
if (length(AVISOS)) { message("\nAVISOS PARA REVISÃO:"); for (a in AVISOS) message(" - ", a) }
