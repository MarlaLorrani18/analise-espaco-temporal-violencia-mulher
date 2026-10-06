# Tendência e padrão espacial das notificações de violência contra a mulher no Brasil, 2009–2025

Código em R do Trabalho de Conclusão de Curso do MBA em Data Science & Analytics (USP/Esalq).

**Autora:** Marla Lorrani Santos de Souza
**Orientador:** Daniel Alvarez Firmino

## O que o código faz

A partir das notificações de violência do Sistema de Informação de Agravos de Notificação (SINAN) e das projeções populacionais do IBGE, o código:

1. audita os 17 arquivos anuais do SINAN (estrutura, ausentes, categorias, duplicatas);
2. seleciona as notificações de mulheres, exceto lesões autoprovocadas;
3. calcula as taxas de notificação por 100 mil mulheres, por Unidade da Federação (UF), região e Brasil;
4. faz a análise descritiva (perfil das vítimas, da ocorrência e do autor; completude; municípios notificantes);
5. avalia a autocorrelação espacial (índice global de Moran, diagrama de espalhamento e LISA) com matriz de contiguidade e de k vizinhos mais próximos;
6. estima a tendência temporal pela regressão de Prais-Winsten e a variação percentual anual (APC), com diagnósticos.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `tcc_viol_mulher_analise.R` | Análise completa (seções 1 a 18) |
| `figuras_tcc.R` | Figuras 1 a 4 do texto do TCC (rodar depois do script principal) |

## Dados (não incluídos no repositório)

Os dados são públicos e devem ser baixados e colocados nas subpastas abaixo, dentro da pasta do projeto:

```
pasta-do-projeto/
├── bases_viol_2009-2025/        VIOLBR09.dbc ... VIOLBR25.dbc
└── bases_estimativa_pop_fem/    projecoes_2024_tab1_idade_simples.xlsx
```

- **SINAN, violência interpessoal e autoprovocada:** arquivos `VIOLBRaa.dbc` do DATASUS (Transferência de Arquivos, fonte SINAN, tipo Violência). O ano de 2025 é preliminar.
- **População:** IBGE, Projeções da População, revisão 2024, tabela 1 por idade simples.

## Como rodar

1. Instale o R (versão usada: 4.5.2). Os pacotes que faltarem são instalados pelo próprio script.
2. Abra o R ou RStudio com a pasta do projeto como pasta de trabalho (ou use `setwd()`).
3. Rode `tcc_viol_mulher_analise.R`. Os parâmetros do estudo estão na seção 1.
4. Rode `figuras_tcc.R` para gerar as figuras do texto.

Os resultados ficam em `resultados_R/` (tabelas, figuras, auditoria e uma planilha consolidada).

## Versões usadas

R 4.5.2; read.dbc 1.2.0; geobr 2.1.0; sf 1.1-3; spdep 1.4-2; prais 1.2.0; ggplot2 4.0.0.
