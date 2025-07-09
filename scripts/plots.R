# ===============================
# Análise dos resultados do projeto de Metabarcoding
# Script reestruturado, comentado e funcional
# ===============================

# --- 1. Bibliotecas necessárias ---
packages <- c(
  "dplyr", "ggplot2", "readr", "lubridate", "forcats", "tidyr",
  "stringr", "ggthemes", "tibble", "broom", "GGally", "purrr", "scales"
)
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# --- 2. Definir diretórios base ---
base_dir <- "/home/scmdcunha/estagio"
paths <- list(
  blast_333 = file.path(base_dir, "blast_333"),
  blast_532 = file.path(base_dir, "blast_532"),
  vsearch_333 = file.path(base_dir, "vsearch_333"),
  vsearch_532 = file.path(base_dir, "vsearch_532")
)

# --- 3. Função para carregar dados principais ---
load_data <- function(path) {
  list(
    scores = read.delim(file.path(path, "scores", "final_scores.csv"), stringsAsFactors = FALSE),
    taxonomy = read.delim(file.path(path, "taxonomy", "taxonomy.csv"), sep = ";", stringsAsFactors = FALSE),
    occurrences = read_csv(file.path(path, "occurrences", "occurrences.csv"), show_col_types = FALSE)
  )
}

# --- 4. Carregar dados para cada cenário ---
data_blast_333 <- load_data(paths$blast_333)
data_blast_532 <- load_data(paths$blast_532)
data_vsearch_333 <- load_data(paths$vsearch_333)
data_vsearch_532 <- load_data(paths$vsearch_532)

# --- 5. Preparar e limpar dados principais (exemplo com blast_333) ---
main_scores <- data_blast_333$scores %>%
  mutate(
    EventDate = as.Date(EventDate),
    event_year = year(EventDate),
    Percent.Identity = as.numeric(Percent.Identity),
    Distance_km = as.numeric(Distance_km),
    YearsSinceEvent = as.numeric(difftime(Sys.Date(), EventDate, units = "days")) / 365.25
  ) %>%
  filter(!is.na(Percent.Identity), !is.na(Distance_km), !is.na(YearsSinceEvent))

occ <- data_blast_333$occurrences %>%
  mutate(eventDate = as.Date(eventDate), event_year = year(eventDate))

# --- 6. Histogramas das variáveis principais ---
hist_data <- main_scores %>%
  select(`Percentagem de Identidade` = Percent.Identity,
         `Distância geográfica (km)` = Distance_km,
         `Anos desde a ocorrência` = YearsSinceEvent) %>%
  pivot_longer(cols = everything(), names_to = "Variável", values_to = "Valor")

ggplot(hist_data, aes(x = Valor)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  facet_wrap(~Variável, scales = "free", ncol = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 8)) +
  labs(
    title = "Distribuição das Variáveis Utilizadas na Validação",
    x = NULL,
    y = "Número de identificações"
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

ggsave("fig2_histogramas.png", width = 7, height = 9)

# --- 7. Número de ocorrências por ano ---
occ_year <- occ %>% count(event_year) %>% filter(!is.na(event_year))

ggplot(occ_year, aes(x = event_year, y = n)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue") +
  labs(
    title = "Número de Ocorrências ao Longo dos Anos",
    x = "Ano da Ocorrência (GBIF)",
    y = "Número Total de Ocorrências"
  ) +
  theme_minimal()

ggsave("fig_ocorrencias_ano.png")

# --- 8. Ocorrências por país ---
occ %>%
  count(country, sort = TRUE) %>%
  filter(!is.na(country)) %>%
  ggplot(aes(x = fct_reorder(country, n), y = n)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Distribuição Geográfica das Ocorrências por País",
    x = "País",
    y = "Número de Registos"
  ) +
  theme_minimal()

ggsave("fig_ocorrencias_pais.png")

# --- 9. Top 15 localidades com mais ocorrências ---
occ %>%
  filter(!is.na(locality), locality != "") %>%
  count(locality, sort = TRUE) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = fct_reorder(locality, n), y = n)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  labs(
    title = "Top 15 Localidades com Mais Ocorrências (GBIF)",
    x = "Localidade",
    y = "Número de Registos"
  ) +
  theme_minimal()

ggsave("fig_ocorrencias_localidade.png")

# --- 10. Distribuição das distâncias geográficas ---
ggplot(main_scores, aes(x = Distance_km)) +
  geom_histogram(bins = 30, fill = "purple", color = "black") +
  labs(
    title = "Distribuição das Distâncias Geográficas das Ocorrências",
    x = "Distância (km) ao ponto de referência",
    y = "Número de Ocorrências"
  ) +
  theme_minimal()

ggsave("fig_distribuicao_distancias.png")

# --- 11. Distribuição das ordens taxonómicas (top 15 + Outros) ---
ordem_counts <- main_scores %>%
  mutate(Order_grouped = fct_lump(Order, n = 15, other_level = "Outros")) %>%
  count(Order_grouped, sort = TRUE)

ggplot(ordem_counts, aes(x = n, y = fct_reorder(Order_grouped, n), fill = Order_grouped)) +
  geom_col() +
  labs(
    title = "Distribuição das Ordens Taxonómicas Identificadas",
    x = "Número de Ocorrências",
    y = "Ordem Taxonómica"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("fig_ordens_taxonomicas.png")

# --- 12. Comparação do número de hits validados: BLAST vs VSEARCH (pesos 3:3:3) ---
count_hits <- function(df, label) {
  df %>%
    distinct(Query.ID, Accession.ID) %>%
    summarise(n_hits = n()) %>%
    mutate(tool = label)
}

blast_hits333 <- count_hits(data_blast_333$scores, "BLAST")
vsearch_hits333 <- count_hits(data_vsearch_333$scores, "VSEARCH")

hits_compare <- bind_rows(blast_hits333, vsearch_hits333)

ggplot(hits_compare, aes(x = tool, y = n_hits, fill = tool)) +
  geom_col() +
  labs(
    title = "Número de hits finais validados: BLAST vs VSEARCH (pesos iguais)",
    x = NULL,
    y = "Número de hits"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("fig_hits_comparacao.png")

# --- 13. Comparação dos scores finais entre pesos 3:3:3 e 5:3:2 (BLAST) ---
blast_scores_333 <- data_blast_333$scores %>%
  select(Query.ID, Accession.ID, Score) %>%
  rename(Score_333 = Score)

blast_scores_532 <- data_blast_532$scores %>%
  select(Query.ID, Accession.ID, Score) %>%
  rename(Score_532 = Score)

scores_join <- left_join(blast_scores_333, blast_scores_532, by = c("Query.ID", "Accession.ID")) %>%
  filter(!is.na(Score_532)) %>%
  mutate(Diff = Score_532 - Score_333)

ggplot(scores_join, aes(x = Score_333, y = Score_532)) +
  geom_point(alpha = 0.5, color = "#2c7fb8") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#de2d26", size = 1) +
  labs(
    title = "Comparação dos Scores entre Pesos (3:3:3) e (5:3:2)",
    x = "Score pesos 3:3:3",
    y = "Score pesos 5:3:2"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

summary(scores_join$Diff)

ggsave("fig_comparacao_scores.png")

# --- 14. Média do score por ordem taxonómica (BLAST 3:3:3) ---
mean_scores <- main_scores %>%
  group_by(Order) %>%
  summarise(Media_Score = mean(Score, na.rm = TRUE))

ggplot(mean_scores, aes(x = fct_reorder(Order, Media_Score), y = Media_Score)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Média do Score Final por Ordem Taxonómica",
    x = "Ordem",
    y = "Score médio"
  ) +
  theme_minimal()

ggsave("fig_media_score_ordem.png")

# --- 15. Modelos GLM para Percent Identity vs Distância e Tempo ---

glm_dist <- glm(Percent.Identity ~ Distance_km, data = main_scores, family = Gamma(link = "inverse"))
glm_time <- glm(Percent.Identity ~ YearsSinceEvent, data = main_scores, family = Gamma(link = "inverse"))

seq_dist <- seq(min(main_scores$Distance_km, na.rm = TRUE), max(main_scores$Distance_km, na.rm = TRUE), length.out = 100)
seq_time <- seq(min(main_scores$YearsSinceEvent, na.rm = TRUE), max(main_scores$YearsSinceEvent, na.rm = TRUE), length.out = 100)

pred_dist <- predict(glm_dist, newdata = data.frame(Distance_km = seq_dist), type = "response")
pred_time <- predict(glm_time, newdata = data.frame(YearsSinceEvent = seq_time), type = "response")

pred_df_dist <- data.frame(Distance_km = seq_dist, Percent.Identity = pred_dist)
pred_df_time <- data.frame(YearsSinceEvent = seq_time, Percent.Identity = pred_time)

p_dist <- ggplot(main_scores, aes(x = Distance_km, y = Percent.Identity)) +
  geom_point(alpha = 0.4, color = "#1f78b4") +
  geom_line(data = pred_df_dist, aes(x = Distance_km, y = Percent.Identity), color = "#e31a1c", size = 1.5) +
  labs(title = "Modelo GLM: Percentagem de Identidade vs Distância",
       x = "Distância geográfica (km)",
       y = "Percentagem de Identidade") +
  scale_x_continuous(breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(85, 100, 1)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

p_time <- ggplot(main_scores, aes(x = YearsSinceEvent, y = Percent.Identity)) +
  geom_point(alpha = 0.4, color = "#33a02c") +
  geom_line(data = pred_df_time, aes(x = YearsSinceEvent, y = Percent.Identity), color = "#e31a1c", size = 1.5) +
  labs(title = "Modelo GLM: Percentagem de Identidade vs Tempo",
       x = "Anos desde a ocorrência",
       y = "Percentagem de Identidade") +
  scale_x_continuous(breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(85, 100, 1)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("fig3a_percent_identity_vs_distance.png", plot = p_dist, width = 8, height = 6, dpi = 300)
ggsave("fig3b_percent_identity_vs_time.png", plot = p_time, width = 8, height = 6, dpi = 300)

print(p_dist)
print(p_time)

# --- 16. Comparação dos scores médios por ordem: pesos 3:3:3 vs 5:3:2 ---

mean_scores_333 <- data_blast_333$scores %>%
  group_by(Order) %>%
  summarise(Media_Score = mean(Score, na.rm = TRUE)) %>%
  mutate(Parametros = "Pesos 3:3:3")

mean_scores_532 <- data_blast_532$scores %>%
  group_by(Order) %>%
  summarise(Media_Score = mean(Score, na.rm = TRUE)) %>%
  mutate(Parametros = "Pesos 5:3:2")

ordem_333 <- unique(mean_scores_333$Order)
ordem_532 <- unique(mean_scores_532$Order)
ordens_comuns <- intersect(ordem_333, ordem_532)

mean_scores_333_filt <- mean_scores_333 %>% filter(Order %in% ordens_comuns)
mean_scores_532_filt <- mean_scores_532 %>% filter(Order %in% ordens_comuns)

mean_scores_all <- bind_rows(mean_scores_333_filt, mean_scores_532_filt)

media_combinada <- mean_scores_all %>%
  group_by(Order) %>%
  summarise(media = mean(Media_Score, na.rm = TRUE)) %>%
  arrange(media) %>%
  pull(Order)

mean_scores_all$Order <- factor(mean_scores_all$Order, levels = media_combinada)

p_fig4 <- ggplot(mean_scores_all, aes(x = Order, y = Media_Score, fill = Parametros)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Pesos 3:3:3" = "#1f78b4", "Pesos 5:3:2" = "#e31a1c")) +
  labs(
    title = "Figura 4 – Efeito da parametrização na média do score por ordem taxonómica (ordens comuns)",
    x = "Ordem Taxonómica",
    y = "Score médio",
    fill = "Parametrização"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("fig4_parametrizacao_ordem_comum.png", plot = p_fig4, width = 9, height = 7, dpi = 300)
print(p_fig4)
