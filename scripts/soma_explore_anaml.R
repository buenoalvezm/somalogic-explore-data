data_previous <- read_adat("data/L0125003149_v5.0_EDTAPlasma.bridge.hybNorm.medNormInt.plateScale.leakDetection.calibrate.anmlQC.qcCheck.adat") # ANML normalized file
data_previous_2 <- read_adat("data/L0125003149_v5.0_EDTAPlasma.adat") # ANML normalized file

library(pheatmap)

final_data <- read_tsv(file = "../soma_data/data_phase2_somalogic_curated_20251021.tsv") 
processed_annots_ensembl <- read_tsv(file = "../soma_data/meta_phase2_somalogic_elod_20251021.tsv") 

# [ ] Correlation before and after ANML?
# [ ] Make PCA - associated with scale factors?
# [ ] s age significantly associated with normalisation scale factor (after bonferroni)
# [ ] Correlation scale factors with disease
# [ ] Amount of values > 6 MAD from population mean

cor_previous_2 <- 
  data_previous_2 |>
  as_tibble() |>
  filter(SampleType == "Sample") |> 
  dplyr::select("SampleId", starts_with("seq"))|> 
  column_to_rownames("SampleId") |> 
  cor()

cor_previous <- 
  data_previous |>
  as_tibble() |>
  filter(SampleType == "Sample") |> 
  dplyr::select("SampleId", starts_with("seq"))|> 
  column_to_rownames("SampleId") |> 
  cor()

cor_final <- 
  final_data |>
  column_to_rownames("DAid") |> 
  cor()


cor_previous |> pheatmap()
cor_final
