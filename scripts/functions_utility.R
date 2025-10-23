
# Functions for visualization
#library(ggrepel)
#library(tidytext)
library(embed)
library(limma)
#library(ggbeeswarm)
#library(patchwork)
#library(ggsci)
#library(eulerr)
#library(ggplotify)
#library(pheatmap)
#library(ggridges)
#library(ggraph)
#library(tidygraph)
#library(ggupset)
library(tidyverse)
library(tidymodels)
library(patchwork)

import_df <- function(file_path) {

  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)

  df <- switch(tolower(file_extension),
               csv = utils::read.csv(file_path, stringsAsFactors = FALSE),
               tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path, guess_max=10000000),
               parquet = arrow::read_parquet(file_path),
               stop("Unsupported file type: ", file_extension))

  df <- tibble::as_tibble(df)
  return(df)
}

savepath <-
  function(savename) {
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)


    return(savename)

  }

savepath_folder <-
  function(folder, savename) {
    result_folder <- paste0("results/", Sys.Date(), "/",folder)
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)


    return(savename)

  }

savepath_data <-
  function(folder, savename) {
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)


    return(savename)

  }

savepath_results <-
  function(folder, savename) {

    dir.create("results/", showWarnings = FALSE)
    result_folder <- paste0("results/", folder)
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)

    return(savename)

  }


do_pca <- function(data,
                   meta = NULL,
                   variable = NULL,
                   wide = T,
                   impute = T,
                   plots = F) {
  if (wide) {
    data_w <-
      data |>
      rename(Sample = DAid)
  } else {
    data_w <-
      data |>
      select(Sample = DAid, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }

  if (impute) {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_pca(all_predictors())

    pca_prep <- prep(pca_rec)

    tidied_pca <- tidy(pca_prep, 3)

  } else {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_pca(all_predictors())

    pca_prep <- prep(pca_rec)

    tidied_pca <- tidy(pca_prep, 2)
  }
  loadings_data <-
    tidied_pca |>
    rename(Assay = terms,
           Value = value,
           PC = component)

  pca_res <-  juice(pca_prep)

  if (plots) {
    # Loadings plot
    loadings_plot <-
      tidied_pca %>%
      filter(component %in% paste0("PC", 1:4)) %>%
      group_by(component) %>%
      top_n(8, abs(value)) %>%
      ungroup() %>%
      mutate(terms = reorder_within(terms, abs(value), component)) %>%
      ggplot(aes(abs(value), terms, fill = value > 0)) +
      geom_col() +
      facet_wrap( ~ component, scales = "free_y") +
      scale_y_reordered() +
      labs(x = "Absolute value of contribution",
           y = NULL, fill = "Positive?") +
      theme_hpa()

    # PCA plot
    pca_plot <-
      pca_res %>%
      left_join(meta |>
                  rename(Sample = DAid), by = "Sample") %>%
      ggplot(aes(PC1, PC2)) +
      geom_point(aes(color = !!sym(variable)), alpha = 0.7, size = 2) +
      labs(color = NULL) +
      theme_hpa() +
      labs(color = variable)

    return(
      list(
        "pca_res" = pca_res,
        "loadings" = loadings_data,
        "pca_plot" = pca_plot,
        "loadings_plot" = loadings_plot
      )
    )
  } else {
    return(list("pca_res" = pca_res,
                "loadings" = loadings_data))
  }

}


do_umap <- function(data,
                    meta = NULL,
                    variable = NULL,
                    wide = T,
                    impute = T,
                    plots = F,
                    seed = 213,
                    n_neighbors = 15) {
  if (wide) {
    data_w <-
      data |>
      rename(Sample = DAid)
  } else {
    data_w <-
      data |>
      select(Sample = DAid, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }

  if (impute) {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)

    set.seed(seed)
    umap_prep <- prep(umap_rec)

  } else {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)

    set.seed(seed)
    umap_prep <- prep(umap_rec)

  }

  umap_res <-  juice(umap_prep)

  if (plots) {
    # Loadings plot
    umap_plot <-
      umap_res |>
      left_join(meta |>
                  rename(Sample = DAid), by = "Sample") |>
      ggplot(aes(UMAP1, UMAP2, color = !!sym(variable))) +
      geom_point(alpha = 0.7, size = 2) +
      theme_hpa()

    return(list("umap_res" = umap_res,
                "umap_plot" = umap_plot))
  } else {
    return(umap_res)
  }

}

# Themes
theme_hpa <-
  function(angled = F,
           axis_x = T,
           axis_y = T,
           facet_title = T) {
    t <-
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        plot.title = element_text(
          face = "bold",
          size = rel(1),
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "bold",
          hjust = 0.5,
          size = rel(1),
          vjust = 1
        ),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "cm"),
        legend.title = element_text(size = rel(1)),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = "grey90", fill = "grey90"),
        strip.text = element_text(face = "bold")
      )

    if (angled) {
      t <-
        t + theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ))
    }

    if (axis_x == F) {
      t <- t +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.x = element_blank()
        )
    }

    if (axis_y == F) {
      t <- t +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank()
        )
    }
    if (facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }


do_limma_disease <-
  function(data_wide,
           metadata,
           disease,
           controls,
           correct = T,
           cutoff = 0) {

    # Select current disease
    dat <-
      data_wide %>%
      inner_join(metadata %>%
                   select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>%
      rename(Group = Disease) %>%
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control"))

    n_males <-
      metadata |>
      filter(Disease == disease,
             Sex == "M") |>
      nrow()

    n_females <-
      metadata |>
      filter(Disease == disease,
             Sex == "F") |>
      nrow()


    if(correct == T) {
      dat <-
        dat |>
        filter(!is.na(Sex),
               !is.na(Age))
    } else {

      dat <- dat
    }

    # Design a model
    if(correct == T) {

      if(n_males == 0 | n_females == 0) {
        design <- model.matrix(~0 + as.factor(dat$Group) + dat$Age)
        colnames(design) <- c("control", "case", "Age")
      } else if(disease %in% pediatric_diseases & controls == "Healthy") {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex))
        colnames(design) <- c("control", "case",  "Sex")
      } else {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age)
        colnames(design) <- c("control", "case",  "Sex", "Age")
      }

    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }

    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)

    # Fit linear model to each protein assay
    dat_fit <-
      dat %>%
      select(-Sex, -Age, -BMI, -Group)  %>%
      column_to_rownames("DAid") %>%
      t()

    fit <- lmFit(dat_fit, design = design, maxit = 100000) #method = "robust",

    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)

    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit, robust = T)

    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n = nrow(ebays_fit$p.value),
               adjust.method = "fdr",
               confint = TRUE)

    DE_res <-
      DE_results %>%
      as_tibble(rownames = "Assay") %>%
      mutate(Disease = disease,
             sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
                             adj.P.Val < 0.05 & logFC > cutoff ~ "significant up",
                             T ~ "not significant"),
             Control = controls)

    return(DE_res)
  }
## Generate volcano plot from differential expression results
plot_volcano <- function(de_results, cutoff = 0, labels_balanced = F) {

  if (labels_balanced) {
    labels <-
      de_results |>
      filter(sig == "significant up") |>
      top_n(n = 10, wt = -log10(adj.P.Val))  |>
      bind_rows(de_results |>
                  filter(sig == "significant down") |>
                  top_n(n = 10, wt = -log10(adj.P.Val))
      )

  } else {
    labels <-
      de_results |>
      top_n(n = 10, wt = -log10(adj.P.Val))
  }


  volcano_plot <-
    de_results |>
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = Assay)) +
    geom_point(size = 1, alpha = 0.4, show.legend = F) +
    geom_text_repel(data = labels, size = 2, show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = -cutoff, linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = cutoff, linetype = 'dashed', color = "darkgrey") +
    scale_color_manual(values = pal_de) +
    theme_hpa() +
    theme(axis.text = element_text(size = 8),
          legend.position = "top")

  return(volcano_plot)
}

# Themes
theme_hpa <-
  function(angled = F, axis_x = T, axis_x_title = F, axis_y = T, facet_title = T) {
    t <-
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )

    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }

    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    }

    if(axis_x_title == T) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }

    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }

pal_de <-
  c("not significant" = "#D3D3D3",
    "significant up" = "#FF7176",
    "significant down" = "#92C9DA")

plot_somamer_boxplots <-
  function(proteins) {

    prots_mapped <-
      annots |>
      filter(AptName %in% proteins) |>
      select(somamer = AptName, Target) |>
      mutate(Protein = paste(Target, somamer, sep = " - "),
             somamer = factor(somamer, levels = proteins)) |>
      arrange(somamer)


    data_boxplot <-
      test_dat |>
      select(DAid, all_of(proteins)) |>
      pivot_longer(cols = -DAid,
                   names_to = "somamer",
                   values_to = "RFU") |>
      left_join(manifest, by = "DAid") |>
      left_join(prots_mapped, by = "somamer") |>
      mutate(Diagnose = ifelse(Diagnose  %in% c("case", "Healthy"), paste(Diagnose, Cohort, sep = " - "), Diagnose),
             Diagnose = factor(Diagnose, levels = order_disease$Diagnose),
             Protein = factor(Protein, levels = prots_mapped$Protein))

    data_boxplot |>
      ggplot(aes(Diagnose,log10(RFU), color = Cohort, fill = Cohort)) +
      geom_quasirandom(size = 0.5) +
      geom_boxplot(outlier.color = NA, alpha = 0.6, show.legend = F) +
      facet_wrap(~Protein, ncol = 1) +
      scale_color_manual(values = pal_cohhorts) +
      scale_fill_manual(values = pal_cohhorts) +
      theme_hpa(angled = T)

  }

plot_de_auto <-
  function(case,
           data) {

    meta <-
      manifest |>
      filter(DAid %in% data$DAid) |>
      select(DAid, Disease, Age, Sex, BMI) |>
      filter(!is.na(Disease))

    test_res <- do_limma_disease(data_wide = data |> mutate(across(-DAid, log10)),
                                 metadata = meta,
                                 disease = case,
                                 controls = "All",
                                 correct = F,
                                 cutoff = 0)

    test_res_mapped <-
      test_res |>
      left_join(annots |>
                  select(AptName, Target),
                by = c("Assay" = "AptName")) |>
      rename(Aptamer = Assay,
             Assay = Target)

    top2 <-
      test_res_mapped |>
      head(2) |>
      pull(Aptamer)

    p1 <- plot_volcano(test_res_mapped) + ggtitle(case)
    p2 <- plot_somamer_boxplots(top2)


    combined_plot <- p1 + p2  + plot_layout(widths = c(1,2))

    return(combined_plot)
  }

do_mixed_effect_model <- function(df,
                                  protein,
                                  type) {

  df <-
    df |>
    filter(Assay == protein) |>
    select(DAid, Age, BMI, womanID, trimester, RFU) |>
    mutate(trimester = as.factor(trimester),
           Age = as.numeric(Age),
           BMI = as.numeric(BMI),
           womanID = as.factor(womanID))


  # Fit the mixed-effects model (no random slopes, random intercepts)
  mixed_model <- lmer(RFU ~ trimester + Age + BMI + (1 | womanID), data = df)

  # Compute R² for fixed and random effects
  r2_values <- r2_nakagawa(mixed_model)

  # Check for singular fit
  if (isSingular(mixed_model, tol = 1e-4)) {
    message(paste("Protein", protein, ": Singular fit detected. Setting random variance explained to 0."))
    marginal_r2 <- r2_values$R2_marginal
    conditional_r2 <- marginal_r2  # No variance from random effects
    random_r2 <- 0  # Explicitly set to 0
  } else {
    marginal_r2 <- r2_values$R2_marginal
    conditional_r2 <- r2_values$R2_conditional
    random_r2 <- conditional_r2 - marginal_r2
  }

  # Store variance explained
  variance_explained <- tibble(
    Component = c("Fixed effects (age & sex)", "Random effects (subject)", "Residual"),
    Variance = c(marginal_r2, random_r2, 1 - conditional_r2)
  ) |> mutate(Protein = protein) |> relocate(Protein)


  # Tidy fixed effects
  fixed_effects <-
    broom.mixed::tidy(mixed_model, effects = "fixed") |>
    mutate(Assay = protein)

  if(type == "fixed_effects") {
    return(fixed_effects)
  } else if(type == "variance_explained") {
    return(variance_explained)
  }

}
