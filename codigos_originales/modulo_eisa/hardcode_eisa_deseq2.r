## ================================================================
##  EISA — Exon/Intron Split Analysis con DESeq2
##  Modelo con términos de interacción
##  + Exportación de genes descartados por filtro mínimo de conteos
## ================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(tools)
  library(svglite)
})

## =========================
## 1) Paths (I/O)
## =========================

PATH_project <- "/mnt/d/work_dir/memoria/obj2/EISA"

# Tablas de conteo por featureCounts
EXON_COUNTS     <- "/mnt/d/work_dir/memoria/obj2/featureCounts/ver_reversely_stranded/exon_counts.txt"
GENEBODY_COUNTS <- "/mnt/d/work_dir/memoria/obj2/featureCounts/ver_reversely_stranded/genebody_counts.txt"

# gene_id -> gene_name
GENE_ID_NAME <- "/mnt/d/work_dir/memoria/obj2/anotacion/genes_id_name.tsv"

NAME_PLOT <- "EISA (DESeq2 interacción)"

# Outputs
if (!exists("OUT_DIR") || is.null(OUT_DIR) || !nzchar(OUT_DIR)) {
  OUT_DIR <- file.path(PATH_project, "results_deseq2_eisa")
}

TABLES_DIR    <- file.path(OUT_DIR, "tables")
PLOTS_DIR     <- file.path(OUT_DIR, "plots")
PLOTS_VEC_DIR <- file.path(PLOTS_DIR, "vectorial")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2) Parámetros generales
## =========================

KOs <- c("RNU5A", "RNU5B", "RNU5D", "RNU5E", "RNU5F")

# Pre-filtro:
# Se exige >= MIN_COUNT reads en >= MIN_REPS_WITH_COUNT réplicas
# en al menos MIN_CONDITIONS_REQUIRED condiciones.
APPLY_MIN_COUNT_FILTER  <- TRUE
MIN_COUNT               <- 10
MIN_REPS_WITH_COUNT     <- 4
MIN_CONDITIONS_REQUIRED <- 1

# Significancia para clasificación interna del gráfico
SIG_MODE      <- "padj"
PVALUE_CUTOFF <- 0.005
PADJ_CUTOFF   <- 0.05

# Regla para efecto PTc/interacción en gráfico 
EFFECT_MODE  <- "interaction"
PTC_MIN_MAG  <- 0.58
PTC_USE_ABS  <- TRUE
LFCI_MIN_MAG <- 0.58
LFCI_USE_ABS <- TRUE

# Configuración del gráfico EISA
PLOT_LIMITS_MODE  <- "fixed"
PLOT_LIMITS_FIXED <- c(-5, 5)
PLOT_AUTO_PAD     <- 0.5
PLOT_AUTO_MAX_ABS <- 12

LABEL_TOP_UP <- 0
LABEL_TOP_DOWN <- 0
LABEL_ONLY_SIGNIFICANT <- TRUE

## =========================
## 3) I/O Utilidades
## =========================

simplify_bam_name <- function(x) {
  nm <- file_path_sans_ext(basename(x))
  sub("\\.sorted$", "", nm)
}

map_condition <- function(id) {
  if (grepl("^CTRL\\d+$", id)) return("CTRL")
  if (grepl("^A\\d+$", id))    return("RNU5A")
  if (grepl("^B\\d+$", id))    return("RNU5B")
  if (grepl("^D\\d+$", id))    return("RNU5D")
  if (grepl("^E\\d+$", id))    return("RNU5E")
  if (grepl("^F\\d+$", id))    return("RNU5F")
  return(NA_character_)
}

read_featurecounts <- function(path) {
  
  df <- read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = "#"
  )
  
  idx_len <- match("Length", colnames(df))
  
  if (is.na(idx_len) || ncol(df) <= idx_len) {
    stop("Archivo inválido: ", path)
  }
  
  cnt <- as.data.frame(df[, (idx_len + 1):ncol(df), drop = FALSE])
  
  colnames(cnt) <- vapply(
    colnames(cnt),
    simplify_bam_name,
    character(1)
  )
  
  for (j in seq_len(ncol(cnt))) {
    storage.mode(cnt[[j]]) <- "integer"
  }
  
  cnt
}

read_gene_names <- function(path) {
  
  if (!file.exists(path)) {
    return(NULL)
  }
  
  g <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    col.names = c("gene_id", "gene_name"),
    stringsAsFactors = FALSE
  )
  
  g$gene_id <- sub("\\.\\d+$", "", g$gene_id)
  
  g
}

build_intronic_counts <- function(cntGB, cntEx) {
  
  if (!identical(colnames(cntGB), colnames(cntEx))) {
    
    common <- intersect(colnames(cntGB), colnames(cntEx))
    
    if (length(common) == 0) {
      stop("No hay muestras comunes entre GB y Exon.")
    }
    
    cntGB <- cntGB[, common, drop = FALSE]
    cntEx <- cntEx[, common, drop = FALSE]
    
    message("[AVISO] Columnas no idénticas; se intersectaron ", length(common), " muestras.")
  }
  
  common_genes <- intersect(rownames(cntGB), rownames(cntEx))
  
  if (length(common_genes) == 0) {
    stop("No hay genes comunes entre GB y Exon.")
  }
  
  cntGB <- cntGB[common_genes, , drop = FALSE]
  cntEx <- cntEx[common_genes, , drop = FALSE]
  
  intr <- cntGB - cntEx
  
  # Evita conteos negativos derivados de diferencias de asignación.
  intr[intr < 0] <- 0L
  
  intr
}

## =========================
## 4) Filtro mínimo de conteos
## =========================

filter_by_min_counts <- function(cntEx, cntIn, group_vec,
                                 min_count, min_reps, min_conditions) {
  
  conds <- unique(group_vec)
  conds <- conds[!is.na(conds)]
  
  # ------------------------------------------------------------
  # Evalúa si cada gen pasa el filtro por condición.
  # Para cada condición:
  # TRUE si el gen tiene >= min_count reads en >= min_reps réplicas.
  # ------------------------------------------------------------
  condition_pass_matrix <- function(mtx) {
    
    ok_list <- lapply(conds, function(co) {
      
      cols <- names(group_vec)[group_vec == co]
      
      if (length(cols) == 0) {
        return(rep(FALSE, nrow(mtx)))
      }
      
      rowSums(mtx[, cols, drop = FALSE] >= min_count) >= min_reps
    })
    
    ok_by_cond <- do.call(cbind, ok_list)
    colnames(ok_by_cond) <- conds
    rownames(ok_by_cond) <- rownames(mtx)
    
    as.data.frame(ok_by_cond, check.names = FALSE)
  }
  
  pass_ex_by_cond <- condition_pass_matrix(cntEx)
  pass_in_by_cond <- condition_pass_matrix(cntIn)
  
  # Pasa en exones si cumple en al menos min_conditions condiciones.
  keep_ex <- rowSums(pass_ex_by_cond) >= min_conditions
  
  # Pasa en intrones si cumple en al menos min_conditions condiciones.
  keep_in <- rowSums(pass_in_by_cond) >= min_conditions
  
  # Regla final original:
  # debe pasar en exones Y en intrones.
  keep <- keep_ex & keep_in
  
  collapse_passing_conditions <- function(pass_df) {
    
    apply(pass_df, 1, function(z) {
      
      cc <- colnames(pass_df)[as.logical(z)]
      
      if (length(cc) == 0) {
        return("")
      }
      
      paste(cc, collapse = ",")
    })
  }
  
  filter_report <- data.frame(
    gene_id_full = rownames(cntEx),
    gene_id = sub("\\.\\d+$", "", rownames(cntEx)),
    
    pass_exon_min_count = keep_ex,
    pass_intron_min_count = keep_in,
    pass_min_count_filter = keep,
    
    n_conditions_passing_exon = rowSums(pass_ex_by_cond),
    n_conditions_passing_intron = rowSums(pass_in_by_cond),
    
    exon_conditions_passing = collapse_passing_conditions(pass_ex_by_cond),
    intron_conditions_passing = collapse_passing_conditions(pass_in_by_cond),
    
    fail_reason = ifelse(
      keep,
      "pass",
      ifelse(
        !keep_ex & !keep_in,
        "fail_exon_and_intron",
        ifelse(!keep_ex, "fail_exon", "fail_intron")
      )
    ),
    
    stringsAsFactors = FALSE
  )
  
  # Añadir columnas TRUE/FALSE por condición.
  pass_ex_by_cond_out <- pass_ex_by_cond
  pass_in_by_cond_out <- pass_in_by_cond
  
  colnames(pass_ex_by_cond_out) <- paste0("exon_pass_", colnames(pass_ex_by_cond_out))
  colnames(pass_in_by_cond_out) <- paste0("intron_pass_", colnames(pass_in_by_cond_out))
  
  filter_report <- cbind(
    filter_report,
    pass_ex_by_cond_out,
    pass_in_by_cond_out
  )
  
  failed_report <- filter_report[
    !filter_report$pass_min_count_filter,
    ,
    drop = FALSE
  ]
  
  list(
    cntEx = cntEx[keep, , drop = FALSE],
    cntIn = cntIn[keep, , drop = FALSE],
    keep = keep,
    filter_report = filter_report,
    failed_report = failed_report
  )
}

## =========================
## 5) EISA con DESeq2
## =========================

compute_eisa_deseq2 <- function(cntEx, cntIn, ko_label) {
  
  stopifnot(identical(colnames(cntEx), colnames(cntIn)))
  
  # ------------------------------------------------------------
  # 5.1 Combinar conteos exónicos e intrónicos
  # ------------------------------------------------------------
  ex_names <- colnames(cntEx)
  in_names <- colnames(cntIn)
  
  colnames(cntEx) <- paste0("ex_", ex_names)
  colnames(cntIn) <- paste0("in_", in_names)
  
  counts_comb <- cbind(cntEx, cntIn)
  
  fraction <- factor(
    c(rep("ex", length(ex_names)), rep("in", length(in_names))),
    levels = c("in", "ex")
  )
  
  cond_ex <- ifelse(grepl("^CTRL\\d+$", ex_names), "CTRL", "KO")
  
  condition <- factor(
    rep(cond_ex, 2),
    levels = c("CTRL", "KO")
  )
  
  colData <- data.frame(
    row.names = colnames(counts_comb),
    fraction = fraction,
    condition = condition
  )
  
  # ------------------------------------------------------------
  # 5.2 Modelo DESeq2
  # ------------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = counts_comb,
    colData = colData,
    design = ~ fraction + condition + fraction:condition
  )
  
  dds <- DESeq(dds)
  
  # ------------------------------------------------------------
  # 5.3 Identificar coeficientes
  # ------------------------------------------------------------
  rn <- resultsNames(dds)
  
  name_cond <- rn[grepl("^condition_KO_vs_CTRL$", rn)]
  name_int  <- rn[grepl("fraction.*ex.*condition.*KO", rn)]
  
  if (length(name_cond) != 1 || length(name_int) != 1) {
    stop(
      "No se pudieron identificar coeficientes únicos. Nombres: ",
      paste(rn, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------
  # 5.4 Extraer efectos
  # ------------------------------------------------------------
  
  # D_in:
  # KO vs CTRL cuando fraction = in, porque in es el baseline.
  res_in <- as.data.frame(results(dds, name = name_cond))
  
  # D_ex:
  # efecto condition + interacción.
  res_ex <- as.data.frame(results(dds, contrast = list(c(name_cond, name_int))))
  
  # PTc:
  # término de interacción.
  res_int <- as.data.frame(results(dds, name = name_int))
  
  stopifnot(
    identical(rownames(res_in), rownames(res_ex)),
    identical(rownames(res_in), rownames(res_int))
  )
  
  out <- data.frame(
    gene_id_full = rownames(res_in),
    
    D_ex = res_ex$log2FoldChange,
    D_in = res_in$log2FoldChange,
    
    PTc = res_ex$log2FoldChange - res_in$log2FoldChange,
    log2FC_interaction_LFC = res_int$log2FoldChange,
    
    FDR_D_ex = res_ex$padj,
    FDR_D_in = res_in$padj,
    FDR_PTc  = res_int$padj,
    
    stat_Wald = res_int$stat,
    pvalue    = res_int$pvalue,
    padj      = res_int$padj,
    
    stringsAsFactors = FALSE
  )
  
  out
}

## =========================
## 6) Clasificación para gráfico
## =========================

classify_rows <- function(df) {
  
  sig_mask <- switch(
    SIG_MODE,
    
    "pvalue" = !is.na(df$pvalue) & df$pvalue <= PVALUE_CUTOFF,
    
    "padj" = !is.na(df$padj) & df$padj <= PADJ_CUTOFF,
    
    "both" = (!is.na(df$pvalue) & df$pvalue <= PVALUE_CUTOFF) &
             (!is.na(df$padj)   & df$padj   <= PADJ_CUTOFF),
    
    stop("SIG_MODE debe ser 'pvalue', 'padj' o 'both'")
  )
  
  pass_ptc_up <- if (PTC_USE_ABS) {
    abs(df$PTc) >= PTC_MIN_MAG & df$PTc > 0
  } else {
    df$PTc >= PTC_MIN_MAG
  }
  
  pass_ptc_down <- if (PTC_USE_ABS) {
    abs(df$PTc) >= PTC_MIN_MAG & df$PTc < 0
  } else {
    df$PTc <= -PTC_MIN_MAG
  }
  
  li <- df$log2FC_interaction_LFC
  
  pass_li_up <- if (LFCI_USE_ABS) {
    abs(li) >= LFCI_MIN_MAG & li > 0
  } else {
    li >= LFCI_MIN_MAG
  }
  
  pass_li_down <- if (LFCI_USE_ABS) {
    abs(li) >= LFCI_MIN_MAG & li < 0
  } else {
    li <= -LFCI_MIN_MAG
  }
  
  pass_up <- switch(
    EFFECT_MODE,
    "PTc"         = pass_ptc_up,
    "interaction" = pass_li_up,
    "both"        = pass_ptc_up & pass_li_up,
    stop("EFFECT_MODE debe ser 'PTc', 'interaction' o 'both'")
  )
  
  pass_down <- switch(
    EFFECT_MODE,
    "PTc"         = pass_ptc_down,
    "interaction" = pass_li_down,
    "both"        = pass_ptc_down & pass_li_down,
    stop("EFFECT_MODE debe ser 'PTc', 'interaction' o 'both'")
  )
  
  cls <- rep("NS", nrow(df))
  cls[sig_mask & pass_up]   <- "PT_up"
  cls[sig_mask & pass_down] <- "PT_down"
  
  df$reg_class <- cls
  
  df
}

## =========================
## 7) Gráfico EISA
## =========================

plot_eisa <- function(df, ko_label,
                      out_png = NULL,
                      out_svg = NULL,
                      width = 6,
                      height = 6,
                      dpi = 300,
                      limits_mode = PLOT_LIMITS_MODE,
                      limits_fixed = PLOT_LIMITS_FIXED,
                      auto_pad = PLOT_AUTO_PAD,
                      auto_max_abs = PLOT_AUTO_MAX_ABS,
                      n_label_up = LABEL_TOP_UP,
                      n_label_down = LABEL_TOP_DOWN,
                      label_only_sig = LABEL_ONLY_SIGNIFICANT) {
  
  lvls <- c("NS", "PT_down", "PT_up")
  
  cls_counts <- table(factor(df$reg_class, levels = lvls))
  total_n <- nrow(df)
  
  labels_vec <- c(
    NS      = paste0("NS (n=", cls_counts["NS"], ")"),
    PT_down = paste0("PT_down (n=", cls_counts["PT_down"], ")"),
    PT_up   = paste0("PT_up (n=", cls_counts["PT_up"], ")"),
    Total   = paste0("Total (n=", total_n, ")")
  )
  
  if (identical(limits_mode, "auto_square")) {
    
    maxabs <- max(abs(c(df$D_in, df$D_ex)), na.rm = TRUE)
    half <- min(auto_max_abs, ceiling(maxabs + auto_pad))
    limits_eff <- c(-half, half)
    
  } else if (identical(limits_mode, "fixed")) {
    
    limits_eff <- limits_fixed
    
  } else {
    
    stop("limits_mode debe ser 'fixed' o 'auto_square'")
  }
  
  breaks_major <- seq(limits_eff[1], limits_eff[2], by = 2)
  breaks_minor <- seq(limits_eff[1], limits_eff[2], by = 1)
  
  if (identical(limits_mode, "fixed")) {
    
    in_x <- df$D_in >= limits_eff[1] & df$D_in <= limits_eff[2]
    in_y <- df$D_ex >= limits_eff[1] & df$D_ex <= limits_eff[2]
    
    df_plot <- df[in_x & in_y, , drop = FALSE]
    
    df_plot$D_in_plot <- df_plot$D_in
    df_plot$D_ex_plot <- df_plot$D_ex
    
  } else {
    
    df_plot <- within(df, {
      D_in_plot <- pmin(pmax(D_in, limits_eff[1]), limits_eff[2])
      D_ex_plot <- pmin(pmax(D_ex, limits_eff[1]), limits_eff[2])
    })
  }
  
  df_plot$label_text <- ifelse(
    is.null(df_plot$gene_name) |
      is.na(df_plot$gene_name) |
      df_plot$gene_name == "",
    df_plot$gene_id,
    df_plot$gene_name
  )
  
  df_lab <- if (isTRUE(label_only_sig)) {
    subset(df_plot, reg_class %in% c("PT_up", "PT_down"))
  } else {
    df_plot
  }
  
  df_up   <- subset(df_lab, PTc > 0)
  df_down <- subset(df_lab, PTc < 0)
  
  sel_up <- head(
    df_up[order(df_up$PTc, decreasing = TRUE), ],
    n = min(n_label_up, nrow(df_up))
  )
  
  sel_down <- head(
    df_down[order(df_down$PTc, decreasing = FALSE), ],
    n = min(n_label_down, nrow(df_down))
  )
  
  lab_df <- rbind(sel_up, sel_down)
  
  legend_total_df <- data.frame(
    D_in_plot = limits_eff[1],
    D_ex_plot = limits_eff[1],
    reg_class = "Total"
  )
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    
    geom_point(
      data = subset(df_plot, reg_class == "NS"),
      aes(x = D_in_plot, y = D_ex_plot, color = reg_class),
      alpha = 0.5,
      size = 1.0
    ) +
    
    geom_point(
      data = subset(df_plot, reg_class == "PT_down"),
      aes(x = D_in_plot, y = D_ex_plot, color = reg_class),
      alpha = 0.7,
      size = 1.0
    ) +
    
    geom_point(
      data = subset(df_plot, reg_class == "PT_up"),
      aes(x = D_in_plot, y = D_ex_plot, color = reg_class),
      alpha = 0.7,
      size = 1.0
    ) +
    
    geom_point(
      data = legend_total_df,
      aes(x = D_in_plot, y = D_ex_plot, color = reg_class),
      inherit.aes = FALSE,
      alpha = 0
    ) +
    
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(x = D_in_plot, y = D_ex_plot, label = label_text),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.3,
      seed = 42
    ) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    
    labs(
      title = paste0(NAME_PLOT, " - ", ko_label),
      x = expression(Log["2"] * "FC (" * D["in"] * ")"),
      y = expression(Log["2"] * "FC (" * D["ex"] * ")"),
      color = NULL
    ) +
    
    scale_color_manual(
      values = c(
        "PT_up" = "#E41A1C",
        "PT_down" = "#0072B2",
        "NS" = "gray50",
        "Total" = "black"
      ),
      breaks = c("NS", "PT_down", "PT_up", "Total"),
      labels = labels_vec
    ) +
    
    scale_x_continuous(
      limits = limits_eff,
      breaks = breaks_major,
      minor_breaks = breaks_minor
    ) +
    
    scale_y_continuous(
      limits = limits_eff,
      breaks = breaks_major,
      minor_breaks = breaks_minor
    ) +
    
    coord_fixed(ratio = 1) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(linewidth = 0.4, colour = "grey85"),
      panel.grid.minor = element_line(linewidth = 0.25, colour = "grey92"),
      legend.background     = element_rect(fill = "white", color = NA),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, size = 2)
      )
    )
  
  if (!is.null(out_png) || !is.null(out_svg)) {
    dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
    dir.create(PLOTS_VEC_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!is.null(out_png)) {
    ggsave(out_png, p, width = width, height = height, dpi = dpi)
  }
  
  if (!is.null(out_svg)) {
    ggsave(
      out_svg,
      p,
      width = width,
      height = height,
      dpi = dpi,
      device = svglite::svglite
    )
  }
  
  invisible(p)
}

## =========================
## 8) Núcleo del análisis
## =========================

EISA_compare_DESeq2 <- function(ko_label,
                                exon_path = EXON_COUNTS,
                                genebody_path = GENEBODY_COUNTS,
                                gene_name_path = GENE_ID_NAME,
                                apply_min_count_filter = APPLY_MIN_COUNT_FILTER,
                                min_count = MIN_COUNT,
                                min_reps = MIN_REPS_WITH_COUNT,
                                min_conditions = MIN_CONDITIONS_REQUIRED,
                                out_dir = OUT_DIR) {
  
  message("== KO: ", ko_label, " ==")
  
  # ------------------------------------------------------------
  # 8.1 Lectura de archivos de conteo
  # ------------------------------------------------------------
  cntEx <- read_featurecounts(exon_path)
  cntGB <- read_featurecounts(genebody_path)
  
  # ------------------------------------------------------------
  # 8.2 Intersección de genes y muestras
  # ------------------------------------------------------------
  common_genes <- intersect(rownames(cntEx), rownames(cntGB))
  common_samples <- intersect(colnames(cntEx), colnames(cntGB))
  
  if (length(common_genes) == 0 || length(common_samples) == 0) {
    stop("No hay genes o muestras comunes.")
  }
  
  cntEx <- cntEx[common_genes, common_samples, drop = FALSE]
  cntGB <- cntGB[common_genes, common_samples, drop = FALSE]
  
  # ------------------------------------------------------------
  # 8.3 Intrones = genebody - exons
  # ------------------------------------------------------------
  cntIn <- build_intronic_counts(cntGB, cntEx)
  
  # ------------------------------------------------------------
  # 8.4 Mapear muestras a condiciones
  # ------------------------------------------------------------
  sample_ids <- colnames(cntEx)
  conditions <- vapply(sample_ids, map_condition, character(1))
  
  if (any(is.na(conditions))) {
    
    bad <- sample_ids[is.na(conditions)]
    
    stop(
      "No se pudo mapear alguna muestra a condición. Ejemplos: ",
      paste(utils::head(bad, 6), collapse = ", "),
      "\nAsegúrate que los IDs de muestra sean como CTRL1, A1, B1, D1, E1 o F1."
    )
  }
  
  # ------------------------------------------------------------
  # 8.5 Seleccionar CTRL + KO actual
  # ------------------------------------------------------------
  keep_cols <- conditions %in% c("CTRL", ko_label)
  
  if (!any(keep_cols)) {
    stop(
      "No hay columnas para CTRL o ",
      ko_label,
      ". IDs vistos: ",
      paste(unique(conditions), collapse = ", ")
    )
  }
  
  cntEx <- cntEx[, keep_cols, drop = FALSE]
  cntIn <- cntIn[, keep_cols, drop = FALSE]
  conditions <- conditions[keep_cols]
  
  # ------------------------------------------------------------
  # 8.6 Filtro mínimo de conteos
  # ------------------------------------------------------------
  if (isTRUE(apply_min_count_filter)) {
    
    fmin <- filter_by_min_counts(
      cntEx,
      cntIn,
      group_vec = setNames(conditions, colnames(cntEx)),
      min_count = min_count,
      min_reps = min_reps,
      min_conditions = min_conditions
    )
    
    # ----------------------------------------------------------
    # Añadir gene_name a reportes del filtro
    # ----------------------------------------------------------
    gn_for_filter <- read_gene_names(gene_name_path)
    
    filter_report <- fmin$filter_report
    failed_min_counts <- fmin$failed_report
    
    if (!is.null(gn_for_filter)) {
      
      filter_report$gene_name <- gn_for_filter$gene_name[
        match(filter_report$gene_id, gn_for_filter$gene_id)
      ]
      
      failed_min_counts$gene_name <- gn_for_filter$gene_name[
        match(failed_min_counts$gene_id, gn_for_filter$gene_id)
      ]
      
    } else {
      
      filter_report$gene_name <- NA_character_
      failed_min_counts$gene_name <- NA_character_
    }
    
    # Reordenar columnas para facilitar revisión.
    filter_report <- filter_report[, c(
      "gene_id",
      "gene_name",
      setdiff(colnames(filter_report), c("gene_id", "gene_name"))
    )]
    
    failed_min_counts <- failed_min_counts[, c(
      "gene_id",
      "gene_name",
      setdiff(colnames(failed_min_counts), c("gene_id", "gene_name"))
    )]
    
    dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
    
    # Reporte completo: genes que pasan y no pasan.
    out_filter_report <- file.path(
      TABLES_DIR,
      paste0("EISA_DESeq2_", ko_label, "_min_count_filter_report.tsv")
    )
    
    write.table(
      filter_report,
      out_filter_report,
      sep = "\t",
      dec = ",",
      quote = FALSE,
      row.names = FALSE
    )
    
    # Reporte específico: solo genes descartados.
    out_failed_min_counts <- file.path(
      TABLES_DIR,
      paste0("EISA_DESeq2_", ko_label, "_failed_min_counts.tsv")
    )
    
    write.table(
      failed_min_counts,
      out_failed_min_counts,
      sep = "\t",
      dec = ",",
      quote = FALSE,
      row.names = FALSE
    )
    
    message(
      "Filtro mínimos conteos: descartados ",
      nrow(failed_min_counts),
      " genes. Lista guardada en: ",
      out_failed_min_counts
    )
    
    message(
      "Reporte completo del filtro guardado en: ",
      out_filter_report
    )
    
    # Continuar con matrices filtradas.
    cntEx <- fmin$cntEx
    cntIn <- fmin$cntIn
    
    message("Filtro mínimos conteos: retenidos ", nrow(cntEx), " genes.")
    
  } else {
    
    message("Filtro mínimos conteos: [DESACTIVADO]")
  }
  
  if (nrow(cntEx) == 0) {
    stop("Sin genes tras el filtro de mínimos conteos.")
  }
  
  # ------------------------------------------------------------
  # 8.7 Modelo DESeq2 y cálculo EISA
  # ------------------------------------------------------------
  res <- compute_eisa_deseq2(
    cntEx,
    cntIn,
    ko_label = ko_label
  )
  
  # ------------------------------------------------------------
  # 8.8 gene_id sin versión + gene_name
  # ------------------------------------------------------------
  res$gene_id <- sub("\\.\\d+$", "", res$gene_id_full)
  res$gene_id_full <- NULL
  
  gn <- read_gene_names(gene_name_path)
  
  if (!is.null(gn)) {
    res <- merge(
      res,
      gn,
      by = "gene_id",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    res$gene_name <- NA_character_
  }
  
  # ------------------------------------------------------------
  # 8.9 Clasificación para gráfico
  # ------------------------------------------------------------
  res <- classify_rows(res)
  
  # ------------------------------------------------------------
  # 8.10 Tabla final EISA
  # ------------------------------------------------------------
  res_out <- res[, c(
    "gene_id",
    "gene_name",
    "D_ex",
    "D_in",
    "PTc",
    "log2FC_interaction_LFC",
    "FDR_D_ex",
    "FDR_D_in",
    "FDR_PTc"
  )]
  
  dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
  
  out_tab <- file.path(
    TABLES_DIR,
    paste0("EISA_DESeq2_", ko_label, ".tsv")
  )
  
  write.table(
    res_out,
    out_tab,
    sep = "\t",
    dec = ",",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # 8.11 Gráficos
  # ------------------------------------------------------------
  out_png <- file.path(
    PLOTS_DIR,
    paste0("EISA_DESeq2_", ko_label, ".png")
  )
  
  out_svg <- file.path(
    PLOTS_VEC_DIR,
    paste0("EISA_DESeq2_", ko_label, ".svg")
  )
  
  plot_eisa(
    res,
    ko_label,
    out_png = out_png,
    out_svg = out_svg
  )
  
  invisible(res_out)
}

## =========================
## 9) Correr todos los KO
## =========================

run_all <- function(kos = KOs) {
  lapply(kos, function(k) EISA_compare_DESeq2(k))
}

if (sys.nframe() == 0) {
  invisible(run_all())
}