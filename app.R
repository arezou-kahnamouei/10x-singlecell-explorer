# app.R — 10x Single-cell Explorer (Seurat + Shiny) + Advanced Violin & Stats
# --------------------------------------------------------------------------

# ===== Packages =====
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(DT)
  library(Seurat); library(dplyr); library(ggplot2); library(Matrix)
  library(openxlsx); library(tidyr); library(reshape2); library(pheatmap)
  library(readxl); library(stringr); library(ggrepel); library(grid)
  suppressWarnings(suppressMessages({
    if (requireNamespace("HGNChelper",    quietly = TRUE)) library(HGNChelper)
    if (requireNamespace("AnnotationDbi", quietly = TRUE)) library(AnnotationDbi)
  }))
})

options(shiny.maxRequestSize = 1024 * 1024^2) # ~1GB

# ===== Helpers =====
get_assay_layer <- function(x, candidates = c("Gene Expression", "RNA")) {
  if (is.list(x)) {
    nm <- intersect(candidates, names(x))
    stopifnot("No valid assay layer found in 10x list!" = length(nm) > 0)
    x[[nm[1]]]
  } else x
}

.get_orgdb_safe <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  tryCatch(getFromNamespace(pkg, pkg), error = function(e) NULL)
}
.alias_select_safe <- function(species = c("Human","Mouse"), keys) {
  species <- match.arg(species)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(NULL)
  pkg <- if (species == "Human") "org.Hs.eg.db" else "org.Mm.eg.db"
  db <- .get_orgdb_safe(pkg); if (is.null(db)) return(NULL)
  tryCatch(AnnotationDbi::select(db, keys = unique(keys), keytype = "ALIAS", columns = "SYMBOL"),
           error = function(e) NULL)
}

clean_symbols <- function(x, species = c("Human","Mouse","Other/Unknown")) {
  species <- match.arg(species)
  x <- unique(na.omit(stringr::str_trim(as.character(x)))); x <- x[nzchar(x)]
  if (!length(x)) return(character(0))
  x2 <- x
  if (species == "Human" && requireNamespace("HGNChelper", quietly = TRUE)) {
    cl <- HGNChelper::checkGeneSymbols(x2)
    x2 <- ifelse(!is.na(cl$Suggested.Symbol), cl$Suggested.Symbol, cl$x)
  }
  if (species %in% c("Human","Mouse")) {
    amap <- .alias_select_safe(species, keys = x2)
    if (!is.null(amap) && nrow(amap)) {
      keep <- !is.na(amap$SYMBOL); amap <- amap[keep, , drop = FALSE]
      idx  <- match(x2, amap$ALIAS); repl <- amap$SYMBOL[idx]
      x2[!is.na(repl)] <- repl[!is.na(repl)]
    }
  }
  unique(x2)
}

read_genes_flex <- function(path, universe, species = c("Human","Mouse","Other/Unknown")) {
  species <- match.arg(species); ext <- tolower(tools::file_ext(path))
  best_col <- function(df) {
    if (is.null(df) || !ncol(df)) return(character(0))
    best <- character(0); best_match <- -1
    for (j in seq_len(ncol(df))) {
      vals <- clean_symbols(df[[j]], species = species)
      m <- sum(toupper(vals) %in% toupper(universe))
      if (m > best_match) { best <- vals; best_match <- m }
    }
    unique(best)
  }
  out <- character(0)
  if (ext == "csv") {
    df <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    out <- best_col(df)
  } else if (ext %in% c("xls","xlsx")) {
    shs <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    tmp <- lapply(shs, function(s){
      d <- tryCatch(as.data.frame(readxl::read_excel(path, sheet = s), stringsAsFactors = FALSE),
                    error = function(e) NULL)
      best_col(d)
    })
    out <- unique(unlist(tmp, use.names = FALSE))
  }
  unique(out)
}

match_genes <- function(gs, universe){
  u <- universe[match(toupper(gs), toupper(universe))]
  unique(u[!is.na(u)])
}
parse_gene_text <- function(txt, species = c("Human","Mouse","Other/Unknown")) {
  species <- match.arg(species)
  if (is.null(txt) || !nzchar(txt)) return(character(0))
  x <- unlist(strsplit(txt, "[,;\\n\\r\\t ]+"))
  clean_symbols(x, species = species)
}

as_tbl <- function(x) datatable(x, extensions = "Buttons",
                                options = list(dom = "Bfrtip", buttons = c("copy","csv","excel"), pageLength = 10))
safe_id   <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
norm_path <- function(p) { if (!nzchar(p)) return(p); p <- gsub("\\\\","/",p); 
tryCatch(normalizePath(p, winslash="/", mustWork=FALSE), error=function(e) p) }

# ---------- Advanced gene-level helpers ----------
fetch_expr_df <- function(so, gene, zscore = FALSE) {
  md <- so@meta.data
  if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
  md$cluster <- factor(md$seurat_clusters)
  expr <- FetchData(so, vars = gene)[,1]
  if (isTRUE(zscore)) expr <- as.numeric(scale(expr))
  cbind(md[, c("cluster","CMO_label")], expr = expr)
}

violin_stats_tbl <- function(so, gene) {
  df <- fetch_expr_df(so, gene, zscore = FALSE)
  tab <- df %>%
    group_by(cluster, CMO_label) %>%
    summarize(n = n(),
              det_pct = 100 * mean(expr > 0, na.rm = TRUE),
              median_pos = suppressWarnings(median(expr[expr > 0], na.rm = TRUE)),
              .groups = "drop") %>%
    tidyr::complete(cluster, CMO_label, fill = list(n=0, det_pct=0, median_pos=NA_real_)) %>%
    arrange(cluster, CMO_label)
  wide <- tidyr::pivot_wider(tab, names_from = CMO_label, values_from = c(n, det_pct, median_pos))
  wide$delta_median <- wide$median_pos_Treated - wide$median_pos_Control
  wide$who_higher <- dplyr::case_when(
    is.na(wide$delta_median)            ~ NA_character_,
    wide$delta_median >  0              ~ "Treated↑",
    wide$delta_median <  0              ~ "Control↑",
    TRUE                                ~ "≈ equal"
  )
  wide
}

one_gene_de_per_cluster <- function(so, gene, min_cells = 10) {
  md <- so@meta.data
  if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
  md$cluster <- factor(md$seurat_clusters)
  out <- list()
  for (cl in levels(md$cluster)) {
    cells <- rownames(md)[md$cluster == cl]
    sub   <- subset(so, cells = cells); Idents(sub) <- sub$CMO_label
    tb <- table(sub$CMO_label)
    if (!all(c("Control","Treated") %in% names(tb)) || any(tb < min_cells)) next
    de <- FindMarkers(sub, ident.1 = "Treated", ident.2 = "Control",
                      features = gene, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
    if (!nrow(de)) next
    de$gene <- rownames(de)
    de$log2FC <- if ("avg_log2FC" %in% names(de)) de$avg_log2FC else
      if ("avg_logFC" %in% names(de)) de$avg_logFC else NA_real_
    de$padj <- if ("p_val_adj" %in% names(de)) de$p_val_adj else p.adjust(de$p_val, "BH")
    de$cluster <- cl
    out[[cl]] <- de[, c("gene","cluster","log2FC","p_val","padj","pct.1","pct.2"), drop=FALSE]
  }
  if (length(out)) dplyr::bind_rows(out) else data.frame()
}

violin_plot_adv <- function(so, gene, show_median = TRUE, show_box = TRUE, zscore = FALSE, trim = TRUE) {
  df <- fetch_expr_df(so, gene, zscore = zscore)
  p <- ggplot(df, aes(x = cluster, y = expr, fill = CMO_label)) +
    geom_violin(position = position_dodge(width = 0.9), trim = trim) +
    labs(title = paste0("Violin: ", gene, " by cluster (split: Treated/Control",
                        if (zscore) ", Z-scored" else "", ")"),
         x = "Cluster", y = ifelse(zscore, "Z-score", "Normalized expression")) +
    theme_minimal()
  if (isTRUE(show_box)) {
    p <- p + geom_boxplot(width = 0.12, outlier.shape = NA, position = position_dodge(width = 0.9))
  }
  if (isTRUE(show_median)) {
    p <- p + stat_summary(fun = median, geom = "crossbar", width = 0.6,
                          position = position_dodge(width = 0.9))
  }
  p
}

# ===== UI =====
ui <- navbarPage(
  title = "10x Single-cell Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  collapsible = TRUE,
  
  tabPanel("Data & QC",
           sidebarLayout(
             sidebarPanel(width = 4,
                          h5("Load 10x folders (absolute path is recommended) or upload a Seurat RDS"),
                          textInput("ctrl_dir", "Control 10x folder",
                                    value = "C:/Users/Admin/Desktop/arezou/ak3/ak2/total4/shctrl-7"),
                          textInput("trt_dir",  "Treated 10x folder",
                                    value = "C:/Users/Admin/Desktop/arezou/ak3/ak2/total4/vps36-7"),
                          tags$small("If these folders differ on your machine, edit them above."),
                          hr(),
                          checkboxInput("has_mux", "Has CMO Multiplex assay?", value = FALSE),
                          textInput("cmo_control", "CMO feature (Control)", value = "CMO310"),
                          textInput("cmo_treated", "CMO feature (Treated)", value = "CMO312"),
                          hr(),
                          numericInput("minfeat", "min features", value = 200, min = 0, step = 50),
                          numericInput("maxfeat", "max features", value = 7000, min = 500, step = 500),
                          numericInput("maxmt",   "max %MT",    value = 15, min = 0, step = 1),
                          hr(),
                          numericInput("dims", "PCA dims", value = 30, min = 10, step = 5),
                          numericInput("res",  "Clustering resolution", value = 0.5, min = 0.1, step = 0.1),
                          actionButton("run_all", "Run pipeline", class = "btn btn-primary"),
                          br(), br(),
                          fileInput("seurat_rds", "Or upload a Seurat .rds (optional)", accept = ".rds")
             ),
             mainPanel(
               h4("QC summaries"),
               verbatimTextOutput("log"),
               fluidRow(
                 column(6, plotOutput("umap_treatment", height = 350)),
                 column(6, plotOutput("umap_cluster",   height = 350))
               ),
               fluidRow(
                 column(6, plotOutput("comp_within_cluster",   height = 320)),
                 column(6, plotOutput("comp_within_condition", height = 320))
               )
             )
           )
  ),
  
  tabPanel("VPS36 & Stats",
           sidebarLayout(
             sidebarPanel(width = 4,
                          textInput("vps_gene", "VPS36 gene symbol (auto-detect if blank)", value = ""),
                          helpText("Human: VPS36 — Mouse: Vps36 (auto-detects if blank)"),
                          actionButton("refresh_vps", "Refresh VPS36 plots")
             ),
             mainPanel(
               fluidRow(
                 column(6, plotOutput("vps_pos",  height = 320)),
                 column(6, plotOutput("vps_mean", height = 320))
               ),
               h5("Tables"),
               DTOutput("tbl_comp_cluster"),
               DTOutput("tbl_comp_condition"),
               DTOutput("tbl_vps")
             )
           )
  ),
  
  tabPanel("DE per Cluster",
           sidebarLayout(
             sidebarPanel(width = 4,
                          numericInput("min_cells", "Min cells per group (cluster × condition)", value = 10, min = 1, step = 1),
                          numericInput("lfc", "log2FC threshold (|.|)", value = 0.25, min = 0, step = 0.05),
                          numericInput("minpct", "min.pct", value = 0.10, min = 0, step = 0.05),
                          actionButton("run_de", "Run per-cluster DE", class = "btn btn-warning")
             ),
             mainPanel(
               h5("Top significant genes per cluster"),
               DTOutput("tbl_top"),
               plotOutput("heatmap_top", height = 700),
               downloadButton("dl_de_excel", "Download DE Excel")
             )
           )
  ),
  
  # ===== Gene-set from files =====
  tabPanel("Gene-set DE",
           sidebarLayout(
             sidebarPanel(width = 4,
                          fileInput("genesets", "Upload gene-set files (.csv/.xlsx)", multiple = TRUE,
                                    accept = c(".csv",".xlsx",".xls")),
                          selectInput("species_gs", "Species (for gene mapping)",
                                      choices = c("Human","Mouse","Other/Unknown"), selected = "Human"),
                          textInput("geneset_filter", "Run only files containing (optional)", value = ""),
                          checkboxInput("combine_sets", "Also create a combined set of ALL uploaded genes", value = TRUE),
                          hr(),
                          h6("Violin display options"),
                          checkboxInput("opt_median_gs", "Show median line", value = TRUE),
                          checkboxInput("opt_box_gs",    "Overlay boxplot",  value = TRUE),
                          checkboxInput("opt_z_gs",      "Z-score (display only)", value = FALSE),
                          checkboxInput("opt_trim_gs",   "Trim tails", value = TRUE),
                          numericInput("opt_min_cells_gs", "Min cells for one-gene DE", value = 10, min = 1, step = 1),
                          actionButton("run_gs", "Run Gene-set DE", class = "btn btn-success")
             ),
             mainPanel(uiOutput("gs_results_ui"))
           )
  ),
  
  # ===== Manual Gene-set =====
  tabPanel("Manual gene-set",
           sidebarLayout(
             sidebarPanel(width = 4,
                          textAreaInput("manual_genes","Enter 10–30 gene symbols (one per line or comma/semicolon/space-separated)",
                                        rows = 8, placeholder = "VPS36\nTP53\nEGFR\n..."),
                          selectInput("species_manual", "Species",
                                      choices = c("Human","Mouse","Other/Unknown"), selected = "Human"),
                          numericInput("manual_min_cells", "Min cells per group (cluster × condition)", value = 10, min = 1, step = 1),
                          hr(),
                          h6("Violin display options"),
                          checkboxInput("opt_median_m", "Show median line", value = TRUE),
                          checkboxInput("opt_box_m",    "Overlay boxplot",  value = TRUE),
                          checkboxInput("opt_z_m",      "Z-score (display only)", value = FALSE),
                          checkboxInput("opt_trim_m",   "Trim tails", value = TRUE),
                          actionButton("run_manual_gs", "Run manual gene-set DE", class = "btn btn-success")
             ),
             mainPanel(
               fluidRow(
                 column(6, DTOutput("manual_diag")),
                 column(6, DTOutput("manual_matched"))
               ),
               h5("Per-cluster DE (only your genes)"),
               DTOutput("manual_tbl_pc"),
               plotOutput("manual_volcano", height = 600),
               hr(),
               uiOutput("manual_violin_ui"),
               downloadButton("dl_manual_gene_stats", "Download selected gene: Stats + DE (Excel)")
             )
           )
  ),
  
  tabPanel("Downloads",
           fluidRow(
             column(4, downloadButton("dl_comp_vps_excel", "Compositions + VPS36 (Excel)")),
             column(4, downloadButton("dl_manifest",        "Manifest CSV")),
             column(4, downloadButton("dl_seurat_rds",      "Download Seurat RDS"))
           ),
           br(),
           verbatimTextOutput("sessioninfo")
  )
)

# ===== Server =====
server <- function(input, output, session) {
  
  log_txt <- reactiveVal("")
  add_log <- function(...) log_txt(paste(log_txt(), paste(..., collapse = " "), sep = "\n"))
  output$log <- renderText(log_txt())
  
  # ---------- Load & Process ----------
  merged_obj <- eventReactive(input$run_all, {
    withProgress(message = "Running pipeline…", value = 0, {
      incProgress(0.05, detail = "Loading data")
      if (!is.null(input$seurat_rds)) {
        so <- readRDS(input$seurat_rds$datapath)
        add_log(sprintf("Loaded Seurat RDS: %d cells", ncol(so)))
      } else {
        ctrl_dir <- norm_path(input$ctrl_dir); trt_dir <- norm_path(input$trt_dir)
        validate(need(dir.exists(ctrl_dir),  paste0("Control folder not found: ", ctrl_dir)))
        validate(need(dir.exists(trt_dir),   paste0("Treated folder not found: ", trt_dir)))
        
        ctrl_data_raw    <- Read10X(data.dir = ctrl_dir)
        treated_data_raw <- Read10X(data.dir = trt_dir)
        ctrl_rna    <- get_assay_layer(ctrl_data_raw)
        treated_rna <- get_assay_layer(treated_data_raw)
        
        raw_rna <- cbind(ctrl_rna, treated_rna)
        so <- CreateSeuratObject(counts = raw_rna, project = "Merged", min.cells = 3, min.features = 200)
        
        has_mux <- isTRUE(input$has_mux) && is.list(ctrl_data_raw) && is.list(treated_data_raw) &&
          ("Multiplexing Capture" %in% names(ctrl_data_raw)) &&
          ("Multiplexing Capture" %in% names(treated_data_raw))
        if (has_mux) {
          raw_mux <- cbind(ctrl_data_raw[["Multiplexing Capture"]], treated_data_raw[["Multiplexing Capture"]])
          so[["Multiplex"]] <- CreateAssayObject(counts = raw_mux)
        }
        
        DefaultAssay(so) <- "RNA"
        if (has_mux && ("Multiplex" %in% names(so@assays))) {
          DefaultAssay(so) <- "Multiplex"
          mux_feats <- rownames(so[["Multiplex"]])
          if (all(c(input$cmo_control, input$cmo_treated) %in% mux_feats)) {
            cA <- GetAssayData(so, slot = "counts")[input$cmo_control, ]
            cB <- GetAssayData(so, slot = "counts")[input$cmo_treated, ]
            so$CMO_ratio <- cB / (cA + cB + 1e-6)
            so$CMO_label <- ifelse(so$CMO_ratio > 0.6, "Treated", "Control")
          } else {
            warning("CMO features not found. Falling back to column membership.")
            so$CMO_label <- ifelse(colnames(so) %in% colnames(treated_rna), "Treated", "Control")
          }
          DefaultAssay(so) <- "RNA"
        } else {
          so$CMO_label <- ifelse(colnames(so) %in% colnames(treated_rna), "Treated", "Control")
        }
        add_log("Merged cells:", ncol(so))
      }
      
      incProgress(0.25, detail = "QC & Normalize")
      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
      so <- subset(so, subset = nFeature_RNA >= input$minfeat & nFeature_RNA <= input$maxfeat & percent.mt <= input$maxmt)
      so <- NormalizeData(so)
      so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
      so <- ScaleData(so, features = VariableFeatures(so))
      
      incProgress(0.60, detail = "PCA/Neighbors/Clusters/UMAP")
      so <- RunPCA(so, features = VariableFeatures(so), verbose = FALSE)
      so <- FindNeighbors(so, dims = 1:input$dims)
      so <- FindClusters(so, resolution = input$res)
      so <- RunUMAP(so, dims = 1:input$dims)
      
      add_log("Done: QC + clustering + UMAP")
      so
    })
  })
  
  # ---------- UMAP & Compositions ----------
  output$umap_treatment <- renderPlot({ req(merged_obj()); DimPlot(merged_obj(), group.by = "CMO_label", pt.size = 0.5) + ggtitle("UMAP by Treatment") })
  output$umap_cluster   <- renderPlot({ req(merged_obj()); DimPlot(merged_obj(), label = TRUE, pt.size = 0.5) + ggtitle("UMAP by Cluster") })
  
  comp_within_cluster <- reactive({
    req(merged_obj())
    meta <- merged_obj()@meta.data
    if (!"seurat_clusters" %in% names(meta)) meta$seurat_clusters <- as.character(Idents(merged_obj()))
    meta$cluster <- factor(meta$seurat_clusters)
    meta %>% count(cluster, CMO_label, name = "Count") %>%
      group_by(cluster) %>% mutate(Percent = 100 * Count / sum(Count)) %>% ungroup() %>%
      tidyr::complete(cluster, CMO_label, fill = list(Count = 0, Percent = 0))
  })
  comp_within_condition <- reactive({
    req(merged_obj())
    meta <- merged_obj()@meta.data
    if (!"seurat_clusters" %in% names(meta)) meta$seurat_clusters <- as.character(Idents(merged_obj()))
    meta$cluster <- factor(meta$seurat_clusters)
    meta %>% count(CMO_label, cluster, name = "Count") %>%
      group_by(CMO_label) %>% mutate(Percent = 100 * Count / sum(Count)) %>% ungroup() %>%
      tidyr::complete(CMO_label, cluster, fill = list(Count = 0, Percent = 0))
  })
  ord_levels <- reactive({
    df <- comp_within_cluster()
    df %>% filter(CMO_label == "Treated") %>% arrange(desc(Percent)) %>% pull(cluster) %>% as.character()
  })
  output$comp_within_cluster <- renderPlot({
    df <- comp_within_cluster(); if (length(ord_levels()) > 0) df$cluster <- factor(df$cluster, levels = ord_levels())
    ggplot(df, aes(x = cluster, y = Percent, fill = CMO_label)) +
      geom_col(position = position_dodge(width = 0.9)) +
      labs(title = "Cluster Composition (% within cluster)", x = "Cluster", y = "Percent") + theme_minimal()
  })
  output$comp_within_condition <- renderPlot({
    df <- comp_within_condition(); if (length(ord_levels()) > 0) df$cluster <- factor(df$cluster, levels = ord_levels())
    ggplot(df, aes(x = CMO_label, y = Percent, fill = cluster)) +
      geom_col(position = "stack") +
      labs(title = "Cluster Distribution (within condition)", x = "Condition", y = "Percent") + theme_minimal()
  })
  output$tbl_comp_cluster   <- renderDT({ as_tbl(comp_within_cluster()) })
  output$tbl_comp_condition <- renderDT({ as_tbl(comp_within_condition()) })
  
  # ---------- VPS36 ----------
  vps_stats <- eventReactive({ input$refresh_vps; merged_obj() }, {
    req(merged_obj()); so <- merged_obj(); DefaultAssay(so) <- "RNA"
    all_genes <- rownames(so); vg <- input$vps_gene
    if (is.null(vg) || vg == "") {
      if ("VPS36" %in% all_genes) vg <- "VPS36" else {
        cand <- all_genes[grepl("^VPS36$", all_genes, ignore.case = TRUE)]
        validate(need(length(cand) > 0, "VPS36 not found in features")); vg <- cand[1]
      }
    } else validate(need(vg %in% all_genes, paste0("Gene ", vg, " not found in object")))
    md <- so@meta.data
    md$cluster <- factor(if (!"seurat_clusters" %in% names(md)) as.character(Idents(so)) else md$seurat_clusters)
    md$VPS36_expr <- FetchData(so, vars = vg)[,1]; md$VPS36_pos <- md$VPS36_expr > 0
    vps_df <- md %>% group_by(cluster, CMO_label) %>%
      summarize(Total = n(), Pos = sum(VPS36_pos, na.rm = TRUE),
                Percent = 100 * Pos / Total, MeanExp = mean(VPS36_expr, na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(cluster, CMO_label, fill = list(Total=0, Pos=0, Percent=0, MeanExp=0))
    list(vps_df = vps_df, gene = vg)
  })
  output$vps_pos  <- renderPlot({ req(vps_stats()); df <- vps_stats()$vps_df; if (length(ord_levels())>0) df$cluster <- factor(df$cluster, levels=ord_levels());
  ggplot(df, aes(x=cluster,y=Percent,fill=CMO_label)) + geom_col(position=position_dodge(width=.9)) +
      labs(title=paste0("VPS36+ rate by cluster (gene: ", vps_stats()$gene, ")"), x="Cluster", y="% positive") + theme_minimal() })
  output$vps_mean <- renderPlot({ req(vps_stats()); df <- vps_stats()$vps_df; if (length(ord_levels())>0) df$cluster <- factor(df$cluster, levels=ord_levels());
  ggplot(df, aes(x=cluster,y=MeanExp,fill=CMO_label)) + geom_col(position=position_dodge(width=.9)) +
      labs(title=paste0("Mean expression of ", vps_stats()$gene, " per cluster × condition"), x="Cluster", y="Mean expr") + theme_minimal() })
  output$tbl_vps <- renderDT({ req(vps_stats()); as_tbl(vps_stats()$vps_df) })
  
  # ---------- per-cluster DE (all genes) ----------
  de_results <- eventReactive(input$run_de, {
    req(merged_obj()); so <- merged_obj(); DefaultAssay(so) <- "RNA"
    md <- so@meta.data; md$Cell <- rownames(md)
    if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
    md$cluster <- factor(md$seurat_clusters); clusters <- levels(md$cluster)
    de_list <- list()
    withProgress(message = "Running per-cluster DE", value = 0, {
      for (i in seq_along(clusters)) {
        cl <- clusters[i]; incProgress(1/length(clusters), detail = paste("cluster", cl))
        cells_cl <- md$Cell[md$cluster == cl]
        subobj   <- subset(so, cells = cells_cl); Idents(subobj) <- subobj$CMO_label
        tb <- table(subobj$CMO_label)
        if (!all(c("Control","Treated") %in% names(tb)) || any(tb < input$min_cells)) next
        de <- FindMarkers(subobj, ident.1 = "Treated", ident.2 = "Control",
                          logfc.threshold = input$lfc, min.pct = input$minpct, test.use = "wilcox")
        if (nrow(de) == 0) next
        de$gene   <- rownames(de)
        de$log2FC <- if ("avg_log2FC" %in% colnames(de)) de$avg_log2FC else if ("avg_logFC" %in% colnames(de)) de$avg_logFC else NA_real_
        de$padj   <- if ("p_val_adj" %in% colnames(de)) de$p_val_adj else p.adjust(de$p_val, "BH")
        de$cluster <- cl
        de$deg    <- (de$padj < 0.05) & (abs(de$log2FC) >= input$lfc)
        de <- de[order(de$padj, -abs(de$log2FC)), ]
        de_list[[cl]] <- de
      }
    })
    if (!length(de_list)) return(NULL)
    de_all <- dplyr::bind_rows(de_list)
    de_sig <- subset(de_all, deg)
    top_by_cluster <- de_sig %>% group_by(cluster) %>% arrange(padj, .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()
    list(de_all = de_all, de_sig = de_sig, top_by_cluster = top_by_cluster)
  })
  output$tbl_top <- renderDT({ req(de_results()); as_tbl(de_results()$top_by_cluster) })
  output$heatmap_top <- renderPlot({
    req(de_results()); de_all <- de_results()$de_all; sel <- unique(de_results()$top_by_cluster$gene)
    validate(need(length(sel) > 1, "Not enough significant genes to draw heatmap"))
    mat_wide <- reshape2::acast(de_all[de_all$gene %in% sel, ], gene ~ cluster, value.var = "log2FC",
                                fun.aggregate = function(z) if (length(z)) mean(z, na.rm = TRUE) else NA_real_)
    mat_wide[is.na(mat_wide)] <- 0; pheatmap(mat_wide)
  })
  output$dl_de_excel <- downloadHandler(
    filename = function() "DE_Treated_vs_Control_byCluster.xlsx",
    content  = function(file) { req(de_results()); openxlsx::write.xlsx(list(
      ALL = de_results()$de_all, SIGNIFICANT = de_results()$de_sig, TOP_perCluster = de_results()$top_by_cluster), file) }
  )
  
  # ---------- Gene-set DE (uploads) ----------
  geneset_results <- eventReactive(input$run_gs, {
    req(merged_obj()); files <- input$genesets
    validate(need(!is.null(files), "Please upload at least one gene-set file"))
    if (nzchar(input$geneset_filter)) files <- files[grepl(input$geneset_filter, files$name, ignore.case = TRUE), , drop = FALSE]
    validate(need(nrow(files) > 0, "No matching files after filter"))
    
    so <- merged_obj(); DefaultAssay(so) <- "RNA"; Idents(so) <- so$CMO_label
    universe <- rownames(so); species <- input$species_gs
    cols_pref <- c("gene","log2FC","p_val","padj","p_val_adj","pct.1","pct.2","deg")
    out <- list(); union_genes <- character(0)
    
    withProgress(message = "Running gene-set DE", value = 0, {
      for (i in seq_len(nrow(files))) {
        incProgress(1/nrow(files), detail = files$name[i])
        fpath <- files$datapath[i]; gname <- tools::file_path_sans_ext(files$name[i]); gid <- safe_id(gname)
        
        req_genes <- read_genes_flex(fpath, universe, species = species); if (!length(req_genes)) next
        map_genes <- match_genes(req_genes, universe); union_genes <- unique(c(union_genes, map_genes))
        if (!length(map_genes)) next
        
        de_g <- FindMarkers(so, ident.1 = "Treated", ident.2 = "Control",
                            features = map_genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
        if (!nrow(de_g)) next
        de_g$gene   <- rownames(de_g)
        de_g$log2FC <- if ("avg_log2FC" %in% names(de_g)) de_g$avg_log2FC else if ("avg_logFC" %in% names(de_g)) de_g$avg_logFC else NA_real_
        de_g$padj   <- if ("p_val_adj" %in% names(de_g)) de_g$p_val_adj else p.adjust(de_g$p_val, "BH")
        de_g$deg    <- (de_g$padj < 0.05) & (abs(de_g$log2FC) >= 0.25)
        de_g        <- de_g[order(de_g$padj, -abs(de_g$log2FC)), ]
        
        vol <- ggplot(de_g, aes(x = log2FC, y = -log10(padj))) +
          geom_point(aes(shape = deg), size = 1.8, alpha = 0.9) +
          geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
          ggrepel::geom_text_repel(
            data = head(subset(de_g, deg)[order(subset(de_g, deg)$padj), ], 20),
            aes(label = gene), size = 3, max.overlaps = Inf
          ) +
          labs(title = paste0("Volcano (", gname, ") — Treated vs Control"), x = "log2FC", y = "-log10(adj.p)") +
          theme_minimal()
        
        # per-cluster on same genes
        md <- so@meta.data
        if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
        md$cluster <- factor(md$seurat_clusters); cl_levels <- levels(md$cluster)
        de_list_gs <- list()
        for (cl in cl_levels) {
          cells_cl <- rownames(md)[md$cluster == cl]
          subobj   <- subset(so, cells = cells_cl); Idents(subobj) <- subobj$CMO_label
          tb <- table(subobj$CMO_label); if (!all(c("Control","Treated") %in% names(tb)) || any(tb < 10)) next
          de_cl <- FindMarkers(subobj, ident.1 = "Treated", ident.2 = "Control",
                               features = map_genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
          if (!nrow(de_cl)) next
          de_cl$gene   <- rownames(de_cl)
          de_cl$log2FC <- if ("avg_log2FC" %in% names(de_cl)) de_cl$avg_log2FC else if ("avg_logFC" %in% names(de_cl)) de_cl$avg_logFC else NA_real_
          de_cl$padj   <- if ("p_val_adj" %in% names(de_cl)) de_cl$p_val_adj else p.adjust(de_cl$p_val, "BH")
          de_cl$cluster <- cl
          de_list_gs[[cl]] <- de_cl
        }
        
        mat_wide_df <- data.frame(); hm <- NULL; de_pc_all <- data.frame()
        if (length(de_list_gs) > 0) {
          de_pc_all <- dplyr::bind_rows(de_list_gs); de_pc_all$absLFC <- abs(de_pc_all$log2FC)
          top_genes <- de_pc_all %>% group_by(gene) %>%
            summarize(score = median(absLFC, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(score)) %>% slice_head(n = 50) %>% pull(gene)
          dd <- de_pc_all[de_pc_all$gene %in% top_genes, ]
          if (nrow(dd)) {
            mat_wide <- reshape2::acast(dd, gene ~ cluster, value.var = "log2FC",
                                        fun.aggregate = function(z) if (length(z)) mean(z, na.rm = TRUE) else NA_real_)
            mat_wide_df <- data.frame(gene = rownames(mat_wide), mat_wide, check.names = FALSE)
            m2 <- mat_wide; m2[is.na(m2)] <- 0; hm <- pheatmap::pheatmap(m2)
          }
        }
        
        out[[gid]] <- list(
          title   = gname, id = gid, genes = unique(de_g$gene),
          global  = de_g, global_tab = dplyr::select(de_g, dplyr::any_of(cols_pref)),
          volcano = vol, percluster = de_pc_all, mat = mat_wide_df, heatmap = hm,
          diag    = data.frame(file = gname, species = species,
                               read_n = length(req_genes), match_n = length(map_genes),
                               tested = nrow(de_g), sig_n  = sum(de_g$deg, na.rm = TRUE))
        )
      }
    })
    
    # union set (optional)
    if (isTRUE(input$combine_sets) && length(union_genes) > 0) {
      gname <- "ALL uploaded genes"; gid <- safe_id(gname); map_genes <- union_genes
      de_g <- FindMarkers(so, ident.1 = "Treated", ident.2 = "Control",
                          features = map_genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
      if (nrow(de_g) > 0) {
        de_g$gene   <- rownames(de_g)
        de_g$log2FC <- if ("avg_log2FC" %in% names(de_g)) de_g$avg_log2FC else if ("avg_logFC" %in% names(de_g)) de_g$avg_logFC else NA_real_
        de_g$padj   <- if ("p_val_adj" %in% names(de_g)) de_g$p_val_adj else p.adjust(de_g$p_val, "BH")
        de_g$deg    <- (de_g$padj < 0.05) & (abs(de_g$log2FC) >= 0.25)
        de_g        <- de_g[order(de_g$padj, -abs(de_g$log2FC)), ]
        vol <- ggplot(de_g, aes(x = log2FC, y = -log10(padj))) +
          geom_point(aes(shape = deg), size = 1.8, alpha = 0.9) +
          geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
          ggrepel::geom_text_repel(
            data = head(subset(de_g, deg)[order(subset(de_g, deg)$padj), ], 20),
            aes(label = gene), size = 3, max.overlaps = Inf
          ) +
          labs(title = paste0("Volcano (", gname, ") — Treated vs Control"), x = "log2FC", y = "-log10(adj.p)") +
          theme_minimal()
        
        md <- so@meta.data
        if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
        md$cluster <- factor(md$seurat_clusters); cl_levels <- levels(md$cluster)
        de_list_gs <- list()
        for (cl in cl_levels) {
          cells_cl <- rownames(md)[md$cluster == cl]
          subobj   <- subset(so, cells = cells_cl); Idents(subobj) <- subobj$CMO_label
          tb <- table(subobj$CMO_label); if (!all(c("Control","Treated") %in% names(tb)) || any(tb < 10)) next
          de_cl <- FindMarkers(subobj, ident.1 = "Treated", ident.2 = "Control",
                               features = map_genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
          if (!nrow(de_cl)) next
          de_cl$gene   <- rownames(de_cl)
          de_cl$log2FC <- if ("avg_log2FC" %in% names(de_cl)) de_cl$avg_log2FC else if ("avg_logFC" %in% names(de_cl)) de_cl$avg_logFC else NA_real_
          de_cl$padj   <- if ("p_val_adj" %in% names(de_cl)) de_cl$p_val_adj else p.adjust(de_cl$p_val, "BH")
          de_cl$cluster <- cl
          de_list_gs[[cl]] <- de_cl
        }
        
        mat_wide_df <- data.frame(); hm <- NULL; de_pc_all <- data.frame()
        if (length(de_list_gs) > 0) {
          de_pc_all <- dplyr::bind_rows(de_list_gs); de_pc_all$absLFC <- abs(de_pc_all$log2FC)
          top_genes <- de_pc_all %>% group_by(gene) %>%
            summarize(score = median(absLFC, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(score)) %>% slice_head(n = 50) %>% pull(gene)
          dd <- de_pc_all[de_pc_all$gene %in% top_genes, ]
          if (nrow(dd)) {
            mat_wide <- reshape2::acast(dd, gene ~ cluster, value.var = "log2FC",
                                        fun.aggregate = function(z) if (length(z)) mean(z, na.rm = TRUE) else NA_real_)
            mat_wide_df <- data.frame(gene = rownames(mat_wide), mat_wide, check.names = FALSE)
            m2 <- mat_wide; m2[is.na(m2)] <- 0; hm <- pheatmap::pheatmap(m2)
          }
        }
        
        out[[gid]] <- list(
          title = gname, id = gid, genes = unique(de_g$gene),
          global = de_g, global_tab = dplyr::select(de_g, dplyr::any_of(cols_pref)),
          volcano = vol, percluster = de_pc_all, mat = mat_wide_df, heatmap = hm,
          diag = data.frame(file = gname, species = species,
                            read_n = length(union_genes), match_n = length(map_genes),
                            tested = nrow(de_g), sig_n = sum(de_g$deg, na.rm = TRUE))
        )
      }
    }
    out
  })
  
  # Dynamic UI per gene-set tab, with advanced violin + stats + one-gene DE + download
  output$gs_results_ui <- renderUI({
    res <- geneset_results()
    if (is.null(res) || !length(res)) {
      return(h5("No gene-set results to show. Make sure you ran the pipeline and uploaded valid gene-set files (first column = gene symbols)."))
    }
    tabs <- lapply(res, function(gs){
      tabPanel(gs$title,
               h5("Global DE"), DTOutput(paste0("gs_tbl_", gs$id)),
               plotOutput(paste0("gs_vol_", gs$id), height = 380),
               h5("Per-cluster DE (if available)"), DTOutput(paste0("gs_tbl_pc_", gs$id)),
               plotOutput(paste0("gs_hm_", gs$id), height = 600),
               hr(), h5("Violin by cluster × condition (pick a gene)"),
               uiOutput(paste0("gs_violin_ui_", gs$id)),
               DTOutput(paste0("gs_violin_stats_", gs$id)),
               DTOutput(paste0("gs_violin_de_", gs$id)),
               downloadButton(paste0("dl_gs_gene_", gs$id), "Download selected gene: Stats + DE (Excel)"),
               hr(), h6("Diagnostics"), DTOutput(paste0("gs_diag_", gs$id)),
               downloadButton(paste0("dl_gs_", gs$id), "Download gene-set Excel")
      )
    })
    do.call(tabsetPanel, c(list(id = "gs_tabs"), unname(tabs)))
  })
  
  observe({
    res <- geneset_results(); if (is.null(res) || !length(res)) return()
    so <- merged_obj()
    for (nm in names(res)) {
      local({
        gs <- res[[nm]]
        output[[paste0("gs_tbl_", gs$id)]]    <- renderDT({ as_tbl(gs$global_tab) })
        output[[paste0("gs_vol_", gs$id)]]    <- renderPlot({ gs$volcano })
        output[[paste0("gs_tbl_pc_", gs$id)]] <- renderDT({ as_tbl(gs$percluster) })
        output[[paste0("gs_hm_", gs$id)]]     <- renderPlot({
          if (!is.null(gs$heatmap)) { grid::grid.newpage(); grid::grid.draw(gs$heatmap$gtable) }
          else { plot.new(); title("No per-cluster DE / heatmap for this set") }
        })
        output[[paste0("gs_diag_", gs$id)]]   <- renderDT({ as_tbl(gs$diag) })
        
        output[[paste0("gs_violin_ui_", gs$id)]] <- renderUI({
          tagList(
            selectizeInput(paste0("gs_gene_", gs$id), "Gene:",
                           choices = sort(gs$genes), multiple = FALSE,
                           options = list(placeholder = 'Select a gene')),
            plotOutput(paste0("gs_violin_", gs$id), height = 420)
          )
        })
        
        output[[paste0("gs_violin_", gs$id)]] <- renderPlot({
          req(input[[paste0("gs_gene_", gs$id)]])
          g <- input[[paste0("gs_gene_", gs$id)]]
          validate(need(g %in% rownames(so), paste0("Gene ", g, " not found in object")))
          violin_plot_adv(
            so, g,
            show_median = isTRUE(input$opt_median_gs),
            show_box    = isTRUE(input$opt_box_gs),
            zscore      = isTRUE(input$opt_z_gs),
            trim        = isTRUE(input$opt_trim_gs)
          )
        })
        
        # stats table + one-gene DE
        output[[paste0("gs_violin_stats_", gs$id)]] <- renderDT({
          req(input[[paste0("gs_gene_", gs$id)]])
          as_tbl(violin_stats_tbl(so, input[[paste0("gs_gene_", gs$id)]]))
        })
        output[[paste0("gs_violin_de_", gs$id)]] <- renderDT({
          req(input[[paste0("gs_gene_", gs$id)]])
          as_tbl(one_gene_de_per_cluster(so, input[[paste0("gs_gene_", gs$id)]],
                                         min_cells = input$opt_min_cells_gs))
        })
        
        # download (selected gene)
        output[[paste0("dl_gs_gene_", gs$id)]] <- downloadHandler(
          filename = function() paste0("SelectedGene_", gs$id, "_", input[[paste0("gs_gene_", gs$id)]], ".xlsx"),
          content  = function(file) {
            g <- input[[paste0("gs_gene_", gs$id)]]
            wb <- list(
              Stats_byCluster = violin_stats_tbl(so, g),
              DE_oneGene_perCluster = one_gene_de_per_cluster(so, g, min_cells = input$opt_min_cells_gs)
            )
            openxlsx::write.xlsx(wb, file)
          }
        )
        
        # download (whole gene-set)
        output[[paste0("dl_gs_", gs$id)]] <- downloadHandler(
          filename = function() paste0("GeneSet_DE_Summary_", gs$id, ".xlsx"),
          content  = function(file) {
            wb_list <- list(Global_DE = gs$global_tab)
            if (nrow(gs$percluster) > 0) wb_list[["PerCluster_DE"]] <- gs$percluster
            if (nrow(gs$mat) > 0)        wb_list[["log2FC_Matrix"]] <- gs$mat
            openxlsx::write.xlsx(wb_list, file)
          }
        )
      })
    }
  })
  
  # ---------- Manual gene-set ----------
  manual_gs_results <- eventReactive(input$run_manual_gs, {
    req(merged_obj()); so <- merged_obj(); DefaultAssay(so) <- "RNA"
    species <- input$species_manual
    req_genes <- parse_gene_text(input$manual_genes, species = species)
    validate(need(length(req_genes) > 0, "No genes provided"))
    universe  <- rownames(so)
    map_genes <- match_genes(req_genes, universe)
    validate(need(length(map_genes) > 0, "None of the provided genes matched features in the object"))
    
    md <- so@meta.data
    if (!"seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(Idents(so))
    md$cluster <- factor(md$seurat_clusters); cl_levels  <- levels(md$cluster)
    
    de_list <- list()
    for (cl in cl_levels) {
      cells_cl <- rownames(md)[md$cluster == cl]
      subobj   <- subset(so, cells = cells_cl); Idents(subobj) <- subobj$CMO_label
      tb <- table(subobj$CMO_label); if (!all(c("Control","Treated") %in% names(tb)) || any(tb < input$manual_min_cells)) next
      de_cl <- FindMarkers(subobj, ident.1 = "Treated", ident.2 = "Control",
                           features = map_genes, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
      if (!nrow(de_cl)) next
      de_cl$gene   <- rownames(de_cl)
      de_cl$log2FC <- if ("avg_log2FC" %in% names(de_cl)) de_cl$avg_log2FC else if ("avg_logFC" %in% names(de_cl)) de_cl$avg_logFC else NA_real_
      de_cl$padj   <- if ("p_val_adj" %in% names(de_cl)) de_cl$p_val_adj else p.adjust(de_cl$p_val, "BH")
      de_cl$cluster <- cl
      de_list[[cl]] <- de_cl
    }
    de_pc_all <- if (length(de_list)) dplyr::bind_rows(de_list) else data.frame()
    if (nrow(de_pc_all)) de_pc_all$deg <- (de_pc_all$padj < 0.05) & (abs(de_pc_all$log2FC) >= 0.25)
    list(requested = req_genes, matched = map_genes, percluster = de_pc_all)
  })
  
  output$manual_diag    <- renderDT({ res <- manual_gs_results(); as_tbl(data.frame(requested_n = length(res$requested), matched_n = length(res$matched))) })
  output$manual_matched <- renderDT({ res <- manual_gs_results(); as_tbl(data.frame(matched_genes = res$matched)) })
  output$manual_tbl_pc  <- renderDT({ res <- manual_gs_results(); as_tbl(res$percluster) })
  output$manual_volcano <- renderPlot({
    res <- manual_gs_results(); df <- res$percluster
    validate(need(nrow(df) > 0, "No per-cluster DE results (insufficient cells or no matches)."))
    top_lab <- df %>% dplyr::filter(deg) %>% dplyr::group_by(cluster) %>% dplyr::arrange(padj, .by_group = TRUE) %>%
      dplyr::slice_head(n = 10) %>% dplyr::ungroup()
    ggplot(df, aes(x = log2FC, y = -log10(padj))) +
      geom_point(aes(shape = deg), size = 1.6, alpha = 0.9) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 3, max.overlaps = Inf) +
      facet_wrap(~ cluster, scales = "free") +
      labs(title = "Volcano per cluster — Treated vs Control (manual gene-set)",
           x = "log2FC", y = "-log10(adj.p)") + theme_minimal()
  })
  
  output$manual_violin_ui <- renderUI({
    res <- manual_gs_results(); validate(need(length(res$matched) > 0, ""))
    tagList(
      selectizeInput("manual_gene_pick", "Gene for violin", choices = sort(res$matched), multiple = FALSE),
      plotOutput("manual_violin", height = 420),
      DTOutput("manual_violin_stats"),
      DTOutput("manual_violin_de")
    )
  })
  output$manual_violin <- renderPlot({
    req(merged_obj()); g <- input$manual_gene_pick; validate(need(!is.null(g) && nzchar(g), "Select a gene"))
    so <- merged_obj(); validate(need(g %in% rownames(so), paste0("Gene ", g, " not found in object")))
    violin_plot_adv(
      so, g,
      show_median = isTRUE(input$opt_median_m),
      show_box    = isTRUE(input$opt_box_m),
      zscore      = isTRUE(input$opt_z_m),
      trim        = isTRUE(input$opt_trim_m)
    )
  })
  output$manual_violin_stats <- renderDT({
    req(merged_obj(), input$manual_gene_pick); so <- merged_obj()
    as_tbl(violin_stats_tbl(so, input$manual_gene_pick))
  })
  output$manual_violin_de <- renderDT({
    req(merged_obj(), input$manual_gene_pick); so <- merged_obj()
    as_tbl(one_gene_de_per_cluster(so, input$manual_gene_pick, min_cells = input$manual_min_cells))
  })
  output$dl_manual_gene_stats <- downloadHandler(
    filename = function() paste0("SelectedGene_manual_", input$manual_gene_pick, ".xlsx"),
    content  = function(file) {
      so <- merged_obj(); g <- input$manual_gene_pick
      wb <- list(
        Stats_byCluster = violin_stats_tbl(so, g),
        DE_oneGene_perCluster = one_gene_de_per_cluster(so, g, min_cells = input$manual_min_cells)
      )
      openxlsx::write.xlsx(wb, file)
    }
  )
  
  # ---------- Downloads ----------
  output$dl_comp_vps_excel <- downloadHandler(
    filename = function() "cluster_composition_and_VPS36.xlsx",
    content  = function(file) {
      req(comp_within_cluster(), comp_within_condition(), vps_stats())
      openxlsx::write.xlsx(list(
        `comp_within_cluster(%)`   = comp_within_cluster(),
        `comp_within_condition(%)` = comp_within_condition(),
        `VPS36_stats`              = vps_stats()$vps_df
      ), file)
    }
  )
  output$dl_manifest <- downloadHandler(
    filename = function() "results_manifest.csv",
    content  = function(file) {
      req(merged_obj())
      man <- data.frame(
        File = c("UMAP_by_treatment (panel)", "UMAP_by_cluster (panel)",
                 "Cluster_Composition_withinCluster (panel)", "Cluster_Distribution_withinCondition (panel)",
                 "VPS36 plots (panel)", "DE tables (panel)", "Geneset tabs (panel)"),
        Description = rep("Rendered in app", 7)
      )
      write.csv(man, file, row.names = FALSE)
    }
  )
  output$dl_seurat_rds <- downloadHandler(
    filename = function() "merged_seurat.rds",
    content  = function(file) { req(merged_obj()); saveRDS(merged_obj(), file) }
  )
  output$sessioninfo <- renderText({ paste(capture.output(sessionInfo()), collapse = "\n") })
}

shinyApp(ui, server)
