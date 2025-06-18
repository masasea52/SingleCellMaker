library(shiny)
library(Seurat)
library(Matrix)
library(ggplot2)
library(DT)
library(dplyr)
library(shinyjs)
library(decoupleR)
library(OmnipathR)
library(stringr)
library(viridis)
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(MAST)

options(shiny.maxRequestSize = 100 * 1024^2)

shinyServer(function(input, output, session) {
  
  initial_seurat_obj <- reactiveVal(NULL)
  seurat_obj <- reactiveVal(NULL)
  marker_result <- reactiveVal(NULL)
  heatmap_obj <- reactiveVal(NULL)
  tf_activities <- reactiveVal(NULL) 
  tf_top_tfs <- reactiveVal(NULL) 
  tf_pheatmap_mat <- reactiveVal(NULL) 
  tf_activities_df <- reactiveVal(NULL)
  
 
  collectri_net <- reactive({
    req(input$target_spiecies) 
    net <- get_collectri(organism = input$target_spiecies, split_complexes = FALSE)
    return(net)
  })
  
  observe({
    
    shinyjs::toggleState("run", !is.null(initial_seurat_obj()))
    shinyjs::toggleState("calc_marker", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_UMAP", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_feature", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_dotplot", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_vlnplot", !is.null(seurat_obj()))
    shinyjs::toggleState("dotsize", !is.null(seurat_obj()))
    shinyjs::toggleState("download_marker", !is.null(marker_result()))
    shinyjs::toggleState("top_n_genes", !is.null(seurat_obj()))
    shinyjs::toggleState("draw_heatmap", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_heatmap", !is.null(heatmap_obj()))
    shinyjs::toggleState("draw_tf_heatmap", !is.null(seurat_obj()) && !is.null(input$target_spiecies))
    shinyjs::toggleState("download_tf_table", !is.null(tf_activities_df()))
    shinyjs::toggleState("top_n_active_TF", !is.null(seurat_obj()))
    shinyjs::toggleState("preview_download_tf", !is.null(tf_pheatmap_mat()))
  })
  
  #Create % Load Demo Data
  observeEvent(input$create_and_load_demo_data, {
    withProgress(message = "Creating and loading demo data...", {
      demo_zip_path <- file.path("www", "demo_10x_data.zip")
      
      if(!file.exists(demo_zip_path)) {
        showModal(modalDialog(
          title = "Error",
          "Demo data zip file not found at 'www/demo_10x_data_zip'. Please ensure it exists",
          footer = modalButton("Close")
        ))
        return(NULL)
      }
      
      extract_dir <- tempfile()
      unzip(demo_zip_path, exdir = extract_dir)
      
      subdirs <- list.dirs(extract_dir, full.names = TRUE, recursive = FALSE)
      data_dir <- if (length(subdirs) == 1) subdirs[1] else extract_dir
      
      expr <- Read10X(data.dir = data_dir)
      seu <- CreateSeuratObject(counts = expr)
      
      DefaultAssay(seu) <- "RNA"
      
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")
      
      initial_seurat_obj(seu)
      showNotification("Demo data created and loaded successfully! Now click 'Run Seurat Analysis'.", type = "message")
    })
  })
  
  #Upload data & Create VlnPlot for QC
  observeEvent(input$run_initial_processing, {
    req(input$zipfile)
    withProgress(message = "Reading data and Preparing for QC...", {
      zip_path <- input$zipfile$datapath
      extract_dir <- tempfile()
      unzip(zip_path, exdir = extract_dir)
      
      subdirs <- list.dirs(extract_dir, full.names = TRUE, recursive = FALSE)
      data_dir <- if (length(subdirs) == 1) subdirs[1] else extract_dir
      
      expr <- Read10X(data.dir = data_dir)
      seu <- CreateSeuratObject(counts = expr)
      
      #if (ncol(seu) > 1000) { 
        #seu <- subset(seu, cells = sample(colnames(seu), size = 500, replace = FALSE))
      #}
      
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")
      
      initial_seurat_obj(seu)
    })
  })
  
  #QC Plot
  output$qc_plot_mt <- renderPlot({
    req(initial_seurat_obj())
    VlnPlot(initial_seurat_obj(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0.1)
  })
  
  #apply QC Check & Run Seurat
  observeEvent(input$run, {
    req(initial_seurat_obj())
    req(input$min_features, input$max_features, input$max_mt_percent)
    
    withProgress(message = "Running Seurat Analysis...", value = 0, {
      seu <- initial_seurat_obj()
      
      #QC filtering
      seu_filtered <- subset(seu,
                             subset = nFeature_RNA > input$min_features &
                               nFeature_RNA < input$max_features &
                               percent.mt < input$max_mt_percent)
      
      
      if (ncol(seu_filtered) == 0) {
        showModal(modalDialog(
          title = "Error",
          "No cells remain after quality filtering. Please adjust your QC parameters.",
          footer = modalButton("Close")
        ))
        seurat_obj(NULL) 
        return(NULL) 
      }
      
      incProgress(0.1, message = "Normalizing data...")
      seu <- NormalizeData(seu_filtered)
      
      incProgress(0.2, message = "Finding Variable Features...")
      seu <- FindVariableFeatures(seu)
      
      incProgress(0.3, message = "Scaling Data...")
      seu <- ScaleData(seu, features = VariableFeatures(seu))
      seu <- RunPCA(seu)
      
      incProgress(0.4, message = "Applying Parameters")
      dims_use <- 1:input$dims_pca
      resolution_use <- input$resolution
      
      seu <- FindNeighbors(seu, dims = dims_use)
      seu <- FindClusters(seu, resolution = resolution_use)
      
      incProgress(0.5, message = "Running UMAP")
      seu <- RunUMAP(seu, dims = dims_use)
      
      incProgress(0.6, message = "Depicting UMAP")
      seurat_obj(seu)
    })
  })
  
  
  observeEvent(seurat_obj(), { # seurat_obj の値が変化したときにのみ実行
    req(seurat_obj()) # Seuratオブジェクトが存在することを要求
    
    current_obj <- seurat_obj()
    
    # Assays(current_obj) の結果を安全に取得
    # エラー発生時のフォールバックも考慮
    assay_choices <- tryCatch(
      {
        Seurat::Assays(current_obj) # Seurat:: を明示的に追加
      },
      error = function(e) {
        warning(paste("Error getting assays from Seurat object:", e$message))
        return(c("RNA")) # エラー時はデフォルトで "RNA" を返す
      }
    )
    
    # assay_choices が確実に文字型ベクトルであることを保証
    if (!is.character(assay_choices) || length(assay_choices) == 0) {
      assay_choices <- c("RNA") # 不正な場合は "RNA" をデフォルトとする
    }
    
    choices_to_use <- c("RNA") # 初期選択肢は "RNA" のみ
    selected_to_use <- "RNA" # 初期選択も "RNA"
    
    # input$selected_assay の現在の値を安全に取得
    current_input_assay <- isolate(input$selected_assay)
    
    # current_input_assay が有効な値（NULLでなく、文字型、長さ1）であれば、それを優先
    if (!is.null(current_input_assay) && is.character(current_input_assay) && length(current_input_assay) == 1) {
      if (current_input_assay %in% assay_choices) { # ここで %in% が評価される
        selected_to_use <- current_input_assay
      }
    }
    
    # "tfsulm" アッセイが存在する場合、選択肢に追加
    if ("tfsulm" %in% assay_choices) {
      choices_to_use <- c("RNA", "tfsulm")
      
      # もし現在選択されている値が新しい選択肢に含まれていない場合、"RNA" に戻す
      if (!(selected_to_use %in% choices_to_use)) {
        selected_to_use <- "RNA"
      }
    }
    
    # selectInput を更新
    updateSelectInput(session, "selected_assay",
                      choices = choices_to_use,
                      selected = selected_to_use)
  })
  
  
  observe({
    
    req(seurat_obj(), input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    genes <- rownames(temp_obj)
    
    updateSelectizeInput(session, 'selected_gene', choices = genes, server = TRUE)
  })
  
  
  observeEvent(input$calc_marker, {
    req(seurat_obj())
    showModal(modalDialog("Calculating marker genes, please wait...", footer = NULL))
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- "RNA"
    markers <- FindAllMarkers(temp_obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.4, test.use = "MAST")
    marker_result(markers)
    removeModal()
  })
  
  output$cluster_selector <- renderUI({
    req(marker_result())
    clusters <- sort(unique(marker_result()$cluster))
    selectInput("selected_cluster", "Select Cluster", choices = clusters, selected = clusters[1])
  })
  
  output$marker_table <- DT::renderDataTable({
    req(marker_result())
    df <- marker_result()
    datatable(df, filter = "top", options = list(pageLength = 10))
  })
  
  output$umapPlot <- renderPlot({
    req(seurat_obj())
    DimPlot(seurat_obj(), reduction = "umap", label = TRUE,
            label.size = 8, pt.size = input$dotsize) +
      ggtitle("UMAP (clusters)") +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 20)) +
      guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
  })
  
  
  observeEvent(input$preview_download_UMAP, {
    req(seurat_obj())
    showModal(modalDialog(
      title = "Download UMAP Plot",
      fluidRow(
        column(9,
               h4("Preview"),
               plotOutput("umap_preview_plot", height = "auto")
        ),
        column(3,
               h4("Settings"),
               numericInput("umap_download_width", "Width (pixel)", value = 500, min = 1),
               numericInput("umap_download_height", "Height (pixel)", value = 500, min = 1),
               textInput("umap_download_filename", "File Name", value = "UMAP_plot.png"),
               downloadButton("do_download_UMAP", "Download Plot")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$umap_preview_plot <- renderPlot({
    req(seurat_obj())
    plt <- DimPlot(seurat_obj(), reduction = "umap", label = TRUE,
                   label.size = 8, pt.size = input$dotsize) +
      ggtitle("UMAP (clusters)") +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 20)) +
      guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
    
    plt
  }, width = function() {
    req(input$umap_download_width)
    return(input$umap_download_width)
  }, height = function() {
    req(input$umap_download_height)
    return(input$umap_download_height)
  })
  
  
  output$do_download_UMAP <- downloadHandler(
    filename = function() {
      input$umap_download_filename
    },
    content = function(file) {
      req(seurat_obj())
      plt <- DimPlot(seurat_obj(), reduction = "umap", label = TRUE,
                     label.size = 8, pt.size = input$dotsize) +
        ggtitle("UMAP (clusters)") +
        theme(axis.text.y = element_text(size = 15, color = "black"),
              axis.text.x = element_text(size = 15, color = "black"),
              axis.line.y = element_line(linewidth = 1),
              axis.line.x = element_line(linewidth = 1),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20),
              axis.ticks.length.y = unit(.0, "cm"),
              axis.ticks.length.x = unit(.0, "cm"),
              legend.text = element_text(size = 20)) +
        guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
      
      ggsave(file, plot = plt,
             width = input$umap_download_width,
             height = input$umap_download_height,
             units = "px",
             dpi = 96)
    }
  )
  
  output$featurePlot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    FeaturePlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 20),
            title = element_text(size = 20, color = "black"),
            legend.key.height = unit(20, "pt")) &
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  })
  
  
  observeEvent(input$preview_download_feature, {
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    showModal(modalDialog(
      title = paste0("Download Feature Plot (", input$selected_gene, ")"),
      fluidRow(
        column(9,
               h4("Preview"),
               plotOutput("feature_preview_plot", height = "auto")
        ),
        column(3,
               h4("Settings"),
               numericInput("feature_download_width", "Width (pixel)", value = 500, min = 1),
               numericInput("feature_download_height", "Height (pixel)", value = 500, min = 1),
               textInput("feature_download_filename", "File Name", value = paste0("FeaturePlot_", input$selected_gene, ".png")),
               downloadButton("do_download_feature", "Download Plot")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$feature_preview_plot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    plt <- FeaturePlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 20),
            title = element_text(size = 20, color = "black"),
            legend.key.height = unit(20, "pt"))
    plt
  }, width = function() {
    req(input$feature_download_width)
    return(input$feature_download_width)
  }, height = function() {
    req(input$feature_download_height)
    return(input$feature_download_height)
  })
  
  
  output$do_download_feature <- downloadHandler(
    filename = function() {
      input$feature_download_filename
    },
    content = function(file) {
      req(seurat_obj(), input$selected_gene, input$selected_assay)
      temp_obj <- seurat_obj()
      DefaultAssay(temp_obj) <- input$selected_assay
      plt <- FeaturePlot(temp_obj, features = input$selected_gene) +
        ggtitle(input$selected_gene) +
        theme(axis.text.y = element_text(size = 15, color = "black"),
              axis.text.x = element_text(size = 15, color = "black"),
              axis.line.y = element_line(linewidth = 1),
              axis.line.x = element_line(linewidth = 1),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20),
              axis.ticks.length.y = unit(.0, "cm"),
              axis.ticks.length.x = unit(.0, "cm"),
              legend.text = element_text(size = 20),
              title = element_text(size = 20, color = "black"))
      ggsave(file, plot = plt,
             width = input$feature_download_width,
             height = input$feature_download_height,
             units = "px",
             dpi = 96)
    }
  )
  
  output$dotPlot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    DotPlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 20, color = "black", angle = 0),
            axis.text.x = element_text(size = 20, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 15),
            title = element_text(size = 20, color = "black"))
  })
  
  
  observeEvent(input$preview_download_dotplot, {
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    showModal(modalDialog(
      title = paste0("Download Dot Plot (", input$selected_gene, ")"),
      fluidRow(
        column(9,
               h4("Preview"),
               plotOutput("dot_preview_plot", height = "auto")
        ),
        column(3,
               h4("Settings"),
               numericInput("dot_download_width", "Width (pixel)", value = 300, min = 1),
               numericInput("dot_download_height", "Height (pixel)", value = 500, min = 1),
               textInput("dot_download_filename", "File Name", value = paste0("DotPlot_", input$selected_gene, ".png")),
               downloadButton("do_download_dotplot", "Download Plot")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$dot_preview_plot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    plt <- DotPlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 20, color = "black", angle = 0),
            axis.text.x = element_text(size = 20, color = "black"),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            legend.text = element_text(size = 15),
            title = element_text(size = 20, color = "black"))
    plt
  }, width = function() {
    req(input$dot_download_width)
    return(input$dot_download_width)
  }, height = function() {
    req(input$dot_download_height)
    return(input$dot_download_height)
  })
  
  
  output$do_download_dotplot <- downloadHandler(
    filename = function() {
      input$dot_download_filename
    },
    content = function(file) {
      req(seurat_obj(), input$selected_gene, input$selected_assay)
      temp_obj <- seurat_obj()
      DefaultAssay(temp_obj) <- input$selected_assay
      plt <- DotPlot(temp_obj, features = input$selected_gene) +
        ggtitle(input$selected_gene) +
        theme(axis.text.y = element_text(size = 20, color = "black", angle = 0),
              axis.text.x = element_text(size = 20, color = "black"),
              axis.line.y = element_line(linewidth = 1),
              axis.line.x = element_line(linewidth = 1),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20),
              axis.ticks.length.y = unit(.0, "cm"),
              axis.ticks.length.x = unit(.0, "cm"),
              legend.text = element_text(size = 15),
              title = element_text(size = 20, color = "black"))
      ggsave(file, plot = plt,
             width = input$dot_download_width,
             height = input$dot_download_height,
             units = "px",
             dpi = 96)
    }
  )
  
  output$vlnPlot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    VlnPlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 25, color = "black", angle = 0),
            axis.text.x = element_text(size = 25, color = "black", angle = 0),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 25, color = "black"),
            axis.title.x = element_text(size = 25, color = "black"),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            axis.ticks.y = element_line(linewidth = 1),
            axis.ticks.x = element_line(linewidth = 1),
            title = element_text(size = 20, color = "black")) +
      NoLegend()
  })
  
  
  observeEvent(input$preview_download_vlnplot, {
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    showModal(modalDialog(
      title = paste0("Download Violin Plot (", input$selected_gene, ")"),
      fluidRow(
        column(9,
               h4("Preview"),
               plotOutput("vln_preview_plot", height = "auto")
        ),
        column(3,
               h4("Settings"),
               numericInput("vln_download_width", "Width (pixel)", value = 1000, min = 1),
               numericInput("vln_download_height", "Height (pixel)", value = 500, min = 1),
               textInput("vln_download_filename", "File Name", value = paste0("VlnPlot_", input$selected_gene, ".png")),
               downloadButton("do_download_vlnplot", "Download Plot")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$vln_preview_plot <- renderPlot({
    req(seurat_obj(), input$selected_gene, input$selected_assay)
    temp_obj <- seurat_obj()
    DefaultAssay(temp_obj) <- input$selected_assay
    plt <- VlnPlot(temp_obj, features = input$selected_gene) +
      ggtitle(input$selected_gene) +
      theme(axis.text.y = element_text(size = 25, color = "black", angle = 0),
            axis.text.x = element_text(size = 25, color = "black", angle = 0),
            axis.line.y = element_line(linewidth = 1),
            axis.line.x = element_line(linewidth = 1),
            axis.title.y = element_text(size = 25, color = "black"),
            axis.title.x = element_text(size = 25, color = "black"),
            axis.ticks.length.y = unit(.0, "cm"),
            axis.ticks.length.x = unit(.0, "cm"),
            axis.ticks.y = element_line(linewidth = 1),
            axis.ticks.x = element_line(linewidth = 1),
            title = element_text(size = 20, color = "black")) +
      NoLegend()
    plt
  }, width = function() {
    req(input$vln_download_width)
    return(input$vln_download_width)
  }, height = function() {
    req(input$vln_download_height)
    return(input$vln_download_height)
  })
  
  
  output$do_download_vlnplot <- downloadHandler(
    filename = function() {
      input$vln_download_filename
    },
    content = function(file) {
      req(seurat_obj(), input$selected_gene, input$selected_assay)
      temp_obj <- seurat_obj()
      DefaultAssay(temp_obj) <- input$selected_assay
      plt <- VlnPlot(temp_obj, features = input$selected_gene) +
        ggtitle(input$selected_gene) +
        theme(axis.text.y = element_text(size = 25, color = "black", angle = 0),
              axis.text.x = element_text(size = 25, color = "black", angle = 0),
              axis.line.y = element_line(linewidth = 1),
              axis.line.x = element_line(linewidth = 1),
              axis.title.y = element_text(size = 25, color = "black"),
              axis.title.x = element_text(size = 25, color = "black"),
              axis.ticks.length.y = unit(.0, "cm"),
              axis.ticks.length.x = unit(.0, "cm"),
              axis.ticks.y = element_line(linewidth = 1),
              axis.ticks.x = element_line(linewidth = 1),
              title = element_text(size = 20, color = "black")) +
        NoLegend()
      ggsave(file, plot = plt,
             width = input$vln_download_width,
             height = input$vln_download_height,
             units = "px",
             dpi = 96)
    }
  )
  
  output$download_marker <- downloadHandler(
    filename = function() {
      paste0("markers_cluster", ".csv")
    },
    content = function(file) {
      df <- marker_result()
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$draw_heatmap, {
    req(seurat_obj(), marker_result())
    withProgress(message = "Generating heatmap...", {
      top_genes <- marker_result() %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = input$top_n_genes, with_ties = FALSE) %>%
        pull(gene) %>%
        unique()
      
      obj <- seurat_obj()
      DefaultAssay(obj) <- "RNA" 
      obj <- ScaleData(obj, features = top_genes)
      
      plt <- DoHeatmap(obj, features = top_genes,
                       group.colors = NULL)+
        theme(axis.text.y = element_text(size = 20, colour = "black")) +
        NoLegend()
      heatmap_obj(plt)
    })
  })
  
  output$heatmapPlot <- renderPlot({
    req(heatmap_obj())
    heatmap_obj()
  })
  
  
  observeEvent(input$preview_download_heatmap, {
    req(heatmap_obj())
    showModal(modalDialog(
      title = "Download Heatmap",
      fluidRow(
        column(10,
               h4("Preview"),
               plotOutput("Heatmap_preview_plot", height = "auto")
        ),
        column(2,
               h4("Settings"),
               numericInput("Heatmap_download_width", "Width (pixel)", value = 600, min = 10),
               numericInput("Heatmap_download_height", "height (pixel)", value = 1000, min = 10),
               textInput("Heatmap_download_filename", "File name", value = "Heatmap_plot.png"),
               downloadButton("do_download_Heatmap", "Download Heatmap")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$Heatmap_preview_plot <- renderPlot({
    req(heatmap_obj())
    heatmap_obj()
  }, width = function() {
    req(input$Heatmap_download_width)
    return(input$Heatmap_download_width)
  }, height = function() {
    req(input$Heatmap_download_height)
    return(input$Heatmap_download_height)
  })
  
  output$do_download_Heatmap <- downloadHandler(
    filename = function() {
      input$Heatmap_download_filename
    },
    content = function(file) {
      req(heatmap_obj())
      plt <- heatmap_obj()
      
      ggsave(file, plot = plt,
             width = input$Heatmap_download_width,
             height = input$Heatmap_download_height,
             units = "px",
             dpi = 96)
    }
  )
  
  
  observeEvent(input$draw_tf_heatmap, {
    req(seurat_obj(), input$target_spiecies)
    
    withProgress(message = "Running TF activity prediction...", value = 0, {
      
      incProgress(0.1, message = "Loading TF network...")
      net <- collectri_net()
      
      
      incProgress(0.2, message = "Extracting expression data...")
      
      mat <- as.matrix(GetAssayData(seurat_obj(), assay = "RNA", slot = "data"))
      
      
      incProgress(0.5, message = "Inferring TF activities with decoupleR...")
      acts <- run_ulm(mat = mat, net = net, .source = 'source', .target = 'target',
                      .mor = 'mor', minsize = 5) # verbose = FALSE 引数を削除
      
      
      incProgress(0.8, message = "Storing TF activities in Seurat object...")
      tf_assay_data <- acts %>%
        pivot_wider(id_cols = 'source', names_from = 'condition',
                    values_from = 'score') %>%
        column_to_rownames('source') %>%
        as.matrix()
      
      temp_seu <- seurat_obj()
      temp_seu[['tfsulm']] <- CreateAssayObject(counts = tf_assay_data)
      
      DefaultAssay(temp_seu) <- "tfsulm"
      temp_seu <- Seurat::ScaleData(temp_seu, features = rownames(tf_assay_data), verbose = FALSE)
      temp_seu@assays$tfsulm@data <- temp_seu@assays$tfsulm@scale.data
      
      seurat_obj(temp_seu) 
      tf_activities(acts) 
      
      
      df_long <- t(as.matrix(temp_seu@assays$tfsulm@data)) %>%
        as.data.frame() %>%
        mutate(cluster = Idents(temp_seu)) %>%
        pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
        group_by(cluster, source) %>%
        summarise(mean = mean(score), .groups = 'drop')
      
      tf_activities_df(df_long) 
      
      incProgress(1, message = "TF Analysis Complete and Heatmap Generating...")
    })
  })
  
  
  output$tfActivityHeatmap <- renderPlot({
    req(tf_activities_df(), input$top_n_active_TF)
    
    n_tfs_val <- input$top_n_active_TF
    df_long <- tf_activities_df()
    
    tfs <- df_long %>%
      group_by(source) %>%
      summarise(std = sd(mean), .groups = 'drop') %>%
      arrange(-abs(std)) %>%
      head(n_tfs_val) %>%
      pull(source)
    
    top_acts_mat_local <- df_long %>%
      filter(source %in% tfs) %>%
      pivot_wider(id_cols = 'cluster', names_from = 'source',
                  values_from = 'mean') %>%
      column_to_rownames('cluster') %>%
      as.matrix()
    
    if (!is.null(seurat_obj())) {
      ordered_clusters <- sort(as.numeric(as.character(rownames(top_acts_mat_local))))
      top_acts_mat_local <- top_acts_mat_local[as.character(ordered_clusters), , drop = FALSE]
    }
    
    tf_pheatmap_mat(top_acts_mat_local)
    
    palette_length <- 100
    my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)
    max_abs_score <- max(abs(top_acts_mat_local), na.rm = TRUE)
    my_breaks <- c(seq(-max_abs_score, 0, length.out = ceiling(palette_length/2) + 1),
                   seq(0.001, max_abs_score, length.out = floor(palette_length/2)))
    
    pheatmap(top_acts_mat_local, border_color = NA,
             color = my_color,
             breaks = my_breaks,
             angle_col = 45,
             cluster_cols = TRUE,
             cluster_rows = FALSE,
             cellwidth = 18,
             cellheight = 18,
             fontsize = 15,
             fontsize_row = 15,
             main = paste0("TF Activity Heatmap (Top ", n_tfs_val, " TFs)")
    )
  })
  
  
  observeEvent(input$preview_download_tf, {
    req(tf_pheatmap_mat())
    showModal(modalDialog(
      title = "Download TF Activity Heatmap",
      fluidRow(
        column(9,
               h4("Preview"),
               plotOutput("tf_preview_heatmap_plot", height = "auto")
        ),
        column(3,
               h4("Settings"),
               numericInput("tf_download_width_px", "Width (pixels)", value = 600, min = 100),
               numericInput("tf_download_height_px", "Height (pixels)", value = 600, min = 100),
               textInput("tf_download_filename", "File name", value = "TF_Activity_Heatmap.png"),
               downloadButton("do_download_tf_heatmap", "Download Heatmap")
        )
      ),
      size = "l",
      footer = modalButton("Close")
    ))
  })
  
  
  output$tf_preview_heatmap_plot <- renderPlot({
    req(tf_pheatmap_mat())
    
    mat_to_plot_preview <- tf_pheatmap_mat()
    
    palette_length <- 100
    my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)
    max_abs_score <- max(abs(mat_to_plot_preview), na.rm = TRUE)
    my_breaks <- c(seq(-max_abs_score, 0, length.out = ceiling(palette_length/2) + 1),
                   seq(0.001, max_abs_score, length.out = floor(palette_length/2)))
    
    pheatmap(mat_to_plot_preview, border_color = NA,
             color = my_color,
             breaks = my_breaks,
             angle_col = 45,
             cluster_cols = TRUE,
             cluster_rows = FALSE,
             cellwidth = 18,
             cellheight = 18,
             fontsize = 15,
             fontsize_row = 15,
             main = paste0("TF Activity Heatmap (Top ", input$top_n_active_TF, " TFs)")
    )
  }, width = function() {
    req(input$tf_download_width_px)
    return(min(input$tf_download_width_px, 800))
  }, height = function() {
    req(input$tf_download_height_px)
    return(min(input$tf_download_height_px, 800))
  })
  
  
  output$do_download_tf_heatmap <- downloadHandler(
    filename = function() {
      input$tf_download_filename
    },
    content = function(file) {
      req(tf_pheatmap_mat())
      
      mat_to_plot_download <- tf_pheatmap_mat()
      
      palette_length <- 100
      my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)
      max_abs_score <- max(abs(mat_to_plot_download), na.rm = TRUE)
      my_breaks <- c(seq(-max_abs_score, 0, length.out = ceiling(palette_length/2) + 1),
                     seq(0.001, max_abs_score, length.out = floor(palette_length/2)))
      
      png(file, width = input$tf_download_width_px, height = input$tf_download_height_px, units = "px", res = 96)
      pheatmap(mat_to_plot_download, border_color = NA,
               color = my_color,
               breaks = my_breaks,
               angle_col = 45,
               cluster_cols = TRUE,
               cluster_rows = FALSE,
               cellwidth = 18,
               cellheight = 18,
               fontsize = 15,
               fontsize_row = 15,
               main = paste0("TF Activity Heatmap (Top ", input$top_n_active_TF, " TFs)")
      )
      dev.off()
    }
  )
  
  output$tf_activity_table <- DT::renderDataTable({
    req(tf_activities_df())
    df <- tf_activities_df() %>%
      dplyr::mutate(mean = round(mean, 4))
    
    datatable(df, 
              filter = "top",
              options = list(pageLength = 10,
                             scrollX = TRUE))
  })
  
  output$download_tf_table <- downloadHandler(
    filename = function() {
      paste0("tf_activity_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(tf_activities_df())
      write.csv(tf_activities_df(), file, row.names = FALSE)
    }
  )
})