library(shiny)
library(shinyjs)

shinyUI(fluidPage(
  titlePanel("Single-cell RNA-seq Viewer with Seurat"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 actionButton("create_and_load_demo_data", "Create & Load Demo Data"),
                 helpText("Click this button to create and load demo data"),
                 hr(),
                 
                 fileInput("zipfile", "Upload 10X folder (as zip file)", accept = ".zip"),
                 
                 useShinyjs(), 
                 actionButton("run_initial_processing", "Quality check"),
                 hr(),
                 
                 h4("Quality Check Parameters"),
                 numericInput("min_features", "Min Features (nFeature_RNA)", value = 30, min = 0),
                 numericInput("max_features", "Max Features (nFeature_RNA)", value = 5000, min = 1),
                 numericInput("max_mt_percent", "Max Mitocondrial gene %", value = 5, min = 0, max = 100),
                 
                 numericInput("dims_pca", "Number of PCA dimensions (FindNeighbors/RunUMAP)",
                              value = 10, min = 1),
                 numericInput("resolution", "Resolution for clustering (FindClusters)",
                              value = 0.5, min = 0.1, step = 0.1),

                 actionButton("run", "Run Seurat Analysis"),

                 hr(),
                 numericInput("dotsize", "Dot size for UMAP",
                              value = 0.5, min = 0.1, step = 0.1),
                 
                 hr(),
                 # ここにアッセイ選択のUIを追加
                 h4("Assay Selection"),
                 selectInput("selected_assay", "Select Assay for Plots",
                             choices = c("RNA"), 
                             selected = "RNA"),
                 hr(),

                 selectizeInput(
                   "selected_gene", "Select Gene / TF",
                   choices = NULL, 
                   options = list(
                     placeholder = 'Search for a gene or TF...',
                     onInitialize = I('function() { this.setValue(""); }') 
                   )
                 ),
                 hr(),
                 
                 actionButton("calc_marker", "Calculate Marker Genes")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Plots",
                 
                 h4("Quality Check Plots"),
                 fluidRow(
                   column(12, plotOutput("qc_plot_mt", height = "400px"))
                 ),
                 hr(),
                 
                 h4("UMAP (Cluster)"),
                 fluidRow(
                   column(6,
                          actionButton("preview_download_UMAP", "Preview & Download PNG"),
                          plotOutput("umapPlot", height = "450px", width = "500px")),
                   
                   column(6,
                          actionButton("preview_download_feature", "Preview & Download PNG"),
                          plotOutput("featurePlot", height = "450px", width = "500px"))
                 ),
                 hr(),
                 
                 h4("Gene Expression Visualizations"),
                 fluidRow(
                   column(9,
                          actionButton("preview_download_vlnplot", "Preview & Download PNG"),
                          plotOutput("vlnPlot", height = "400px", width = "800px"))
                 ),
                 
                 fluidRow(
                   column(3,
                          actionButton("preview_download_dotplot", "Preview & Download PNG"),
                          plotOutput("dotPlot", height = "450px", width = "350px"))
                 ),
        ),
        
        tabPanel("DEG Heatmap",
                 h4("Heatmap of each cluster marker"),
                 numericInput("top_n_genes", "Top N marker genes per cluster", value = 5, min = 1),
                 actionButton("draw_heatmap", "Generate Heatmap"),
                 br(), br(),
                 actionButton("preview_download_heatmap", "Preview & Download PNG"),
                 plotOutput("heatmapPlot", height = "800px")
        ),
        
        tabPanel("DEG Table",
                 h4("Marker Genes for Each Cluster"),
                 downloadButton("download_marker", "Download Diffrentially Expressed Gene (CSV)"),
                 DT::dataTableOutput("marker_table")
        ),

        tabPanel("TF",
                 h4("Estimation of Transcription Factor Activity"),
                 selectInput("target_spiecies", "Species",
                             choices = c("human" = "human",
                                         "mouse" = "mouse")),
                 numericInput("top_n_active_TF", "Top N most active TFs for Heatmap", value = 20, min = 1),
                 actionButton("draw_tf_heatmap", "Run TF Activity Prediction"),
                 br(), br(),
                 actionButton("preview_download_tf", "Preview & Download PNG"),
                 plotOutput("tfActivityHeatmap", height = "800px", width = "100%") 
        ),
        
        tabPanel("TF table",
                 h4("Transcription Factor Activity Results"),
                 downloadButton("download_tf_table", "Download TF Activity Table (CSV)"),
                 DT::dataTableOutput("tf_activity_table")
        ),
        
        tabPanel("Pathway Activity",
                 h4("Pathway Activity Analysis with Progeny"),
                 selectInput("target_species_progeny", "Species",
                             choices = c("Human" = "Human",
                                         "Mouse" = "Mouse")),
                 helpText("Estimate the activity of major signaling pathway based on gene expression."),
                 hr(),
                 
                 actionButton("run_progeny_analysis", "Run Pathway Activity Analysis"),
                 br(), br(),
                 
                 h5("Pathway Activity Heatmap"),
                 column(12, align = "center",
                        plotOutput("progeny_heatmap", height = "600px", width = "100%")
                 ),
                 downloadButton("download_progeny_heatmap", "Download Heatmap (PNG)"),
                 
                 hr(),
                 h5("Pathway Activity Table"),
                 downloadButton("download_progeny_table", "Download Table (CSV)"),
                 DT::dataTableOutput("progeny_table")
                 )
      )
    )
  )
))