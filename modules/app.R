library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(bslib)
library(visNetwork)
library(shinyFeedback)
library(DT)
library(colourpicker)
library(ggplot2)
library(ggpubr)
library(vegan)
library(tidyr)
library(ape) # pcoa
library(uwot) # umap
library(Rtsne) # tsne
library(igraph) # network analysis (cluster)
library(RColorBrewer)
library(patchwork)
library(PreLectR)


# Rscript -e "shiny::runApp('${app_dir}', launch.browser=function(url) { utils::browseURL(url, browser='/usr/bin/google-chrome') })"

# env var
working_dir <- Sys.getenv("WORKING_DIR")
script_dir <- Sys.getenv("SCRIPT_DIR")
threads <- Sys.getenv("THREADS")
# working_dir <- "~/Documents/YFC_MCI_project/96mci250normal_bacnex"
# script_dir <- "~/Documents/vscode/files/BacNex/modules"
# threads <- 5

# check working_dir exit
if(!dir.exists(working_dir)) {
  stop("The directory does not exist.")
}


# app onstart
tmp_dir <- file.path(working_dir, "tmp_files")
onstart <- function() {
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }
}


# taxa table
data <- read.table(file.path(working_dir, 'taxa_table.csv'), sep = ',', header = T, row.names = 1, stringsAsFactors = FALSE)
split_result <- strsplit(row.names(data), "\\.")
genus <- sapply(split_result, function(x) sub("g__", "", x[1]))
species <- sapply(split_result, function(x) sub("s__", "", x[2]))
taxa <- data.frame(Genus=genus, Species=species)
rownames(taxa) <- row.names(data)

# 72 color map
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Functions
# network analysis (cluster)
make_net <- function(df) { # make the adjacency matrix for loading graph to igraph
  name = unique(c(as.character(df$from),as.character(df$to)))
  mat <- matrix(0, nrow = length(name), ncol = length(name))
  rownames(mat) <- name
  colnames(mat) <- name
  for(i in c(1:dim(df)[1])){
    anchor1 <- as.character(df$from[i])
    anchor2 <- as.character(df$to[i])
    w <- df$weight[i]
    mat[anchor1, anchor2] <- w
    mat[anchor2, anchor1] <- w
  }
  return(as.matrix(mat))
}


assignModule <- function(nodeT, module, Method, cutoff) { # assign the graph partition for each node
  nodeT$module <- "unclustered"
  for(x in nodeT$id){
    if(x %in% module$names) {
      nodeT[nodeT$id == x, "module"] <- paste0(c("module", membership(module)[x]), collapse = "_")
    } else {
      next
    }
  }
  sortModule <- sort(table(nodeT$module),decreasing = T)   # name the module from the largest group to the smallest one
  count <- 1
  for(x in names(sortModule)){
    if(x == "uncluster"){next}
    nodeT$module[nodeT$module == x] <- paste0(c(Method,count), collapse = "_")
    count <- count + 1
  }
  sortModule <- sort(table(nodeT$module),decreasing = T)   # if the member number of module is less than cutoff, discard it
  nodeT$module[nodeT$module %in% names(sortModule[sortModule < cutoff])] <- "unclustered"
  colnames(nodeT)[ncol(nodeT)] <- Method
  return(nodeT)
}






# ===== UI Main =====
ui <- dashboardPage(
  dashboardHeader(title = "BacNex",
                  tags$li(
                    class = "dropdown",
                    style = "padding: 1px; margin-right: 4px;",
                    tags$span("Beta 2.1", style = "font-size: 10px; color: black;"))
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("home")),
      menuItem("Analysis pipeline", tabName = "Pipeine", icon = icon("chart-simple"),
        menuSubItem("Metadata submission", tabName = "Metadata", icon = icon("angle-right")),
        menuSubItem("Diversity analysis", tabName = "Diversity", icon = icon("angle-right")),
        menuSubItem("Differental taxa", tabName = "PreLect", icon = icon("angle-right")),
        menuSubItem("Network construction", tabName = "Network", icon = icon("angle-right")),
        menuSubItem("Network filter", tabName = "Networkfilter", icon = icon("angle-right")),
        menuSubItem("Pathway analysis", tabName = "PathwayAnalysis", icon = icon("angle-right")),
        menuSubItem("Network analysis", tabName = "NetworkAnalysis", icon = icon("angle-right"))
      )
      # menuItem("Antibiotic genes", tabName = "ARG", icon = icon("dna"))
    )
  ),
  dashboardBody(
    # add_busy_spinner(spin = "atom", position="bottom-right", color="indianred"), # busy spinner
    tabItems(
      # home
      tabItem("Home",
              fluidRow(
                box(title="BacNex", width=12, status="primary", solidHeader=TRUE,
                    h4("Welcome to BacNex"),
                    p(HTML('Watch tutorial from <a href="https://github.com/RynoLiu/BacNex.git" target="_blank">BacNex github</a>')),
                    HTML('<br>'),
                    p("Please select the Analysis pipeline from the sidebar menu.")
                )
              )
      ),

      # module 1
      tabItem("Metadata",
              fluidRow(
                box(title = "Metadata Input Format", width =12, status = "primary", solidHeader = TRUE,
                    tabBox(width = 9,
                           tabPanel("Binary classification",
                                    helpText("Please ensure the column names match the demo table"),
                                    HTML("<ol>\
                                          <li>Accession_ID : The unique identifier for each sample</li>\
                                          <li>Labels : A character vector assigned to each sample</li>\
                                          </ol>"),
                                    div(DT::dataTableOutput("demo_meta1"))),

                           tabPanel("Multi classification",
                                    helpText("Please ensure the column names match the demo table"),
                                    HTML("<ol>\
                                          <li>Accession_ID : The unique identifier for each sample</li>\
                                          <li>Labels : A character vector assigned to each sample</li>\
                                          </ol>"),
                                    div(DT::dataTableOutput("demo_meta2"))),

                           tabPanel("Regression",
                                    helpText("Please ensure the column names match the demo table"),
                                    HTML("<ol>\
                                          <li>Accession_ID : The unique identifier for each sample</li>\
                                          <li>Labels : A numeric vector assigned to each sample</li>\
                                          </ol>"),
                                    div(DT::dataTableOutput("demo_meta3"))),

                           tabPanel("Time to event",
                                    helpText("Please ensure the column names match the demo table"),
                                    HTML("<ol>\
                                          <li>Accession_ID : The unique identifier for each sample</li>\
                                          <li>Event : Numeric vector. Binary indicator for event status: 1 for event occurrence (e.g., death), and 0 for censored (e.g., alive).</li>\
                                          <li>Duration : Numeric vector. Duration of follow-up time for each sample.</li>\
                                          </ol>"),
                                    div(DT::dataTableOutput("demo_meta4")))
                    )
                )
              ),
              fluidRow(
                box(title = "Upload Metadata", width =12, status = "primary", solidHeader = TRUE,
                    fileInput(inputId = "meta_file",label = HTML("Input your metadata with csv format :"), accept = c(".csv"), placeholder = "metadata.csv", multiple = F, width="40%"),
                    radioButtons("object_type", "Object type", c("Binary classification"='BC', "Multi classification"="MC", "Regression"="Rg", "Time to event"="Cox"), selected="BC"),
                    actionButton("submit_meta_button", "Upload"),
                    helpText("Format checking"),
                    textOutput("meta_check1"),
                    textOutput("meta_check2"),
                    uiOutput("meta_check3"),
                    HTML('<br>'),
                    div(DT::dataTableOutput("submit_meta"))
                )
              )
      ),

      # module 2
      tabItem("Diversity",
              card(input_switch("prelect_switch", "Use PreLect features (In the diversity analysis)", value=FALSE, width="700px")),
              fluidRow(
                box(title = "Alpha Diversity", width=12, status="primary", solidHeader = TRUE,
                    fluidRow(
                      column(2,
                        fluidRow(
                          column(10,
                            loadingButton("alpha_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                            HTML('<br><br>'))),
                        fluidRow(column(8, uiOutput("a_color_input"))),
                        fluidRow(column(8, uiOutput("a_size_input"))) ),
                      column(10, plotOutput("alpha_plot", height="500px")),
                    ),
                    HTML('<br><br>'),
                    conditionalPanel(
                      condition = "output.alpha_plot !== null && output.alpha_plot !== undefined",
                      actionButton("alpha_save", "Save Results", icon=icon("download"))
                    )
                ) # box
              ),
              fluidRow(
                box(title = "Beta Diversity", width =12, status = "primary", solidHeader = TRUE,
                    # selectInput("dimension", label="Dimension reduction", choices=c("PCoA"="pcoa", "NMDS"="nmds", "t-SNE"="tsne", "UMAP"="umap"), selected='pcoa', width="15%"),
                    # loadingButton("beta_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    # HTML('<br><br>'),
                    fluidRow(
                      column(2,
                        fluidRow(
                          column(10,
                          selectInput("dimension", label="Dimension reduction", choices=c("PCoA"="pcoa", "NMDS"="nmds", "t-SNE"="tsne", "UMAP"="umap"), selected='pcoa'),
                          loadingButton("beta_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                          HTML('<br><br>') )),
                        fluidRow(column(8, uiOutput("b_color_input"))),
                        fluidRow(column(8, uiOutput("b_size_input"))) ),
                      column(10, plotOutput("beta_plot", height="500px"))
                    ),
                    HTML('<br><br>'),
                    conditionalPanel(
                      condition = "output.beta_plot !== null && output.beta_plot !== undefined",
                      actionButton("beta_save", "Save Results", icon=icon("download"))
                    )
                ) # box
              ),
              fluidRow(
                box(title = "Composition", width =12, status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(2,
                        fluidRow(column(10,
                          selectInput("label_selector_com", label="Cohort", choices=""),
                          selectInput("target_level_com", label="Level", choices=c("Phylum"="phylum", "Class"="class",
                                      "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"), selected="phylum"),
                        loadingButton("composed_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                        HTML('<br><br>'),)),
                        fluidRow(column(8, uiOutput("com_size_input"))) ),
                      column(10, plotOutput("com_plot", height="500px"))
                    ),
                    HTML('<br><br>'),
                    conditionalPanel(
                      condition = "output.com_plot !== null && output.com_plot !== undefined",
                      actionButton("com_save", "Save Results", icon=icon("download")))
                ) # box
              ) # row
      ), # tabItem

      # module 3
      tabItem("PreLect",
              fluidRow(
                box(title = "Lambda Tuning", width=12, status="primary", solidHeader = TRUE,
                    textOutput("model_usage"),
                    HTML('<br>'),
                    fluidRow(
                      column(3, selectInput("norm_method", label="Normalization method", choices=c('Z-norm'='z','Z-norm + min-max scaling'='zm', 'CLR'='clr'), selected='z')),
                      column(3, numericInput('auto_step', label='Number of examine lambda', value=30, min=10, max=50, step=1)),
                      column(3, numericInput('max_iter', label='Maximum iterations for fitting', value=10000, min=1000, max=100000, step=1)),
                      column(3, numericInput('lr', label='Learning rate', value=0.001, min=0.1, max=0.0001))
                    ),
                    loadingButton("tuning_run", "Run lambda tuning", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    ) # box
              ), # row

              fluidRow(
                box(title="Lambda Decision", width=12,status="primary", solidHeader=TRUE,
                    fluidRow(
                      column(3, numericInput('max_depth', label='Maximum depth for dTree', value=3, min=1, max=10, step=1)),
                      column(3, numericInput('min_bucket', label='Minimun point for each sequement', value=3, min=1, max=50, step=1))
                    ),
                    loadingButton("opt_lambda", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML('<br><br>'),
                    plotOutput("lambda_dec_plot"),
                    h4(textOutput("opt_lambda_value"))
                    ) # box
              ), # row

              fluidRow(
                box(title = "Feature Visualization", width=12, status="primary", solidHeader = TRUE,
                    # selectInput("PL_task", label="PreLect task", choices=c('classification'='c', 'regression'='r'), selected='c'),
                    numericInput('opt_lmbd', label='optimized lambda', value=1e-5, min=0, max=10, step=1, width="25%"),
                    loadingButton("run_prelect", "Run PreLect", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML('<br><br>'),
                    div(DT::dataTableOutput("PLres_table")),
                    HTML('<br>'),
                    textOutput("PLres_note1"),
                    textOutput("PLres_note2"),
                    HTML('<br>'),
                    conditionalPanel(
                      condition = "output.PLres_table !== null && output.PLres_table !== undefined",
                      actionButton("PLres_save", "Save Results", icon=icon("download")))
                    ) # box
              ) # row
      ), # tabItem

      # module 4
      tabItem("Network",
               fluidRow(
                 box(title="Network construction", width=12, status="primary", solidHeader = TRUE,
                    fluidRow(
                      column(3, selectInput("label_selector", label="Cohort", choices="")),
                      column(3, numericInput('adj_p', label='p-value to decide taxa connections', value=0.05, min=0, max=0.05, step=0.01))
                    ),
                    loadingButton("make_network_button", "Generate network", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML('<br><br>'),
                    h4(textOutput("graph_title")),
                    div(DT::dataTableOutput("edge_table")),
                    conditionalPanel(
                      condition = "output.edge_table !== null && output.edge_table !== undefined",
                      actionButton("network_save", "Save Results", icon=icon("download")))
                    ) # box
              ) # row
      ),

      # module 5
      tabItem("Networkfilter",
              fluidRow(
                box(title="Network filter", width=12, status="primary", solidHeader = TRUE,
                    fluidRow(
                      column(3, selectInput("label_selector_f", label="Cohort", choices="")),
                      column(3, selectInput("node_size", label="Node size", choices=c('Fold-change'='FC', 'Degree'='degree', 'KO counts'='KO'), selected='FC')),
                      column(3, selectInput("select_rank", label="colored by group", choices=c('Phylum'='phylum', 'Class'='class',
                                'Order'='order', "Family"='family', 'Genus'='genus', 'No color'='none'), selected='phylum'))
                    ),
                    fluidRow(
                      column(3, sliderInput("size_range", "Range of node size", min = 0.1, max = 100, value = c(10,40), step=1)),
                      column(3, sliderInput("width_range", "Range of edge width", min = 0.1, max = 100, value = c(0.1,10), step=1)),
                    ),
                    fluidRow(
                      column(3, numericInput('PL_w', label='PreLect model : Weight', value=0.5, min=0, max=100, step=0.1)),
                      column(3, numericInput('edge_w', label='Edge weight : Odds ratio', value=1, min=0, max=100, step=1))
                    ),
                    loadingButton("filter_network_button", "Generate network", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML('<br><br>'),
                    h4(textOutput("graph_f_title")),
                    textOutput("filter_note"),
                    fluidRow(
                      column(6, uiOutput("edge_f_ui")),
                      column(6, uiOutput("node_f_ui"))
                    ),
                    HTML('<br>'),
                    visNetworkOutput("network_graph", width = "100%", height = "800px"),
                    conditionalPanel(
                      condition = "output.network_graph !== null && output.network_graph !== undefined",
                      actionButton("network_f_save", "Save Results", icon=icon("download")))
                  ) # box
            ) # row
      ),

      # module 6
      tabItem("PathwayAnalysis",
              fluidRow(
                box(
                  title="Pathway analysis", width=12, status="primary", solidHeader = TRUE,
                  fluidRow(
                    column(3, selectInput("label_selector_p", label="Cohort", choices="" )),
                    column(3, numericInput('pathway_p', label='Adjusted p_value threshold', value=0.05, min=0, max=100, step=0.01))
                  ),
                  loadingButton("pathway_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                  HTML('<br><br>'),
                  h4(textOutput("pathway_title")),
                  textOutput("pathway_note"),
                  HTML('<br>'),
                  div(DT::dataTableOutput("vital_pathway")),
                  HTML('<br>'),
                  fluidRow(
                    column(2, uiOutput("path_bar_size_input")),
                    column(10, plotOutput("path_plot", inline=TRUE), br())
                  ),
                  uiOutput("bar_note"),
                  HTML('<br>'),
                  conditionalPanel(
                    condition = "output.vital_pathway !== null && output.vital_pathway !== undefined",
                    actionButton("pathway_save", "Save Results", icon=icon("download")))
                ), # box

                box(
                  title="Pathway sub-network visualization", width=12, status="primary", solidHeader=TRUE,
                  fluidRow(
                    column(3, selectInput("label_selector_path_net", label="Cohort", choices="")),
                    column(3, tagList(
                              textInput("map_id", label="Pathway ID", value=""),
                              div(style = "font-size: 10px; color: gray; margin-top: -15px;",
                                  "The kegg map id. e.g., map00010"))),
                  ),
                  loadingButton("pathway_net_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                  HTML('<br><br>'),
                  uiOutput("pathway_net_title"),
                  HTML('<br>'),
                  fluidRow(
                    column(6, uiOutput("pathway_net_edgeUI")),
                    column(6, uiOutput("pathway_net_nodeUI"))
                  ),
                  HTML('<br>'),
                  visNetworkOutput("pathway_net_graph", width = "100%", height = "800px"),
                  conditionalPanel(
                    condition = "output.pathway_net_graph !== null && output.pathway_net_graph !== undefined",
                    actionButton("pathway_net_save", "Save Results", icon=icon("download")))
                ) # box
            ) # row
      ),

      # module 7
      tabItem("NetworkAnalysis",
              fluidRow(
                box(
                  title="Network cluster", width=12, status="primary", solidHeader = TRUE,
                  fluidRow(
                    column(3, selectInput("label_selector_netA", label="Cohort", choices="" )),
                    column(3,
                           tagList(
                                numericInput('clust_cutoff', label='Threshold', value=3, min=0, max=100, step=1),
                                div(style = "font-size: 10px; color: gray; margin-top: -15px;",
                                    "If the number of members in cluster less than threshold, the cluster will be discarded.")
                            ))
                  ),
                  HTML("<br>"),
                  column(2, input_switch("B_switch", "Edge betweenness", value = TRUE)),
                  column(2, input_switch("W_switch", "Random walks", value = TRUE)),
                  column(2, input_switch("L_switch", "Louvain", value = TRUE)),
                  column(2, input_switch("G_switch", "Greedy", value = TRUE)),
                  column(4, input_switch("LE_switch", "Non-negative eigenvector", value = TRUE)),
                  loadingButton("cluster_out", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                  HTML('<br><br>'),
                  h4(uiOutput("modularity_title")),
                  div(style = "width: 30%;", DT::dataTableOutput("modularity_table")),
                  HTML('<br>'),
                  h4(uiOutput("clustNode_title")),
                  div(DT::dataTableOutput("clustNode_table")),
                  HTML('<br>'),
                  h4(uiOutput("clust_stat_title")),
                  plotOutput("clust_stat_plot", height = "800px"),
                  conditionalPanel(
                    condition = "output.clustNode_table !== null && output.clustNode_table !== undefined",
                    actionButton("clust_stat_save", "Save Results", icon=icon("download")))
                ), # box

                box(
                  title="Cluster visualization", width=12, status="primary", solidHeader = TRUE,
                  selectInput("clust_select", label="Method to visualize",
                              choices=c("Edge betweenness"="Betweenness", "Random walks"="RandomWalks",
                              "Louvain"="Louvain","Greedy" ="Greedy", "Non-negative eigenvector"="Leading_eigen"), selected="Betweenness", width="20%"),
                  loadingButton("make_clust_network", "Generate network", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                  HTML('<br><br>'),
                  uiOutput("forest_plot_ui"),
                  HTML('<br>'),
                  uiOutput("forest_df_ui"),
                  HTML('<br><br>'),
                  fluidRow(
                    column(6, uiOutput("edge_clust_ui")),
                    column(6, uiOutput("node_clust_ui"))
                  ),
                  HTML('<br><br>'),
                  visNetworkOutput("clust_graph", width = "100%", height = "800px"),
                  conditionalPanel(
                    condition = "output.clust_graph !== null && output.clust_graph !== undefined",
                    actionButton("clust_graph_save", "Save Results", icon=icon("download")))
                ), # box

              box(
                title="Clustered pathway analysis", width=12, status="primary", solidHeader = TRUE,
                fluidRow(
                  column(3, selectInput("label_netA_path", label="Cohort", choices="" )),
                  column(3, selectInput("netA_path_clust_id", label="Clustered group ID", choices="")),
                  column(3, numericInput('netA_path_p', label='Adjusted p_value threshold', value=0.05, min=0, max=100, step=0.01))
                ),
                loadingButton("netA_path_button", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                HTML('<br><br>'),
                h4(textOutput("netA_path_title")),
                textOutput("netA_pathway_note"),
                HTML('<br>'),
                div(DT::dataTableOutput("path_clust_table")),
                fluidRow(
                  column(2, uiOutput("path_clust_bar_size_input")),
                  column(10, plotOutput("path_clust_plot", inline=TRUE), br())
                ),
                uiOutput("netA_path_bar_note"),
                HTML('<br>'),
                conditionalPanel(
                  condition = "output.path_clust_table !== null && output.path_clust_table !== undefined",
                  actionButton("netA_path_save", "Save Results", icon=icon("download"))),
                HTML('<br>'),
                HTML('<hr />'),
                fluidRow(
                  column(3, tagList(
                    textInput("netA_path_map_id", label="Pathway ID to visualize", value=""),
                    div(style = "font-size: 10px; color: gray; margin-top: -15px;",
                    "The kegg map id. e.g., map00010")))
                  ),
                HTML('<br>'),
                loadingButton("netA_path_vis", "Get the result", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                HTML('<br><br>'),
                uiOutput("netA_path_vis_title"),
                HTML('<br>'),
                fluidRow(
                  column(6, uiOutput("netA_path_vis_edgeUI")),
                  column(6, uiOutput("netA_path_vis_nodeUI"))
                ),
                HTML('<br>'),
                visNetworkOutput("netA_path_graph", width = "100%", height = "800px"),
                conditionalPanel(
                  condition = "output.netA_path_graph !== null && output.netA_path_graph !== undefined",
                  actionButton("netA_path_graph_save", "Save Results", icon=icon("download")))
              ) # box
          ) # row
      ) # Item


    ) # items
  ) # body
) # boardpage

# ===== Server Main =====
server <- function(input, output, session) {

  output$demo_meta1 <- DT::renderDT({
    demoMeta <- data.frame('Accession_ID'=paste0('NIM00',1:5), 'Labels'=c('control', 'CRC', 'control', 'CRC', 'control'))
    demoMeta
  })

  output$demo_meta2 <- DT::renderDT({
    demoMeta <- data.frame('Accession_ID'=paste0('NIM00',1:5), 'Labels'=c('control', 'CRC', 'Adenoma', 'Adenoma', 'control'))
    demoMeta
  })

  output$demo_meta3 <- DT::renderDT({
    demoMeta <- data.frame('Accession_ID'=paste0('NIM00',1:5), 'Labels'=c(82, 64, 72, 83, 57))
    demoMeta
  })

  output$demo_meta4 <- DT::renderDT({
    demoMeta <- data.frame('Accession_ID'=paste0('NIM00',1:5), 'Event'=c(1,0,1,1,0), 'Duration'=c(81, 721, 531, 461, 970))
    demoMeta
  })


  # ===== All Reactive values =====
  # meta values
  values_meta <- reactiveValues(select_task=NULL, meta_v=NULL, groups=NULL)

  # diversity values
  values_diversity <- reactiveValues(alpha_p=NULL, beta_p=NULL, beta_method=NULL,
                                      com_cohort=NULL, com_level=NULL, com_p=NULL,
                                      compair=NULL)

  # auto scanning values
  values_auto <- reactiveValues(x_scaled=NULL)

  # opt lambda values
  values_opt_lmbd <- reactiveValues(opt_lmbd=NULL)

  # run_PreLect values
  values_runprelect <- reactiveValues(PLres=NULL)

  # Network_c values
  values_network_c <- reactiveValues(edge=NULL)

  # Network_f values
  values_network_f <- reactiveValues(network=NULL, node=NULL, edge=NULL)

  # pathway values
  values_pathway <- reactiveValues(target_label_p=NULL, pathway_table=NULL, pathway_plot=NULL,
                                   label_selector_path_net=NULL, path_edge=NULL, path_node=NULL, path_visnet=NULL)

  # network_A values
  values_network_netA <- reactiveValues(target_label=NULL, module_eva=NULL, method_opts=NULL,
                                        clust_forest_df=NULL, clust_forest_plot=NULL, clust_edge=NULL, clust_node=NULL,
                                        clust_edge_vis=NULL, clust_node_vis=NULL, clust_method_vis=NULL,
                                        clust_stat_p=NULL, clust_graph=NULL, # netA
                                        clust_path_label=NULL, clust_id=NULL, clust_id_pathway=NULL, netA_clust_path_plot=NULL, # netA_pathway
                                        clust_vis_mapid=NULL, clust_path_edge=NULL, clust_path_node=NULL, clust_path_visnet=NULL) # netA_pathway_vis


  # ===== main Event 1 (meta) =====
  observeEvent(input$submit_meta_button, {
    meta <- reactive({
      req(input$meta_file)
      meta <- read.csv(input$meta_file$datapath, header = T, stringsAsFactors = F)
      meta
    })

    # save task obj
    obj_type <- input$object_type
    values_meta$select_task <- obj_type # stored in reactive values
    values_meta$meta_v <- meta() # assign to reactiveValues
    # save meta in tmp_files
    write.csv(values_meta$meta_v, file=file.path(tmp_dir, 'meta.csv'), row.names=T)

    output$meta_check1 <- renderText({
      req(input$meta_file)

      if(values_meta$select_task != "Cox"){
        required_columns <- c("Accession_ID", "Labels")
      } else {
        required_columns <- c("Accession_ID", "Event", "Duration")
      }

      submitted_columns <- colnames(meta())
      missing_col <- setdiff(required_columns, submitted_columns)
      message <- "Column names check passed!"
      if (length(missing_col) > 0) {
        message <- paste("The following required column names are missing:", paste(missing_col, collapse = ", "))
      }
      return(message)
    })

    output$meta_check2 <- renderText({
      req(input$meta_file)
      data_sample <- colnames(data)
      submitted_sample <- meta()$Accession_ID
      message <- 'All sample are mapped!'
      missing_sample <- setdiff(data_sample, submitted_sample)
      if(length(missing_sample) > 0){
        message <- paste(
          "The following required sample are missing:",
          paste(missing_sample, collapse = ", ")
        )
      }
      return(message)
    })

    # check if Labels match classification task
    output$meta_check3 <- renderUI({
      req(input$meta_file)
      if(values_meta$select_task == "BC") { # binary classification
        unique_labels <- unique(meta()$Labels)
        if(length(unique_labels) != 2) {
          return(tags$span(style = "color: red;", "Warning: Binary classification task requires 2 unique labels."))
        } else {
          return(tags$span(style = "color: #333;", "The labels are suitable for binary classification."))
        }
      }

      if(values_meta$select_task == "MC") { # multi classification
        unique_labels <- unique(meta()$Labels)
        if( length(unique_labels) < 3 || (length(unique_labels)>10 && is.numeric(meta()$Labels)) ) {
          return(tags$span(style = "color: red;", "Warning: Multi classification task requires at least 3 unique labels or labels belong to regression task."))
        } else {
          return(tags$span(style = "color: #333;", "The labels are suitable for multi classification."))
        }
      }

      if(values_meta$select_task == "Rg") { # regression task
        unique_labels <- unique(meta()$Labels)
        if(length(unique_labels)>10 && is.numeric(meta()$Labels)) {
          return(tags$span(style = "color: #333;", "The labels are suitable for regression task."))
        } else {
          return(tags$span(style = "color: red;", "Warning: Regression task requires numeric labels."))
        }
      }

      if(values_meta$select_task == "Cox") { # time to event
        unique_labels <- unique(meta()$Event)
        if(is.null(unique_labels) || length(unique_labels) != 2) {
          return(tags$span(style = "color: red;", "Warning: Time to event task requires event labels (0 and 1)."))
        } else {
          return(tags$span(style = "color: #333;", "The labels are suitable for time to event."))
        }

      }

    })


    # update label select in network modules
    groups <- unique(values_meta$meta_v$Labels)
    values_meta$groups <- groups # reactive values
    if(values_meta$select_task %in% c("BC")) {
      groups <- c("All", groups)
      updateSelectInput(session, "label_selector_com", choices = groups ) # diversity composition
      updateSelectInput(session, "label_selector", choices = groups ) # network construction
      updateSelectInput(session, "label_selector_f", choices = groups ) # network filter
      updateSelectInput(session, "label_selector_p", choices = groups ) # pathway analysis
      updateSelectInput(session, "label_selector_path_net", choices = groups ) # pathway sub-network
      updateSelectInput(session, "label_selector_netA", choices = groups ) # network analysis cluster
      updateSelectInput(session, "label_netA_path", choices = groups ) # network analysis clustered pathway analysis
    }
    else if(values_meta$select_task == "MC") { # multi classification (not support All)
      updateSelectInput(session, "label_selector_com", choices = c("All", groups) ) # diversity composition
      updateSelectInput(session, "label_selector", choices = groups ) # network construction
      updateSelectInput(session, "label_selector_f", choices = groups ) # network filter
      updateSelectInput(session, "label_selector_p", choices = groups ) # pathway analysis
      updateSelectInput(session, "label_selector_path_net", choices = groups ) # pathway sub-network
      updateSelectInput(session, "label_selector_netA", choices = groups ) # network analysis cluster
      updateSelectInput(session, "label_netA_path", choices = groups ) # network analysis clustered pathway analysis
    }
    else if(values_meta$select_task == "Rg") {
      groups <- c("All", "high", "low")
      updateSelectInput(session, "label_selector_com", choices = c("All") ) # diversity composition
      updateSelectInput(session, "label_selector", choices = groups ) # network construction
      updateSelectInput(session, "label_selector_f", choices = groups ) # network filter
      updateSelectInput(session, "label_selector_p", choices = groups ) # pathway analysis
      updateSelectInput(session, "label_selector_path_net", choices = groups ) # pathway sub-network
      updateSelectInput(session, "label_selector_netA", choices = groups ) # network analysis cluster
      updateSelectInput(session, "label_netA_path", choices = groups ) # network analysis clustered pathway analysis
    }
    else if(values_meta$select_task == "Cox") {
      cox_groups <- as.character(unique(values_meta$meta_v$Event))
      values_meta$groups <- cox_groups # reactive values
      cox_groups <- c("All", cox_groups)
      task_groups <- c("All", "harmful", "protective")
      updateSelectInput(session, "label_selector_com", choices = cox_groups ) # diversity composition
      updateSelectInput(session, "label_selector", choices = task_groups ) # network construction
      updateSelectInput(session, "label_selector_f", choices = task_groups ) # network filter
      updateSelectInput(session, "label_selector_p", choices = task_groups ) # pathway analysis
      updateSelectInput(session, "label_selector_path_net", choices = task_groups ) # pathway sub-network
      updateSelectInput(session, "label_selector_netA", choices = task_groups ) # network analysis cluster
      updateSelectInput(session, "label_netA_path", choices = task_groups ) # network analysis clustered pathway analysis
    }


    # layout
    output$submit_meta <- DT::renderDT({
      req(input$meta_file)
      datatable(values_meta$meta_v, options=list(scrollX=T))
    })

    output$model_usage <- renderText({
      usage <- c("Binary classification", "Multi classification", "Regression", "Time to event")
      names(usage) <- c("BC", "MC", "Rg", "Cox")
      paste('Current model is :', usage[values_meta$select_task])
    })

  }) # end of main Event 1


  # ===== main Event 2 (Diversity) =====

  # alpha
  output$a_color_input <- renderUI({
    # process color related input
    if(is.null(values_meta$select_task)) {
      return(NULL)
    }

    if(values_meta$select_task == "Rg") {
        return(
          HTML("<b>Cannot support to regression task.</b>")
        )
      }

    if(values_meta$select_task != "Rg") {
      color_inputs <- lapply(values_meta$groups, function(group) {
      colourInput(inputId = paste0("color_", group), label = paste("Color for", group), value="white")
      })
      do.call(tagList, color_inputs)
    }

    })

    # size scaler
    output$a_size_input <- renderUI({
      if(is.null(values_meta$select_task)) {
        return(NULL)
      }
      if(values_meta$select_task == "Rg") {
        return(NULL)
      }
      tagList(
        numericInput("a_width", "Width", value = 1100, min = 1, max = 2000),
        numericInput("a_height", "Height", value = 500, min = 1, max = 2000)
      )
    })

  observeEvent(input$alpha_button, {

    # cannot supports to regression task warning msg
    if(values_meta$select_task == "Rg") {
      showNotification("Error: Cannot support to regression task.")
      resetLoadingButton("alpha_button")
      req(FALSE)
    }

    # check prelect_swtich and which data to use
    div_data <- NULL
    if(input$prelect_switch) {
      if(!file.exists(file.path(working_dir, "prelect_dir", "PLres.csv"))) {
        showNotification("PreLect features are not found. Please run Differential taxa first.")
        resetLoadingButton("alpha_button")
        req(FALSE)
      }

      if(values_meta$select_task %in% c("BC", "MC", "Cox")) {
        PL_feat <- read.csv(file.path(working_dir, "prelect_dir", "PLres.csv"), header=T, row.names=1)
        PL_feat <- PL_feat[PL_feat$selected == "Selected",]
        div_data <- data[PL_feat$FeatName,]
      }

    } else {
      div_data <- data
    } # end of prelect switch

    method <- input$alpha_method
    # costum color
    colors <- sapply(values_meta$groups, function(group) {
      input[[paste0("color_", group)]]
    })
    names(colors) <- values_meta$groups

    # alpha df
    if(values_meta$select_task %in%  c("BC", "MC")) { # binary classification, multi classification
      AlphaIndex <- data.frame(shannon = vegan::diversity(t(div_data), index = "shannon"),
                               simpson = vegan::diversity(t(div_data), index = "simpson"),
                               invsimpson = vegan::diversity(t(div_data), index = "invsimpson"),
                               group = values_meta$meta_v$Labels)
    }
    else if(values_meta$select_task == "Cox") { # cox task
      AlphaIndex <- data.frame(shannon = vegan::diversity(t(div_data), index = "shannon"),
                               simpson = vegan::diversity(t(div_data), index = "simpson"),
                               invsimpson = vegan::diversity(t(div_data), index = "invsimpson"),
                               group = factor(values_meta$meta_v$Event))
    }


    # boxplot
    AlphaIndex_long <- AlphaIndex %>%
      pivot_longer(cols = c("shannon", "simpson", "invsimpson"), names_to="method", values_to="value")

    alpha_plot <- ggplot(AlphaIndex_long, aes(x = group, y = value, fill = group)) +
                          geom_boxplot() +
                          geom_signif(comparisons = list(values_meta$groups),
                                      test="wilcox.test", test.args=list(paired=F),
                                      map_signif_level = function(x) paste("p =", scales::pvalue(x)),
                                      textsize=2.5) +
                          scale_fill_manual(values = colors) +
                          labs(x = "Cohort", y = NULL, fill = "Cohort") +
                          theme_minimal() +
                          theme(axis.line = element_line(color = "black", linewidth = 0.5),
                                strip.text = element_text(size = 10, face = "bold"),
                                panel.spacing = unit(4, "lines")) +
                          facet_wrap(~ method, scales = "free_y")


    # reactive values
    values_diversity$alpha_p <- alpha_plot

    if(!is.null(values_diversity$alpha_p)){
      output$alpha_plot <- renderPlot(
        values_diversity$alpha_p,
        width = function() input$a_width,
        height = function() input$a_height,
        res = 100
      )
    }

    resetLoadingButton("alpha_button")
  }) # end of alpha

  # alpha save
  observeEvent(input$alpha_save, {
    a_width <- as.numeric(input$a_width)
    a_height <- as.numeric(input$a_height)
    if(!dir.exists(file.path(working_dir, "diversity"))) {
      dir.create(file.path(working_dir, "diversity"))
    }
    if(file.exists(file.path(working_dir, "diversity", "alpha.png"))) {
      showNotification("The plot is exist.")
      req(FALSE)
    }
    save_path <- file.path(working_dir, "diversity", "alpha.png")
    ggsave(save_path, values_diversity$alpha_p, width=a_width, height=a_height,
            units="px", bg="white", dpi=100)
    showNotification("The plot is saved successfully.")
  })


  # beta
  output$b_color_input <- renderUI({

    if(is.null(values_meta$select_task)) {
      return(NULL)
    }

    if(values_meta$select_task == "Rg") {
      return(
        HTML("<b>Cannot support to regression task.</b>")
      )
    }

    if(values_meta$select_task %in% c("BC", "MC", "Cox")) {
      color_inputs <- lapply(values_meta$groups, function(group) {
      colourInput(inputId = paste0("color_", group), label = paste("Color for", group), value="white")
      })
      do.call(tagList, color_inputs)
    }

  })

  # size scaler
  output$b_size_input <- renderUI({
    if(is.null(values_meta$select_task)) {
      return(NULL)
    }
    if(values_meta$select_task == "Rg") {
      return(NULL)
    }
    tagList(
      numericInput("b_width", "Width", value = 1100, min = 1, max = 2000),
      numericInput("b_height", "Height", value = 500, min = 1, max = 2000)
    )
  })

  observeEvent(input$beta_button, {
    # cannot supports to regression task warning
    if(values_meta$select_task == "Rg") {
      showNotification("Error: Cannot support to regression task.")
      resetLoadingButton("beta_button")
      req(FALSE)
    }

    # check prelect_swtich and which data to use
    div_data <- NULL
    if(input$prelect_switch) {
      if(!file.exists(file.path(working_dir, "prelect_dir", "PLres.csv"))) {
        showNotification("PreLect features are not found. Please run Differential taxa first.")
        resetLoadingButton("beta_button")
        req(FALSE)
      }

      if(values_meta$select_task %in% c("BC", "MC", "Cox")) {
        PL_feat <- read.csv(file.path(working_dir, "prelect_dir", "PLres.csv"), header=T, row.names=1)
        PL_feat <- PL_feat[PL_feat$selected == "Selected",]
        div_data <- data[PL_feat$FeatName,]
      }

    } else {
      div_data <- data
    }

    dim_method <- input$dimension
    values_diversity$beta_method <- dim_method # stored in reactive values
    # custom color
    colors <- sapply(values_meta$groups, function(group) {
      input[[paste0("color_", group)]]
    })
    names(colors) <- values_meta$groups

    # beta df
    beta_plot <- NULL
    data_dist <- vegdist(t(div_data), method="bray")

    # time2event switch
    t2e <- FALSE
    if(values_meta$select_task == "Cox") {
      t2e <- TRUE
    } else { t2e <- FALSE }


    # permanova p-val
    perm_result <- NULL
    if(t2e == FALSE) {
      perm_result <- adonis2(data_dist ~ Labels, data=values_meta$meta_v, permutations = 999) # p_val
    }
    else if(t2e == TRUE) {
      perm_result <- adonis2(data_dist ~ Event, data=values_meta$meta_v, permutations = 999) # p_val
    }
    p_val <- perm_result$`Pr(>F)`

    # pcoa
    if(dim_method == "pcoa") {
      pcoa_result <- pcoa(data_dist)
      pcoa_df <- data.frame(PC1=as.numeric(pcoa_result$vectors[,1]),
                            PC2=as.numeric(pcoa_result$vectors[,2]))

      if(t2e) { # t2e check
        pcoa_df$group <- factor(values_meta$meta_v$Event)
      } else {
        pcoa_df$group <- factor(values_meta$meta_v$Labels)
      }
      pcoa_prop <- round(pcoa_result$values$Relative_eig * 100, 1)

      # plot
      beta_plot <- ggscatter(pcoa_df, x="PC1", y="PC2", combine=T, color="group",
                            ellipse.type = "norm", ellipse = T,ellipse.level = 0.5, ellipse.alpha = 0.5, repel = TRUE) +
                            scale_color_manual(values = colors)+
                            scale_fill_manual(values = colors) +
                            xlab(paste0(c('PC1 (', pcoa_prop[1],'% var.explained)'), collapse = "")) +
                            ylab(paste0(c('PC2 (', pcoa_prop[2],'% var.explained)'), collapse = "")) +
                            labs(tag = paste0("p-val: ", p_val)) +
                            theme(panel.background = element_rect(fill = 'transparent'),
                                  panel.grid = element_blank(),
                                  axis.ticks.length = unit(0.4,"lines"),
                                  axis.ticks = element_line(color='black'),
                                  axis.line = element_line(colour = "black"),
                                  legend.title=element_blank(),
                                  legend.position  = 'right',
                                  plot.tag = element_text(size = 9),
                                  plot.tag.position=c(0.95, 0.02))
    }
    else if(dim_method == "tsne") {
      tsne_res <- Rtsne(t(div_data))
      tsne_plot <- data.frame(t_SNE1 = tsne_res$Y[,1],
                              t_SNE2 = tsne_res$Y[,2])

      if(t2e) { # t2e check
        tsne_plot$group <- factor(values_meta$meta_v$Event)
      } else {
        tsne_plot$group <- factor(values_meta$meta_v$Labels)
      }

      # plot
      beta_plot <- ggscatter(tsne_plot, x="t_SNE1", y="t_SNE2", combine=T, color="group",
                            ellipse.type = "norm", ellipse = T, ellipse.level = 0.5, ellipse.alpha = 0.5, repel = TRUE) +
                            scale_color_manual(values = colors)+
                            scale_fill_manual(values = colors) +
                            xlab('t-SNE1') + ylab("t-SNE2") +
                            labs(tag = paste0("p-val: ", p_val)) +
                            theme(panel.background = element_rect(fill = 'transparent'),
                                  panel.grid = element_blank(),
                                  axis.ticks.length = unit(0.4,"lines"),
                                  axis.ticks = element_line(color='black'),
                                  axis.line = element_line(colour = "black"),
                                  legend.title=element_blank(),
                                  legend.position  = 'right',
                                  plot.tag = element_text(size = 9),
                                  plot.tag.position=c(0.95, 0.02))
    }
    else if(dim_method == "umap") {
      U_plot <- umap(t(div_data), n_components=2, metric="euclidean")
      U_plot <- data.frame(UMAP1=as.numeric(U_plot[,1]),
                          UMAP2=as.numeric(U_plot[,2]),
                          row.names=rownames(U_plot))

      if(t2e) { # t2e check
        U_plot$group <- factor(values_meta$meta_v$Event)
      } else {
        U_plot$group <- factor(values_meta$meta_v$Labels)
      }

      # plot
      beta_plot <- ggscatter(U_plot, x="UMAP1", y="UMAP2", combine=T, color="group",
                            ellipse.type = "norm", ellipse = T,ellipse.level = 0.5, ellipse.alpha = 0.5, repel = TRUE) +
                            scale_color_manual(values = colors)+
                            scale_fill_manual(values = colors) +
                            xlab('UMAP1') + ylab("UMAP2") +
                            labs(tag = paste0("p-val: ", p_val)) +
                            theme(panel.background = element_rect(fill = 'transparent'),
                                  panel.grid = element_blank(),
                                  axis.ticks.length = unit(0.4,"lines"),
                                  axis.ticks = element_line(color='black'),
                                  axis.line = element_line(colour = "black"),
                                  legend.title=element_blank(),
                                  legend.position  = 'right',
                                  plot.tag = element_text(size = 9),
                                  plot.tag.position=c(0.95, 0.02))
    }
    else if(dim_method == "nmds") {
      NMDS <- metaMDS(t(div_data), distance = "bray")
      NMDSplot <- as.data.frame(NMDS$points)

      if(t2e) { # t2e check
        NMDSplot$group <- factor(values_meta$meta_v$Event)
      } else {
        NMDSplot$group <- factor(values_meta$meta_v$Labels)
      }

      # plot
      beta_plot <- ggscatter(NMDSplot, x="MDS1", y="MDS2", combine=T, color="group",
                            ellipse.type = "norm", ellipse = T,ellipse.level = 0.5, ellipse.alpha = 0.5, repel = TRUE) +
                            scale_color_manual(values = colors)+
                            scale_fill_manual(values = colors) +
                            xlab("MDS1") + ylab("MDS2") +
                            labs(tag = paste0("p-val: ", p_val)) +
                            theme(panel.background = element_rect(fill = 'transparent'),
                                  panel.grid = element_blank(),
                                  axis.ticks.length = unit(0.4,"lines"),
                                  axis.ticks = element_line(color='black'),
                                  axis.line = element_line(colour = "black"),
                                  legend.title=element_blank(),
                                  legend.position  = 'right',
                                  plot.tag = element_text(size = 9),
                                  plot.tag.position=c(0.95, 0.02))
    } # end of if else

    # reactive values
    values_diversity$beta_p <- beta_plot
    # ouput plot
    if(!is.null(values_diversity$beta_p)){
      output$beta_plot <- renderPlot(
        values_diversity$beta_p,
        width = function() input$b_width,
        height = function() input$b_height,
        res = 100
      )
    }

    resetLoadingButton("beta_button")
  }) # end of beta

  # beta save
  observeEvent(input$beta_save, {
    b_width <- as.numeric(input$b_width)
    b_height <- as.numeric(input$b_height)
    select_method <- values_diversity$beta_method
    if(!dir.exists(file.path(working_dir, "diversity"))) {
      dir.create(file.path(working_dir, "diversity"))
    }
    if(file.exists(file.path(working_dir, "diversity", paste0(select_method, "_beta.png")))) {
      showNotification("The plot is exist.")
      req(FALSE)
    }
    save_path <- file.path(working_dir, "diversity", paste0(select_method, "_beta.png"))
    ggsave(save_path, values_diversity$beta_p, width=b_width, height=b_height, units="px", bg="white", dpi=100)
    showNotification("The plot is saved successfully.")
  })

  # composition
  # size scaler
  output$com_size_input <- renderUI({
    if(is.null(values_meta$select_task)) {
      return(NULL)
    }
    tagList(
    numericInput("com_width", "Width", value = 1100, min = 1, max = 2000),
    numericInput("com_height", "Height", value = 500, min = 1, max = 2000)
    )
  })

  observeEvent(input$composed_button, {
    target_label <- input$label_selector_com
    sp_map <- read.csv(file.path(script_dir, "taxa_lineage", "species_map.csv"), header=T)
    target_level <- input$target_level_com
    # stored in reactive values
    values_diversity$com_cohort <- target_label
    values_diversity$com_level <- target_level

    # convert to relative abundance percentage
    for(i in 1:ncol(data)){
      data[,i] <- data[,i]/colSums(data)[i]*100
    }

    # make a taxa target level map
    taxa_df <- taxa
    rank_col <- c()
    for(i in seq(dim(taxa_df)[1])) {
      genus <- taxa$Genus[i]
      target <- sp_map[[target_level]] [sp_map$genus == genus][1]
      if(is.null(target)) {
        species <- taxa_df$Species[i]
        cut <- strsplit(species, "_")[1]
        if(is.null(target)) {
          target <- sp_map[[target_level]] [sp_map$genus == cut][1]
        }
      }
      rank_col <- c(rank_col, target)
    }
    taxa_df$target <- rank_col
    taxa_df <- taxa_df %>% mutate(target = ifelse(is.na(target), "Other", target))
    # make a plot df
    uniqtaxa <- unique(taxa_df[,"target"])
    taxaRA <- c()
    for(i in 1:length(uniqtaxa)){
      taxaRA <- c(taxaRA, sum(data[rownames(taxa_df)[taxa_df[,"target"] == uniqtaxa[i]],]))
    }
    names(taxaRA) <- uniqtaxa
    sort_list <- sort(taxaRA, decreasing = T)
    sort_list <- sort_list[!names(sort_list) == "Other"]
    sort_list <- names(sort_list[1:10])


    # make plot dataframe
    compair <- data.frame()
    compair[1,1:4] <- 0
    colnames(compair) <- c("sample", "taxa", "percentage", "group")
    cnt <- 1
    for(i in 1:ncol(data)) {
      for(j in sort_list) {
        compair[cnt, 1] <- colnames(data)[i]
        compair[cnt, 2] <- j
        compair[cnt, 3] <- sum(data[rownames(taxa_df)[taxa_df[,"target"] == j], i])
        compair[cnt, 4] <- ifelse(values_meta$select_task == "Cox", values_meta$meta_v[i, "Event"], values_meta$meta_v[i, "Labels"])
        cnt <- cnt + 1
      }
      compair[cnt, 1] <- colnames(data)[i]
      compair[cnt, 2] <- "Other"
      compair[cnt, 3] <- sum(data[rownames(taxa_df)[!taxa_df[,"target"] %in% sort_list], i])
      compair[cnt, 4] <- ifelse(values_meta$select_task == "Cox", values_meta$meta_v[i, "Event"], values_meta$meta_v[i, "Labels"])
      cnt <- cnt + 1
    }
    compair$taxa <- factor(compair$taxa, levels = c(sort_list, "Other"))
    sorder <- compair[compair$taxa == sort_list[1],]
    sorder <- sorder$sample[order(sorder$percentage, decreasing = T)]
    compair$sample <- factor(compair$sample, levels = sorder)

    # plot
    target_cohort_compair <- NULL
    if(target_label != "All") {
      target_cohort_compair <- compair[compair$group == target_label,]
    } else {
      target_cohort_compair <- compair
    }

    # if sample is too many, do not show geom-line
    com_plot <- NULL
    if(nrow(target_cohort_compair) > 200) {
      com_plot <- ggplot(target_cohort_compair, aes(x = sample, y = percentage, fill = taxa)) +
                    geom_bar(stat="identity") +
                    labs(fill = target_level) + ylab("Relative abundance (%)") +
                    scale_fill_brewer(palette = "Paired") +
                    theme(axis.text.x = element_blank(),
                          axis.line = element_line(linetype = 1,colour = 'black'),
                          panel.background = element_rect(I(0)),
                          panel.grid.major = element_line(colour = NA),
                          panel.grid.minor = element_line(colour = NA),
                          text = element_text(size=16))
    } else {
      com_plot <- ggplot(target_cohort_compair, aes(x = sample, y = percentage, fill = taxa)) +
                    geom_bar(stat="identity", colour = "black") +
                    labs(fill = target_level) + ylab("Relative abundance (%)") +
                    scale_fill_brewer(palette = "Paired") +
                    theme(axis.text.x = element_blank(),
                          axis.line = element_line(linetype = 1,colour = 'black'),
                          panel.background = element_rect(I(0)),
                          panel.grid.major = element_line(colour = NA),
                          panel.grid.minor = element_line(colour = NA),
                          text = element_text(size=16))
    }
    values_diversity$com_p <- com_plot # stored in reactive values

    # output
    if(!is.null(values_diversity$com_p)) {
      output$com_plot <- renderPlot(
        values_diversity$com_p,
        width = function() input$com_width,
        height = function() input$com_height,
        res = 100
      )
    }

    resetLoadingButton("composed_button")
  }) # end of composition

  # composition save
  observeEvent(input$com_save, {
    cohort <- values_diversity$com_cohort
    taxa_level <- values_diversity$com_level
    com_width <- as.numeric(input$com_width)
    com_height <- as.numeric(input$com_height)
    if(!dir.exists(file.path(working_dir, "diversity"))) {
      dir.create(file.path(working_dir, "diversity"))
    }
    if(file.exists(file.path(working_dir, "diversity", paste(c(cohort, taxa_level, "composition.png"), collapse="_")))) {
      showNotification("The plot is exist.")
      req(FALSE)
    }
    save_path <- file.path(working_dir, "diversity", paste(c(cohort, taxa_level, "composition.png"), collapse="_"))
    ggsave(save_path, values_diversity$com_p, width=com_width, height=com_height, units="px", bg="white", dpi=100)
    showNotification("The plot is saved successfully.")
  })


  # ===== main Event 3 (PreLect) =====
  # PreLect auto scanning
  observeEvent(input$tuning_run, {

    meta_v <- values_meta$meta_v
    sp <- intersect(colnames(data),  meta_v$Accession_ID)
    X_raw <- data[, sp]

    meta_ <- meta_v
    rownames(meta_) <- meta_$Accession_ID
    sub_meta <- meta_[sp, ]

    if(input$norm_method == 'z'){
      X_scaled <- t(scale(t(X_raw)))
    } else if(input$norm_method == 'zm') {
      X_scaled <- t(scale(t(X_raw)))
      values_auto$x_scaled <- X_scaled # stored in reactiveValues
      X_scaled <- apply(X_scaled, 2, function(col) {(col - min(col)) / (max(col) - min(col))})
    } else if(input$norm_method == 'clr'){
      # Centered log ratio transform
      clr_function <- function(x) {
        return( log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))) )
      }
      X_scaled <- clr_function(X_raw)
    }


    step_ = as.integer(input$auto_step)
    mi_ = as.integer(input$max_iter)
    lr_ = input$lr
    # create dir
    save_dir = file.path(working_dir, 'prelect_dir')


    if(values_meta$select_task == 'BC'){
      lrange <- AutoScanning(X_scaled, X_raw, sub_meta$Labels, task="classification", step=step_, max_iter=mi_, lr=lr_)
      tuning_res <- LambdaTuning(X_scaled, X_raw, sub_meta$Labels, lrange, outpath=save_dir, max_iter=mi_, lr=lr_)
    }

    if(values_meta$select_task == 'MC') {
      lrange <- AutoScanningMultiClass(X_scaled, X_raw, sub_meta$Labels, step=step_, max_iter=mi_, lr=lr_)
      tuning_res <- LambdaTuningMultiClass(X_scaled, X_raw, sub_meta$Labels, lrange, outpath=save_dir, max_iter=mi_, lr=lr_)
    }

    if(values_meta$select_task == 'Rg') {
      lrange <- AutoScanning(X_scaled, X_raw, sub_meta$Labels, task="regression", step=step_, max_iter=mi_, lr=lr_)
      tuning_res <- LambdaTuning(X_scaled, X_raw, sub_meta$Labels, task="regression", lrange, outpath=save_dir, max_iter=mi_, lr=lr_)
    }

    if(values_meta$select_task == 'Cox') {
      lrange <- AutoScanningCoxPH(X_scaled, X_raw, sub_meta$Event, sub_meta$Duration, step=step_, max_iter=mi_, lr=lr_)
      tuning_res <- LambdaTuningCoxPHParallel(X_scaled, X_raw, sub_meta$Event, sub_meta$Duration, lrange,
                                              outpath=save_dir, max_iter=mi_, lr=lr_, n_cores=threads)
    }

    # saved X_scaled in prelect_dir
    write.csv(X_scaled, file=file.path(save_dir, "dataStd.csv"), row.names=T)
    showNotification("PreLect lambda tuning finished!")

    resetLoadingButton("tuning_run")
  }) # end of auto_scan



  # PreLect optimized lambda
  observeEvent(input$opt_lambda, {

    maxdepth_ = as.integer(input$max_depth)
    minbucket_ = as.integer(input$min_bucket)
    save_dir = file.path(working_dir, 'prelect_dir')

    if(dir.exists(save_dir)){
      d1 <- read.csv(file.path(save_dir, 'TuningResult.csv'))
      d2 <- read.csv(file.path(save_dir,'Pvl_distribution.csv'))
      lmbd_picking <- LambdaDecision(d1, d2, maxdepth=maxdepth_, minbucket=minbucket_)

      # stored in values
      values_opt_lmbd$opt_lmbd <- lmbd_picking$opt_lmbd

      output$lambda_dec_plot <- renderPlot({
        lmbd_picking$selected_lmbd_plot/lmbd_picking$pvl_plot
      })

      output$opt_lambda_value <- renderText({
        paste0('optimized lambda is ', lmbd_picking$opt_lmbd)
      })

    } else {
      output$lambda_dec_plot <- renderPlot({
        NULL
      })

      output$opt_lambda_value <- renderText({
        'No result to present. Please run the Lambda Tuning'
      })

    }

    # update opt_lmbd selector
    updateNumericInput(session, "opt_lmbd", label='optimized lambda', value=lmbd_picking$opt_lmbd, min=0, max=10, step=1)

    resetLoadingButton("opt_lambda")
  }) # end of opt_lambda


  # PreLect visualize
  observeEvent(input$run_prelect, {

    # X_scaled
    PL_dir <- file.path(working_dir, "prelect_dir")
    X_scaled <- read.csv(file.path(PL_dir, "dataStd.csv"), row.names=1)
    # meta
    meta_ <- values_meta$meta_v
    sp <- intersect(colnames(X_scaled),  meta_$Accession_ID)
    X_raw <- data[, sp]
    X_scaled <-  X_scaled[, sp]
    sub_meta <- meta_[meta_$Accession_ID %in% sp, ]
    # opt lambda
    lmbd <- input$opt_lmbd
    # task
    task <- values_meta$select_task

    # output reactive
    visual_table <- reactive({

      pvlvec <- GetPrevalence(X_raw)

      if(task == 'BC'){
        PL_out <- PreLect(X_scaled, pvlvec, sub_meta$Labels, lambda=lmbd, task="classification")
        visual_table <- FeatureProperty(X_raw, sub_meta$Labels, PL_out, task="classification")
      }

      if(task == 'Rg') {
        PL_out <- PreLect(X_scaled, pvlvec, sub_meta$Labels, lambda=lmbd, task="regression")
        visual_table <- FeatureProperty(X_raw, sub_meta$Labels, PL_out, task="regression")
      }

      if(task == 'MC') {
        PL_out <- PreLectMultiClass(X_scaled, X_raw, sub_meta$Labels, lambda=lmbd)
        visual_table <- FeatureProperty(X_raw, sub_meta$Labels, PL_out, task="multi-class classification")
      }

      if(task == "Cox") {
        PL_out <- PreLectCoxPH(X_scaled, pvlvec, sub_meta$Event, sub_meta$Duration, lambda=lmbd)
        visual_table <- FeatureProperty(X_raw, sub_meta$Event, PL_out)
      }

      visual_table
    }) # end of reactive

    # stored reactive values
    values_runprelect$PLres <- visual_table()
    # saved res in tmp_files
    write.csv(values_runprelect$PLres, file=file.path(tmp_dir, 'PLres.csv'), row.names=T)

    numeric_cols <- colnames(values_runprelect$PLres)[sapply(values_runprelect$PLres, is.numeric)]
    output$PLres_table <- DT::renderDT({
      datatable(values_runprelect$PLres, options=list(scrollX=T)) %>%
                formatSignif(columns = numeric_cols, digits = 3)
    })

    # render text
    num_total <- nrow(values_runprelect$PLres)
    num_feature <- sum(values_runprelect$PLres$selected == "Selected")
    selected_percent <- num_feature/num_total*100

    output$PLres_note1 <- renderText({
      paste0("The number of total features : ", num_total)
    })

    output$PLres_note2 <- renderText({
      paste0("The number of selected features : ", num_feature, " ( " , round(selected_percent, 2), "% )")
    })

    showNotification("PreLect feature selection is finished!")

    resetLoadingButton("run_prelect")
  }) # end of run_PreLect


  # save PLres.csv
  observeEvent(input$PLres_save, {
    pl_dir <- file.path(working_dir, "prelect_dir")
    write.csv(values_runprelect$PLres, file=file.path(pl_dir, "PLres.csv"), row.names=T)
    showNotification("Results are saved in prelect_dir successfully!")
  })



  # ===== main Event 4 (Network) =====
  # Network construction
  observeEvent(input$make_network_button, {
    select_task <- values_meta$select_task
    target_label <- input$label_selector
    adj_p <- input$adj_p
    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "network_construction.py"),
                              working_dir, target_label, adj_p, threads, select_task))
    if( exit_code != 0){
      showNotification("Error occurred in network construction. Please check log in terminal.")
      resetLoadingButton("make_network_button")
      req(FALSE)
    }

    # read table
    network_con_dir <- file.path(working_dir, "tmp_files")
    edge <- read.csv(file.path(tmp_dir, paste0(target_label, "_edge.csv")))
    # Reactive values
    values_network_c$edge <- edge

    # visualize
    output$graph_title <- renderText({
      paste("Whole", target_label, "tendency network")
    })

    output$edge_table <- DT::renderDT({
      datatable(values_network_c$edge, options=list(scrollX=T, columnDefs =list(list(visible=FALSE, targets=c(-1))))) %>%
        formatSignif(columns = c('weight', 'p', 'p_adj'), digits = 3)
    })

    resetLoadingButton("make_network_button")
  }) # end of Network_construction

  # network saved button
  observeEvent(input$network_save, {
    target_label <- input$label_selector
    net_dir <- file.path(working_dir, "network_construction")
    if(!dir.exists(net_dir)){
      dir.create(net_dir)
    }

    if(!file.exists(file.path(net_dir, paste0(target_label, "_edge.csv")))) {
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_edge.csv")), to=net_dir)
      showNotification("Network is saved in network_construction successfully!")
    } else {
      showNotification(paste0(target_label, " network already exists. Cannot save the results."))
      req(FALSE)
    }
  }) # end of network_save


  # Network filter
  observeEvent(input$filter_network_button, {
    pl_weight <- input$PL_w
    edge_weight <- input$edge_w
    size <- input$node_size
    rank <- input$select_rank
    target_label_f <- input$label_selector_f
    size_range_min <- input$size_range[1]
    size_range_max <- input$size_range[2]
    width_range_min <- input$width_range[1]
    width_range_max <- input$width_range[2]
    select_task <- values_meta$select_task



    if(!file.exists(file.path(working_dir, "network_construction", paste0(target_label_f, "_edge.csv")))){
      showNotification(
        paste0("Cannot find the networt construction result of ", target_label_f, ". Please run network construction first."))
        resetLoadingButton("filter_network_button")
        req(FALSE)
      }

    # main
    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "network_filter.py"),
                          working_dir, target_label_f, size, pl_weight, edge_weight, rank,
                          size_range_min, size_range_max, width_range_min, width_range_max, select_task))
    if( exit_code != 0){
      showNotification("Error occurred in network filter. Maybe too large threshold leads to the empty network. Please check log in terminal.")
      resetLoadingButton("filter_network_button")
      req(FALSE)
    }

    # read table
    network_f_dir <- file.path(working_dir, "tmp_files")
    node_f <- read.csv(file.path(tmp_dir, paste0(target_label_f, "_node_f.csv")))
    edge_f <- read.csv(file.path(tmp_dir, paste0(target_label_f, "_edge_f.csv")))
    # Reactive values
    values_network_f$node <- node_f
    values_network_f$edge <- edge_f

    if(rank == 'none') {
      visnet <- visNetwork(node_f, edge_f, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            background = "#0085AF",
            border = "#013848",
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

    } else {
      # deal with colored by groups
      lnodes <- data.frame(label=unique(node_f$group), shape=c("dot"), color=col_vector[1:length(unique(node_f$group))])
      color_v <- lnodes$color[match(node_f$group, lnodes$label)]
      node_to_plot <- node_f; node_to_plot$color <- color_v
      visnet <- visNetwork(node_to_plot, edge_f, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLegend(useGroups = F, width=0.2, addNodes = lnodes, position = "right")
    }
    # reactive Values
    values_network_f$network <- visnet


    # visualize
    output$graph_f_title <- renderText({
      paste("Filtered", target_label_f, "tendency network")
    })

    output$filter_note <- renderText({
      paste0("( by PreLect weight : ", pl_weight, ", edge weight : ", edge_weight, " )")
    })

    # left panel
    output$edge_f_ui <- renderUI({
      tagList(
        HTML("<h4><b>Edge</b></h4>"),
        DT::dataTableOutput("edge_f_table")
      )
    })

    # right panel
    output$node_f_ui <- renderUI({
      tagList(
        HTML("<h4><b>Node</b></h4>"),
        DT::dataTableOutput("node_f_table")
      )
    })

    output$edge_f_table <- DT::renderDT({
      datatable(values_network_f$edge, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(6,7))))) %>%
                  formatSignif(columns = c('weight', 'p', 'p_adj'), digits = 3)
    })

    output$node_f_table <- DT::renderDT({
      datatable(values_network_f$node, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,7))))) %>%
                  formatSignif(columns = c('raw_value'), digits = 3)
    })

    output$network_graph <- renderVisNetwork({
      visnet
    }) # end of renderVisNetwork

    resetLoadingButton("filter_network_button")
  }) # end of Network filter

  # network saved button
  observeEvent(input$network_f_save, {
    target_label_f <- input$label_selector_f
    net_f_dir <- file.path(working_dir, "network_filter")
    if(!dir.exists(net_f_dir)){
      dir.create(net_f_dir)
    }
    save_dir <- file.path(net_f_dir, target_label_f)

    if (!dir.exists(save_dir)){
      dir.create(save_dir)
      visSave(values_network_f$network, file = file.path(save_dir, paste0(target_label_f, ".html")))
      file.copy(from=file.path(tmp_dir, paste0(target_label_f, "_edge_f.csv")), to=save_dir)
      file.copy(from=file.path(tmp_dir, paste0(target_label_f, "_node_f.csv")), to=save_dir)
      showNotification("Results are saved in network_filter successfully!")
    } else {
      showNotification(paste0("Directory of ", target_label_f, " already exists. Cannot save the result."))
      req(FALSE)
    }

  }) # end of network_save




  # ===== main Event 5 (pathway) =====
  output$path_bar_size_input <- renderUI({
    if (is.null(values_pathway$pathway_plot)) {
      return(NULL)
    }
    tagList(
      numericInput("path_bar_width", "Width", value = 1100, min = 1, max = 2000),
      numericInput("path_bar_height", "Height", value = 500, min = 1, max = 2000)
    )
  })


  observeEvent(input$pathway_button, {
    target_label_p <- input$label_selector_p
    values_pathway$target_label_p <- target_label_p
    adj_p <- input$pathway_p

    # check network filter table
    if(!file.exists(file.path(working_dir, "network_filter", target_label_p, paste0(target_label_p, "_edge_f.csv")))){
      showNotification(
        paste0("Cannot find the network filter result of ", target_label_p, ". Please run network filter first."))
        resetLoadingButton("pathway_button")
        req(FALSE)
    }

    # main
    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "pathway_analysis.py"),
                          working_dir, target_label_p, adj_p, threads))
    if( exit_code != 0){
      showNotification("Error occurred in pathway analysis. Please check log in terminal.")
      resetLoadingButton("pathway_button")
      req(FALSE)
    }

    # read table
    pathway_table <- read.csv(file.path(tmp_dir, paste0(target_label_p, "_pathway.csv")))
    # Reactive values
    values_pathway$pathway_table <- pathway_table

    if(nrow(pathway_table) != 0) {

      # bar plots
      # node
      node_plot <- data.frame(Pathway=values_pathway$pathway_table$annotation,
                              OR=values_pathway$pathway_table$step1_stats,
                              p=values_pathway$pathway_table$step1_p_adj,
                              Type="Node")

      # edge
      edge_plot <- data.frame(Pathway=values_pathway$pathway_table$annotation,
                              OR=values_pathway$pathway_table$step2_stats,
                              p=values_pathway$pathway_table$step2_p_adj,
                              Type="Edge")

      combined_plot <- rbind(node_plot, edge_plot)
      combined_plot <- combined_plot[order(combined_plot$Pathway, combined_plot$OR), ]
      p <- ggplot(combined_plot, aes(x = reorder(Pathway, OR), y = OR, fill = p, group = Type)) +
        scale_fill_gradient(low = "red", high = "pink") +
        coord_flip() +
        labs(x = "Pathway", y = "Odds Ratio", fill = "p-value") +
        theme_minimal() +
        theme(axis.line = element_line(color = "black", linewidth = 0.5),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10))

      for (i in 1:(length(node_plot$Pathway))) {
        p <- p + annotate('rect', xmin = i-0.5, xmax = i+0.5, ymin = -Inf, ymax = Inf,
                          fill = ifelse(i %% 2 == 0, 'white', 'gray95'))}

      p <- p + geom_bar(stat = "identity", position = position_dodge(width = 0.6), width=0.6, colour="black") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      values_pathway$pathway_plot <- p # stored in reactive values

    } else {
      values_pathway$pathway_plot <- NULL
    }

    # output
    output$pathway_title <- renderText({
      paste("Pathway analysis of", target_label_p)
    })

    output$pathway_note <- renderText({
      paste0("( adjusted p_value threshold : ", adj_p, " )")
    })

    output$vital_pathway <- DT::renderDT({
      datatable(values_pathway$pathway_table, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,6))))) %>%
                  formatSignif(columns = c('step1_stats', 'step1_p_adj', 'step2_stats', 'step2_p_adj'), digits = 3)
    })

    output$path_plot <- renderPlot(
      values_pathway$pathway_plot,
      width = function() input$path_bar_width,
      height = function() input$path_bar_height,
      res = 100
    )

    output$bar_note <- renderUI({
      if (is.null(values_pathway$pathway_plot)) {
        return(NULL)
      }

      if(!is.null(values_pathway$pathway_plot)) {
        return(HTML("<ul>
          <li>The first bar: odds ratio of node.</li>
          <li>The second bar: odds ratio of edge.</li>
        </ul>"))
      }
    })

    resetLoadingButton("pathway_button")
  }) # end of pathway_button


  # save pathway save button
  observeEvent(input$pathway_save, {
    target_label_p <- values_pathway$target_label_p
    pathway_dir <- file.path(working_dir, "pathway_analysis")
    if(!dir.exists(pathway_dir)){
      dir.create(pathway_dir)
    }

    if(!file.exists(file.path(pathway_dir, paste0(target_label_p, "_pathway.csv")))){
      file.copy(from=file.path(tmp_dir, paste0(target_label_p, "_pathway.csv")), to=pathway_dir)
      ggsave(filename=file.path(pathway_dir, paste0(target_label_p, "_pathway_plot.png")),
              plot=values_pathway$pathway_plot, bg="white",
              width=input$path_bar_width, height=input$path_bar_height, units="px", dpi=100)
      showNotification("Result is saved in pathway_analysis successfully!")
    } else {
      showNotification(paste0("Directory of ", target_label_p, " already exists. Cannot save the result."))
      req(FALSE)
    }

  }) # end of pathway save button

  # pathway sub-network
  observeEvent(input$pathway_net_button, {
    target_label <- input$label_selector_path_net
    values_pathway$label_selector_path_net <- target_label
    target_pathway <- input$map_id

    # check network filter results
    if(!dir.exists(file.path(working_dir, "network_filter", target_label))) {
      showNotification(
        paste0("Cannot find the filtered sub-network result of ", target_label, ". Please run network filter first."))
        resetLoadingButton("pathway_net_button")
        req(FALSE)
    }

    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "path_net.py"),
                          working_dir, target_pathway, target_label))
    if(exit_code != 0){
      showNotification("Error occurred in pathway sub-network visualization. Please check log in terminal.")
      resetLoadingButton("pathway_net_button")
      req(FALSE)
    }

    # read table
    path_edge <- read.csv(file.path(tmp_dir, paste0(target_label, "_", target_pathway, "_edge.csv")))
    path_node <- read.csv(file.path(tmp_dir, paste0(target_label, "_", target_pathway, "_node.csv")))
    # rds file
    ko_map <- readRDS(file.path(script_dir, "whole_komap", "total_KO_pathway.rds"))
    annotated_name <- ko_map[ko_map$mapid == target_pathway, "pathway"][1]

    # stored in reactive values
    values_pathway$path_edge <- path_edge
    values_pathway$path_node <- path_node

    # visualize network
    if(!"group" %in% colnames(path_node)) {
      visnet <- visNetwork(path_edge, path_node, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            background = "#0085AF",
            border = "#013848",
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

    } else {
      # deal with colored by groups
      lnodes <- data.frame(label=unique(path_node$group), shape=c("dot"), color=col_vector[1:length(unique(path_node$group))])
      color_v <- lnodes$color[match(path_node$group, lnodes$label)]
      node_to_plot <- path_node; node_to_plot$color <- color_v    # assgin color column
      visnet <- visNetwork(node_to_plot, path_edge, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLegend(useGroups = F, width=0.2, addNodes = lnodes, position = "right")
    }

    values_pathway$path_visnet <- visnet  # stored in reactive values

    # output
    output$pathway_net_title <- renderUI({
      tagList(
        HTML(paste0("The selected pathway is : <b>", annotated_name,"</b><br>")),
        HTML(paste0("The selected cohort is : <b>", target_label, "</b>"))
      )
    })

    output$pathway_net_edgeUI <- renderUI({
      tagList(
        HTML("<h4><b>Edge</b></h4>"),
        DT::dataTableOutput("pathway_net_edge")
      )
    })

    output$pathway_net_nodeUI <- renderUI({
      tagList(
        HTML("<h4><b>Node</b></h4>"),
        DT::dataTableOutput("pathway_net_node")
      )
    })

    output$pathway_net_edge <- DT::renderDT({
      datatable(values_pathway$path_edge, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(6,7))))) %>%
                  formatSignif(columns = c('weight', 'p', 'p_adj'), digits = 3)
    })

    output$pathway_net_node <- DT::renderDT({
      datatable(values_pathway$path_node, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,7))))) %>%
                  formatSignif(columns = c('raw_value'), digits = 3)
    })

    output$pathway_net_graph <- renderVisNetwork({
      values_pathway$path_visnet
    }) # end of renderVisNetwork

    resetLoadingButton("pathway_net_button")
  }) # end of pathway sub-network

  # save pathway sub-network
  observeEvent(input$pathway_net_save, {
    target_label <- values_pathway$label_selector_path_net
    target_pathway <- input$map_id
    path_net_dir <- file.path(working_dir, "pathway_analysis")
    if(!dir.exists(path_net_dir)){
      dir.create(path_net_dir)
    }
    save_dir <- file.path(path_net_dir, paste0(target_label, "_", target_pathway))

    if (!dir.exists(save_dir)){
      dir.create(save_dir)
      visSave(values_pathway$path_visnet, file = file.path(save_dir, paste0(target_label, "_", target_pathway, ".html")))
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", target_pathway, "_edge.csv")), to=save_dir)
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", target_pathway, "_node.csv")), to=save_dir)
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", target_pathway, "_koset.txt")), to=save_dir)
      showNotification("Results are saved in pathway_analysis successfully!")
    } else {
      showNotification(paste0("Directory of ", target_pathway, " already exists. Cannot save the result."))
      req(FALSE)
    }

  }) # end of pathway save button


  # ===== main Event 6 (Network analysis) =====
  # vis func
  clust_cnt_vis <- function(cnt_table, clust_method) {
    Xaxis <- ifelse(cnt_table$Var1 == "unclustered", "unclustered", as.character(1:length(cnt_table$Var1)))
    plot <- ggplot(cnt_table, aes(x = Var1, y = Freq, fill = Method)) +
      geom_bar(stat = "identity", colour = "black") +
      geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
      labs(x = "Cluster ID", y = "Taxa counts", fill = "Method", title = clust_method) +
      theme_minimal() +
      theme(axis.line = element_line(color = "black", linewidth = 0.5),
            axis.text.x = element_text(size=10), legend.position="none",
            plot.title = element_text(hjust = 0.5, color = "dimgrey", size = 12, face = "bold")) +
      scale_x_discrete(labels = Xaxis)
    return(plot)
  }

  observeEvent(input$cluster_out, {
    target_label_netA <- input$label_selector_netA
    values_network_netA$target_label <- target_label_netA
    cutoff <- input$clust_cutoff

    # main
    if(!file.exists(file.path(working_dir, "network_filter"))) {
      showNotification("Cannot find the directory of networt filter. Please run network filter first.")
      resetLoadingButton("cluster_out")
      req(FALSE)
    }
    else if(!file.exists(file.path(working_dir, "network_filter", target_label_netA))) {
      showNotification(paste0("Cannot find the result of ", target_label_netA,". Please run network filter first."))
      resetLoadingButton("cluster_out")
      req(FALSE)
    }

    edge_NetA <- read.csv(file.path(working_dir, "network_filter", target_label_netA, paste0(target_label_netA,"_edge_f.csv")))
    node_NetA <- read.csv(file.path(working_dir, "network_filter", target_label_netA, paste0(target_label_netA,"_node_f.csv")))
    # store edge_NetA in reactive values
    values_network_netA$clust_edge <- edge_NetA

    net <- make_net(edge_NetA) # adj. matrix
    net_pack <- graph_from_biadjacency_matrix(net, mode="total", weighted=TRUE) # igraph obj

    # cluster
    clust_methods <- c()
    modularities <- c()
    stat_tables <- list()
    plots <- list()
    method_opts <- c() # for cluster pathway analysis id selectors

    if(input$B_switch) {
      clustB <- cluster_edge_betweenness(net_pack)
      clust_methods <- c(clust_methods, "Edge betweenness")
      modularities <- c(modularities, modularity(clustB))
      method_opts <- c(method_opts, "Betweenness")
      # statistic
      node_NetA <- assignModule(node_NetA,clustB,"Betweenness", cutoff) # add columns in exist node table
      B_cnt <- as.data.frame(table(node_NetA$Betweenness)) # statistic table
      B_cnt$Method <- "Betweenness"
      stat_tables <- append(stat_tables, list(B_cnt))
      B_plot <- clust_cnt_vis(B_cnt, "Betweenness")
      plots <- append(plots, list(B_plot))
    }

    if(input$W_switch) {
      clustW <- cluster_walktrap(net_pack)
      clust_methods <- c(clust_methods, "Random walks")
      modularities <- c(modularities, modularity(clustW))
      method_opts <- c(method_opts, "RandomWalks")
      # statistic
      node_NetA <- assignModule(node_NetA,clustW,"RandomWalks", cutoff) # add columns in exist node table
      W_cnt <- as.data.frame(table(node_NetA$RandomWalks)) # statistic table
      W_cnt$Method <- "RandomWalks"
      stat_tables <- append(stat_tables, list(W_cnt))
      W_plot <- clust_cnt_vis(W_cnt, "RandomWalks")
      plots <- append(plots, list(W_plot))
    }

    if(input$L_switch) {
      clustL <- cluster_louvain(net_pack)
      clust_methods <- c(clust_methods, "Louvain")
      modularities <- c(modularities, modularity(clustL))
      method_opts <- c(method_opts, "Louvain")
      # statistic
      node_NetA <- assignModule(node_NetA,clustL,"Louvain", cutoff) # add columns in exist node table
      L_cnt <- as.data.frame(table(node_NetA$Louvain)) # statistic table
      L_cnt$Method <- "Louvain"
      stat_tables <- append(stat_tables, list(L_cnt))
      L_plot <- clust_cnt_vis(L_cnt, "Louvain")
      plots <- append(plots, list(L_plot))
    }

    if(input$G_switch) {
      clustG <- cluster_fast_greedy(net_pack)
      clust_methods <- c(clust_methods, "Greedy")
      modularities <- c(modularities, modularity(clustG))
      method_opts <- c(method_opts, "Greedy")
      # statistic
      node_NetA <- assignModule(node_NetA,clustG,"Greedy", cutoff) # add columns in exist node table
      G_cnt <- as.data.frame(table(node_NetA$Greedy)) # statistic table
      G_cnt$Method <- "Greedy"
      stat_tables <- append(stat_tables, list(G_cnt))
      G_plot <- clust_cnt_vis(G_cnt, "Greedy")
      plots <- append(plots, list(G_plot))
    }

    if(input$LE_switch) {
      clustLE <- cluster_leading_eigen(net_pack)
      clust_methods <- c(clust_methods, "Non-negative eigenvector")
      modularities <- c(modularities, modularity(clustLE))
      method_opts <- c(method_opts, "Leading_eigen")
      # statistic
      node_NetA <- assignModule(node_NetA,clustLE,"Leading_eigen", cutoff) # add columns in exist node table
      LE_cnt <- as.data.frame(table(node_NetA$Leading_eigen)) # statistic table
      LE_cnt$Method <- "Leading_eigen"
      stat_tables <- append(stat_tables, list(LE_cnt))
      LE_plot <- clust_cnt_vis(LE_cnt, "Leading_eigen")
      plots <- append(plots, list(LE_plot))
    }

    Modularity_data <- data.frame(Method = clust_methods, Modularity = modularities)
    # store modularity in reactive values
    values_network_netA$module_eva <- Modularity_data

    # store updated clustered result overview in reactive values
    values_network_netA$clust_node <- node_NetA

    # combine the bar plot of each method
    combined_plot <- wrap_plots(plots) + plot_layout(ncol=3, nrow=2)
    values_network_netA$clust_stat_p <- combined_plot

    # store user method options
    values_network_netA$method_opts <- method_opts

    # update clustered pathway ID_label and stored cluster info in tmp_files
    id_vec <- c()
    for(m in values_network_netA$method_opts){
      all_id <- unique(values_network_netA$clust_node[[m]])
      id_vec <- c(id_vec, all_id)
    }
    id_vec <- unique(id_vec[! id_vec %in% c("unclustered")]) # remove unclustered
    updateSelectInput(session, "netA_path_clust_id", choices=id_vec)


    # output
    output$modularity_title <- renderUI({
      HTML(paste0("<b>Modularity of ", target_label_netA, "</b>"))
    })

    output$modularity_table <- DT::renderDT({
      datatable(values_network_netA$module_eva, options=list(scrollX=T)) %>%
                  formatSignif(columns = c('Modularity'), digits = 3)
    })

    output$clustNode_title <- renderUI({
      HTML(paste0("<b>Cluster results of ", target_label_netA, "</b>"))
    })

    output$clustNode_table <- DT::renderDT({
      datatable(values_network_netA$clust_node, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,7))))) %>%
                  formatSignif(columns = c('raw_value'), digits = 3)
    })

    output$clust_stat_title <- renderUI({
      HTML(paste0("<b>Cluster stastistic of ", target_label_netA, "</b>"))
    })

    output$clust_stat_plot <- renderPlot({
      values_network_netA$clust_stat_p
    })

    resetLoadingButton("cluster_out")
  }) # end of network analysis button


  # save cluster stat results
  observeEvent(input$clust_stat_save, {
    target_label <- values_network_netA$target_label

    netA_dir <- file.path(working_dir, "network_analysis")
    if(!dir.exists(netA_dir)){
      dir.create(netA_dir)
    }

    save_dir <- file.path(netA_dir, paste0(target_label, "_cluster_overview"))
    if(!dir.exists(save_dir)) {
      dir.create(save_dir)
      write.csv(values_network_netA$clust_node, file=file.path(save_dir, paste0(target_label, "_clust_node.csv")), row.names=F)
      write.csv(values_network_netA$module_eva, file=file.path(save_dir, paste0(target_label, "_modularity.csv")), row.names=F)
      ggsave(filename=file.path(save_dir, paste0(target_label, "_clust_stat_plot.png")),
              plot=values_network_netA$clust_stat_p, bg="white", width=10, height=8, units="in")
      showNotification("Results are saved in network_analysis successfully!")
    } else {
      showNotification(paste0("Directory of ", target_label, " already exists. Cannot save the result."))
      req(FALSE)
    }

  }) # end of clust_stat_save

  # generate clustered network
  observeEvent(input$make_clust_network, {
    clust_select <- input$clust_select
    values_network_netA$clust_method_vis <- clust_select # reactive values
    clust_edge <- values_network_netA$clust_edge
    clust_node <- values_network_netA$clust_node

    # forest plot (fisher's exact test)
    make_forest <- function(node_table, tend, total_cluster) {
      test_list <- list()
      for(cl in total_cluster) {
        tend_vec <- ifelse(node_table$tendency == tend, 1, 0)
        cl_vec <- ifelse(node_table[[clust_select]] == cl, 1, 0)
        contingency_table <- table(tend_vec, cl_vec)
        fisher_result <- fisher.test(contingency_table)
        conf_int <- fisher_result$conf.int
        OR <- fisher_result$estimate
        p_val <- fisher_result$p.value
        # add in list
        test_list[[cl]] <- list(OR = OR, p_val = p_val, conf_int = conf_int)
      }
      forest <- do.call(rbind, lapply(names(test_list), function(name) {
        data.frame(
          clust_id = name,
          OR = test_list[[name]]$OR,
          p_val = test_list[[name]]$p_val,
          CI_lower = test_list[[name]]$conf_int[1],
          CI_upper = test_list[[name]]$conf_int[2],
          group = tend
        )
      }))
      return (forest)
    }  # end of function

    # call function
    total_tend <- unique(clust_node$tendency)
    total_clust <- unique(clust_node[[clust_select]])

    # check the number of clusters
    clust_check <- FALSE
    if(length(total_clust) < 2) {
      showNotification("Warning: The number of clusters is less than 2. Cannot perform correlation results.")
      clust_check <- FALSE
    } else {
      clust_check <- TRUE
    }
    # check the number of tendency
    tend_check <- FALSE
    if(length(total_tend) < 2) {
      showNotification("Warning: The number of tendency is less than 2. Cannot perform correlation results.")
      tend_check <- FALSE

    } else if(clust_check == TRUE) {

      tend_forest <- list()
      for(t in total_tend) {
        tend_clust_df <- make_forest(clust_node, t, total_clust)
        tend_forest[[t]] <- tend_clust_df # add in list
      }

      # combine the two datasets
      forest_df <- do.call("rbind", tend_forest)
      forest_df$clust_id <- factor(forest_df$clust_id)
      rownames(forest_df) <- 1:nrow(forest_df)

      # plot
      # define colours for dots and bars
      num_groups <- length(unique(forest_df$group))
      if(num_groups == 2) {
        dotCOLS = c("#a6d8f0","#f9b282")
        barCOLS = c("#008fd5","#de6b35")
      } else {
        barCOLS <- scales::hue_pal()(num_groups)
        dotCOLS <- scales::hue_pal()(num_groups)
      }

      forest_plot <- ggplot(forest_df, aes(x=clust_id, y=OR, ymin=CI_lower, ymax=CI_upper, col=group, fill=group)) +
        #specify position here
        geom_linerange(linewidth=3, position=position_dodge(width = 0.5)) +
        geom_hline(yintercept=1, lty=2) +
        #specify position here too
        geom_point(size=3, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
        scale_fill_manual(values=barCOLS)+
        scale_color_manual(values=dotCOLS)+
        scale_x_discrete(name="Cluster ID") +
        scale_y_continuous(name="Odds ratio", breaks=sort(c(seq(0, 10, length.out=5), 1))) +
        coord_flip(ylim=c(0,10)) +
        theme_minimal() +
        theme(axis.line = element_line(color = "black", linewidth = 0.5)) +
        labs(fill="Cohort", color="Cohort")

      # store in reactive values
      values_network_netA$clust_forest_df <- forest_df
      values_network_netA$clust_forest_plot <- forest_plot

      tend_check <- TRUE
    }

    # check clust method have been utilized
    if(!clust_select %in% colnames(clust_node)){
      showNotification(paste0("Cannot find the ", clust_select," clusters."))
      resetLoadingButton("make_clust_network")
      req(FALSE)
    }

    # visual network
    lnodes <- data.frame(label=unique(clust_node[[clust_select]]), shape=c("dot"), color=col_vector[1:length(unique(clust_node[[clust_select]]))])
    color_v <- lnodes$color[match(clust_node[[clust_select]], lnodes$label)]
    node_to_plot <- clust_node; node_to_plot$color <- color_v
    lnodes$color[lnodes$label == "unclustered"] <- "gray" # change the color of unclustered to gray
    node_to_plot$color[node_to_plot[[clust_select]] == "unclustered"] <- "gray"
    visnet_clust <- visNetwork(node_to_plot, clust_edge, height="800px") %>%
      visIgraphLayout(layout="layout_with_fr") %>%
      visInteraction(navigationButtons = TRUE) %>%
      visNodes(
        shape = "dot",
        color = list(
          highlight = "yellow"),
          font = list(size = 20)
      ) %>%
      visEdges(
        shadow = FALSE,
        color = list(color = "#0085AF", highlight = "#C62F4B")
      ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visLegend(useGroups = F, width=0.2, addNodes = lnodes, position = "right")

    # store in reactive values
    values_network_netA$clust_graph <- visnet_clust

    # keep target label in node table
    clust_node <- clust_node[, c(1:7, which(names(clust_node) == clust_select))]

    # stored visualize edge and node
    values_network_netA$clust_edge_vis <- clust_edge
    values_network_netA$clust_node_vis <- clust_node

    # render output
    if(tend_check == TRUE) {
      output$forest_plot_ui <- renderUI({
        tagList(
          HTML("<h4><b>Correlation between clusters and cohorts</b></h4>"),
          plotOutput("clust_forest_plot")
        )
      })

      output$clust_forest_plot <- renderPlot({
        values_network_netA$clust_forest_plot
      })

      output$forest_df_ui <- renderUI({
        tagList(
          div(style = "width: 60%;", DT::dataTableOutput("forest_df")),
          HTML("<br>"),
          HTML("<ul style='text-align: left; color: gray;'><li> The blank values are the <b>Inf</b> due to zero cell in matrix. </li></ul>")
          )
      })

      output$forest_df <- DT::renderDT({
        datatable(values_network_netA$clust_forest_df, options=list(scrollX=T)) %>%
                    formatSignif(columns = c('OR', 'p_val', 'CI_lower', 'CI_upper'), digits = 3)
      })
    } else {
      output$forest_plot_ui <- renderUI({ NULL })
      output$clust_forest_plot <- renderPlot({ NULL })
      output$forest_df_ui <- renderUI({ NULL })
      output$forest_df <- DT::renderDT({ NULL })
    }

    output$edge_clust_ui <- renderUI({
      tagList(
        HTML("<h4><b>Edge</b></h4>"),
        DT::dataTableOutput("edge_clust_table")
      )
    })

    output$node_clust_ui <- renderUI({
      tagList(
        HTML("<h4><b>Node</b></h4>"),
        DT::dataTableOutput("node_clust_table")
      )
    })

    output$edge_clust_table <- DT::renderDT({
      datatable(clust_edge, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(6,7))))) %>%
                  formatSignif(columns = c('weight', 'p', 'p_adj'), digits = 3)
    })

    output$node_clust_table <- DT::renderDT({
      datatable(clust_node, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,7))))) %>%
                  formatSignif(columns = c('raw_value'), digits = 3)
    })

    output$clust_graph <- renderVisNetwork({
      visnet_clust
    })

    resetLoadingButton("make_clust_network")
  }) # end of make_clust_network


  # save clustered network
  observeEvent(input$clust_graph_save, {
    clust_method <- values_network_netA$clust_method_vis
    target_label <- values_network_netA$target_label

    # network analysis directory
    netA_dir <- file.path(working_dir, "network_analysis")
    if(!dir.exists(netA_dir)){
      dir.create(netA_dir)
    }

    save_dir <- file.path(netA_dir, paste(c(target_label, clust_method, "network"), collapse="_"))
    if(!dir.exists(save_dir)) {
      dir.create(save_dir)
      if(target_label == "All") {
        ggsave(filename=file.path(save_dir, paste0(target_label, "_correlation.png")),
                plot=values_network_netA$clust_forest_plot, bg="white", width=10, height=8, units="in")
        write.csv(values_network_netA$clust_forest_df, file=file.path(save_dir, paste0(target_label, "_correlation.csv")), row.names=F)
      }

      write.csv(values_network_netA$clust_node_vis, file=file.path(save_dir, paste(c(target_label, clust_method, "node.csv"), collapse="_")), row.names=F)
      write.csv(values_network_netA$clust_edge_vis, file=file.path(save_dir, paste(c(target_label, clust_method, "edge.csv"), collapse="_")), row.names=F)
      visSave(values_network_netA$clust_graph, file=file.path(save_dir, paste0(target_label, "_", clust_method, ".html")))
      showNotification("Results are saved in network_analysis successfully!")
    } else {
      showNotification(paste0("Directory of ", target_label, "_", clust_method, " already exists. Cannot save the result."))
      req(FALSE)
    }
  }) # end of clust_graph_save


  # ===== main Event 7 (Clustered pathway analysis) =====
  output$path_clust_bar_size_input <- renderUI({
    if (is.null(values_network_netA$netA_clust_path_plot)) {
      return(NULL)
    }
    tagList(
      numericInput("path_clust_bar_width", "Width", value = 1100, min = 1, max = 2000),
      numericInput("path_clust_bar_height", "Height", value = 500, min = 1, max = 2000)
    )
  })

  observeEvent(input$netA_path_button, {
    clust_label <- input$label_netA_path
    values_network_netA$clust_path_label <- clust_label
    clust_method <- values_network_netA$clust_method_vis
    clust_id <- input$netA_path_clust_id
    values_network_netA$clust_id <- clust_id
    clust_p <- input$netA_path_p

    # check whether the cohort clustering task has been performed
    if(!file.exists(file.path(working_dir, "network_analysis", paste0(clust_label, "_cluster_overview")))) {
      showNotification("Cannot find the current cohort clustered results. Please save the Network Cluster results first.")
      resetLoadingButton("netA_path_button")
      req(FALSE)
    }

    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "clust_pathway_analysis.py"),
                          working_dir, clust_label, clust_id, clust_p, threads))
    if( exit_code != 0){
      showNotification("Error occurred in Clustered pathway analysis. Please check the log in terminal.")
      resetLoadingButton("netA_path_button")
      req(FALSE)
    }

    # read table
    clust_pathway_res <- read.csv(file.path(tmp_dir, paste(clust_label, clust_id, "pathway.csv", sep="_")))
    # Reactive values
    values_network_netA$clust_id_pathway <- clust_pathway_res

    if(nrow(clust_pathway_res) != 0) { # check whether the table is empty

      # bar plots
      # node
      node_plot <- data.frame(Pathway=values_network_netA$clust_id_pathway$annotation,
                              OR=values_network_netA$clust_id_pathway$step1_stats,
                              p=values_network_netA$clust_id_pathway$step1_p_adj,
                              Type="Node")

      # edge
      edge_plot <- data.frame(Pathway=values_network_netA$clust_id_pathway$annotation,
                              OR=values_network_netA$clust_id_pathway$step2_stats,
                              p=values_network_netA$clust_id_pathway$step2_p_adj,
                              Type="Edge")

      combined_plot <- rbind(node_plot, edge_plot)
      combined_plot <- combined_plot[order(combined_plot$Pathway, combined_plot$OR), ]
      p <- ggplot(combined_plot, aes(x = reorder(Pathway, OR), y = OR, fill = p, group = Type)) +
        scale_fill_gradient(low = "red", high = "pink") +
        coord_flip() +
        labs(x = "Pathway", y = "Odds Ratio", fill = "p-value") +
        theme_minimal() +
        theme(axis.line = element_line(color = "black", linewidth = 0.5),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10))

      for (i in 1:(length(node_plot$Pathway))) {
        p <- p + annotate('rect', xmin = i-0.5, xmax = i+0.5, ymin = -Inf, ymax = Inf,
                          fill = ifelse(i %% 2 == 0, 'white', 'gray95'))}

      p <- p + geom_bar(stat = "identity", position = position_dodge(width = 0.6), width=0.6, colour="black") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      values_network_netA$netA_clust_path_plot <- p # stored in reactive values

    } else {
      values_network_netA$netA_clust_path_plot <- NULL
    }

    # ouput
    output$netA_path_title <- renderText({
      paste("Pathway analysis of", clust_id, "in", clust_label)
    })

    output$netA_pathway_note <- renderText({
      paste0("( adjusted p_value threshold : ", clust_p, " )")
    })

    output$path_clust_table <- DT::renderDT({
      datatable(values_network_netA$clust_id_pathway, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,6))))) %>%
                  formatSignif(columns = c('step1_stats', 'step1_p_adj', 'step2_stats', 'step2_p_adj'), digits = 3)
    })

    output$path_clust_plot <- renderPlot(
      values_network_netA$netA_clust_path_plot,
      width = function() input$path_clust_bar_width,
      height = function() input$path_clust_bar_height,
      res = 100
    )

    output$netA_path_bar_note <- renderUI({
      if (is.null(values_network_netA$netA_clust_path_plot)) {
        return(NULL)
      }

      if(!is.null(values_network_netA$netA_clust_path_plot)) {
        return(HTML("<ul>
          <li>The first bar: odds ratio of node.</li>
          <li>The second bar: odds ratio of edge.</li>
        </ul>"))
      }
    })

    resetLoadingButton("netA_path_button")
  }) # end of path_clust_button


  # save button
  observeEvent(input$netA_path_save, {
    target_label <- values_network_netA$clust_path_label
    target_id <- values_network_netA$clust_id
    netA_dir <- file.path(working_dir, "network_analysis")
    if(!dir.exists(netA_dir)){
      dir.create(netA_dir)
    }

    res_dir <- file.path(netA_dir, paste(target_label, target_id, "pathway", sep="_"))
    if(!dir.exists(res_dir)){
      dir.create(res_dir)
      file.copy(from=file.path(tmp_dir, paste(target_label, target_id, "pathway.csv", sep="_")), to=res_dir)
      ggsave(filename=file.path(res_dir, paste(target_label, target_id, "pathway_plot.png", sep="_")),
             plot=values_network_netA$netA_clust_path_plot, bg="white",
             width=input$path_clust_bar_width, height=input$path_clust_bar_height, units="px", dpi=100)
      showNotification("Result is saved in network_analysis successfully!")
    } else {
      showNotification(paste0("Directory of ", res_dir, " already exists. Cannot save the result."))
      req(FALSE)
    }
  }) # end of netA_path_save

  # netA_pathway_visualize
  observeEvent(input$netA_path_vis, {
    target_pathway <- input$netA_path_map_id
    values_network_netA$clust_vis_mapid <- target_pathway
    target_label <- values_network_netA$clust_path_label
    clust_id <- values_network_netA$clust_id

    # check whether the cohort clustering task has been performed
    if(!file.exists(file.path(working_dir, "network_analysis", paste0(target_label, "_cluster_overview")))) {
      showNotification("Cannot find the current cohort clustered results. Please save the Network Cluster results first.")
      resetLoadingButton("netA_path_button")
      req(FALSE)
    }

    # check network filter results
    if(!dir.exists(file.path(working_dir, "network_filter", target_label))) {
      showNotification(
        paste0("Cannot find the filtered sub-network result of ", target_label, ". Please run network filter first."))
        resetLoadingButton("pathway_net_button")
        req(FALSE)
    }

    # call py script
    exit_code <- system(paste("python", file.path(script_dir, "clust_path_net.py"),
                          working_dir, target_pathway, target_label, clust_id))
    if(exit_code != 0){
      showNotification("Error occurred in clustered pathway sub-network visualization. Please check the log in terminal.")
      resetLoadingButton("netA_path_vis")
      req(FALSE)
    }

    # read table
    path_edge <- read.csv(file.path(tmp_dir, paste(target_label, clust_id, target_pathway, "edge.csv", sep="_")))
    path_node <- read.csv(file.path(tmp_dir, paste(target_label, clust_id, target_pathway, "node.csv", sep="_")))
    # rds file
    ko_map <- readRDS(file.path(script_dir, "whole_komap", "total_KO_pathway.rds"))
    annotated_name <- ko_map[ko_map$mapid == target_pathway, "pathway"][1]

    # stored in reactive values
    values_network_netA$clust_path_edge <- path_edge
    values_network_netA$clust_path_node <- path_node

    # visualize network
    if(!"group" %in% colnames(path_node)) {
      visnet <- visNetwork(path_edge, path_node, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            background = "#0085AF",
            border = "#013848",
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

    } else {
      # deal with colored by groups
      lnodes <- data.frame(label=unique(path_node$group), shape=c("dot"), color=col_vector[1:length(unique(path_node$group))])
      color_v <- lnodes$color[match(path_node$group, lnodes$label)]
      node_to_plot <- path_node; node_to_plot$color <- color_v    # assgin color column
      visnet <- visNetwork(node_to_plot, path_edge, height="800px") %>%
        visIgraphLayout(layout="layout_with_fr") %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(
          shape = "dot",
          color = list(
            highlight = "yellow"),
            font = list(size = 20)
        ) %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#0085AF", highlight = "#C62F4B")
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLegend(useGroups = F, width=0.2, addNodes = lnodes, position = "right")
    }

    values_network_netA$clust_path_visnet <- visnet  # stored in reactive values

    # output
    output$netA_path_vis_title <- renderUI({
      tagList(
        HTML(paste0("The selected pathway is : <b>", annotated_name,"</b><br>")),
        HTML(paste0("The selected cohort is : <b>", target_label, "</b>"))
      )
    })

    output$netA_path_vis_edgeUI <- renderUI({
      tagList(
        HTML("<h4><b>Edge</b></h4>"),
        DT::dataTableOutput("netA_path_edge")
      )
    })

    output$netA_path_vis_nodeUI <- renderUI({
      tagList(
        HTML("<h4><b>Node</b></h4>"),
        DT::dataTableOutput("netA_path_node")
      )
    })

    output$netA_path_edge <- DT::renderDT({
      datatable(values_network_netA$clust_path_edge, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(6,7))))) %>%
                  formatSignif(columns = c('weight', 'p', 'p_adj'), digits = 3)
    })

    output$netA_path_node <- DT::renderDT({
      datatable(values_network_netA$clust_path_node, options=list(scrollX=T, columnDefs=list(list(visible=FALSE, targets=c(3,7))))) %>%
                  formatSignif(columns = c('raw_value'), digits = 3)
    })

    output$netA_path_graph <- renderVisNetwork({
      values_network_netA$clust_path_visnet })

    resetLoadingButton("netA_path_vis")
  }) # end of netA_path_vis

  # save button
  observeEvent(input$netA_path_graph_save, {
    target_label <- values_network_netA$clust_path_label
    target_pathway <- values_network_netA$clust_vis_mapid
    cluster_id <- values_network_netA$clust_id
    netA_dir <- file.path(working_dir, "network_analysis")
    if(!dir.exists(netA_dir)){
      dir.create(netA_dir)
    }
    save_dir <- file.path(netA_dir, paste0(target_label, "_", cluster_id, "_", target_pathway))

    if (!dir.exists(save_dir)){
      dir.create(save_dir)
      visSave(values_network_netA$clust_path_visnet, file = file.path(save_dir, paste0(target_label, "_", cluster_id, "_", target_pathway, ".html")))
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", cluster_id, "_", target_pathway, "_edge.csv")), to=save_dir)
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", cluster_id, "_", target_pathway, "_node.csv")), to=save_dir)
      file.copy(from=file.path(tmp_dir, paste0(target_label, "_", cluster_id, "_", target_pathway, "_koset.txt")), to=save_dir)
      showNotification("Results are saved in Clustered pathway analysis visualized network successfully!")
    } else {
      showNotification(paste0("Directory of ", target_pathway, " already exists. Cannot save the result."))
      req(FALSE)
    }
  }) # end of netA_path_graph_save


  # Delete the tmp_files when the session ends
  session$onSessionEnded(function() {
    if (file.exists(file.path(tmp_dir))) {
      unlink(tmp_dir, recursive = TRUE)
      stopApp()
    }
  })

}


# Run the application
shinyApp(ui=ui, server=server, onStart=onstart)

