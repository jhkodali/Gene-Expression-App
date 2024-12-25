library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library(colourpicker)
library(DT)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
library(ggbeeswarm)

ui <- fluidPage(
  
  titlePanel(h1("Gene Expression Explorer",
                h3("Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"))),
  
  # navigation bar
  navbarPage("Data Exploration",
             
             # samples tab
             tabPanel("Samples",
                      sidebarLayout(
                        sidebarPanel(
                          # accept = ".csv" only allows for CSV uploads
                          fileInput("samples", "Upload CSV file", accept = ".csv"),
                          selectInput("samples_variable", "Select variable to plot",
                                      c("Post-mortem Interval" = "post_mortem_interval",
                                        "Age of Death" = "age_of_death",
                                        "RNA Integrity Number" = "RNA_integrity_number",
                                        "mRNAseq Reads" = "mRNAseq_reads"))
                        ),
                        mainPanel(
                          h4("Exploring Sample Data"),
                          p("In this section, you can upload the GSE64810_samples_input.csv file."),
                          p("This file includes various metadata information such as age of death and diagnosis."),
                          tabsetPanel(
                            tabPanel("Summary", 
                                     p("Here is a summary of different annotations found in the metadata"),
                                     tableOutput("summary_table")),
                            tabPanel("Table", 
                                     p("Here is the table of the sample matrix. You can sort columns and search for text"),
                                     DTOutput("sample_table")),
                            tabPanel("Plot", 
                                     p("You can choose a variable from the dropdown menu to plot."),
                                     plotOutput("sample_plot"))
                          )
                        )
                      )),
             
             # counts tab
             tabPanel("Counts", 
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("counts", "Upload CSV file", accept = ".csv"),
                          sliderInput("filter1", "Select gene percentile variance treshold:",
                                      min = 0, max = 100, value = 50),
                          sliderInput("filter2", "Select minimum number of non-zero samples to filter genes:",
                                      # total number of samples is 69
                                      min = 0, max = 69, value = 0),
                          checkboxInput("log_transform", "Log Transform Counts for Heatmap", value = TRUE)
                        ),
                        mainPanel(
                          h4("Exploring Counts Matrix"),
                          p("Here, you can upload the GSE64810_counts_input.csv to explore the counts matrix"),
                          p("You can use the slider to choose thresholds for the gene percentile of variance and number of non-zero samples"),
                          p("Here, you can find different visualizations based on filtering metrics."),
                          tabsetPanel(
                            tabPanel("Summary", tableOutput("counts_summary")),
                            tabPanel("Scatterplot", plotOutput("counts_scatter")),
                            tabPanel("Heatmap", plotOutput("counts_heat")),
                            tabPanel("PCA", 
                                     p("You can choose which PC is plotted for each axis."),
                                     sidebarLayout(
                                       # max 69 PCs
                                       selectInput("pc_x", "Select PC for the X axis", choices = 1:69, selected = 1),
                                       selectInput("pc_y", "Select PC for the Y axis", choices = 1:69, selected = 2)
                                     ),
                                     mainPanel(
                                       plotOutput("counts_PCA")
                                     ))
                            
                          )
                        )
                      )),
             
             # differential expression tab
             tabPanel("Differential Expression",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("deseq", "Upload CSV file", accept = ".csv"),
                          radioButtons("dist_x", "Choose the column for the x-axis",
                                       c("baseMean" = "baseMean",
                                         "log2FoldChange" = "log2FoldChange",
                                         "lfcSE" = "lfcSE",
                                         "stat" = "stat",
                                         "pvalue" = "pvalue",
                                         "padj" = "padj"),
                                       selected = "log2FoldChange"),
                          radioButtons("dist_y", "Choose the column for the y-axis",
                                       c("baseMean" = "baseMean",
                                         "log2FoldChange" = "log2FoldChange",
                                         "lfcSE" = "lfcSE",
                                         "stat" = "stat",
                                         "pvalue" = "pvalue",
                                         "padj" = "padj"),
                                       selected = "padj"),
                          colourInput("col_1", "Base point color", "maroon"),
                          colourInput("col_2", "Highlight point color", "cadetblue3"),
                          sliderInput("obs", "Select the magnitude of the p adjusted coloring:",
                                      # lowest padj was 10^-34
                                      min = -35, max = 0, value = -2),
                          actionButton("submitbutton", "Plot it!", class = "btn-success", icon = icon("wand-sparkles"), width = '100%')
                        ),
                        mainPanel(
                          h4("Exploring Differential Expression Results"),
                          p("Here, you can upload the GSE64810_DESeq2_input.csv file to explore the results from DESEq analysis."),
                          tabsetPanel(
                            tabPanel("Table", 
                                     p("This table allows you to search through the results of DESeq analysis"),
                                     DTOutput("deseq_table")),
                            tabPanel("Plot", 
                                     p("Here you can choose the parameters for the x and y columns and plot a volcano plot."),
                                     p("You can choose the magnitude of padj value as a filtering threshold."),
                                     plotOutput("deseq_plot"))
                          )
                        )
                      )),
             
             # individual gene expression tab
             tabPanel("Individual Gene Expression", 
                      sidebarLayout(
                        sidebarPanel(
                          # no need for csv inputs, should already be loaded
                          # could easily add if user prefers to use the tab individually
                          # only one categorical variable in sample matrix
                          selectInput("gene_exp", "Choose a categorical variable",
                                      c("Diagnosis" ="diagnosis")),
                          selectizeInput("select_gene", "Select Gene Symbol:", choices = NULL, multiple = FALSE),
                          radioButtons("type_plot", "Select plot:",
                                       c("Bar Plot" = "barplot",
                                         "Boxplot" = "boxplot",
                                         "Violin Plot" = "violin",
                                         "Beeswarm Plot" = "beeswarm")),
                          actionButton("submitbutton2", "Plot it!", class = "btn-success", icon = icon("wand-sparkles"), width = '100%')
                        ),
                        mainPanel(
                          h4("Exploring Individual Gene Expression"),
                          p("This page allows you to visualize the gene expression of an individual gene"),
                          p("This page uses the input files from the Samples and Counts files."),
                          p("You can select a category to group gene counts by (in this case, there is only one category: diagnosis."),
                          p("You can choose a gene by typing into the search box or finding one in dropdown menu."),
                          p("You can also choose which type of plot."),
                          plotOutput("exp_plot")
                        )
                      )
             )
  ))

# setting max upload size (counts csv went over default max)
options(shiny.maxRequestSize = 30 * 1024^2)

server <- function(input, output, session) {
  
  # samples tab
  # read in samples csv
  sample_data <- reactive({
    req(input$samples)
    data <- read_csv(input$samples$datapath)
    
    # error if input file is empty
    if (nrow(data) == 0) {
      showNotification("Error: The CSV file is empty.", type = "error")
      return(NULL)  
    }
    
    return(data)
  })
  
  # function to summarize columns from samples data
  get_sum_table <- function(df) {
    
    sum_table <- tibble(
      "Column Name" = c(
        "Diagnosis", 
        "Post-mortem Interval", 
        "Age of Death", 
        "RNA Integrity Number",
        "mRNA-seq Reads"),
      
      "Type" = c(
        paste0(class(df$diagnosis)), 
        paste0(class(df$post_mortem_interval)), 
        paste0(class(df$age_of_death)), 
        paste0(class(df$RNA_integrity_number)), 
        paste0(class(df$mRNAseq_reads))
      ),
      
      "Mean (sd) or Distinct Values" = c(
        paste(unique(df$diagnosis), collapse = ", "),
        paste0(round(mean(df$post_mortem_interval, na.rm = TRUE), 2), " (+/- ",
               round(sd(df$post_mortem_interval, na.rm = TRUE), 2), ")"),
        paste0(round(mean(df$age_of_death, na.rm = TRUE), 2), " (+/- ",
               round(sd(df$age_of_death, na.rm = TRUE), 2), ")"),
        paste0(round(mean(df$RNA_integrity_number, na.rm = TRUE), 2), " (+/- ",
               round(sd(df$RNA_integrity_number, na.rm = TRUE), 2), ")"),
        paste0(round(mean(df$mRNAseq_reads, na.rm = TRUE), 2), " (+/- ",
               round(sd(df$mRNAseq_reads, na.rm = TRUE), 2), ")")
      )
    )
  }
  
  # samples: summary tab output
  output$summary_table <- renderTable({
    sum_df <- sample_data()
    table <- get_sum_table(sum_df)
    return(table)
  })
  
  
  # samples: table tab output
  output$sample_table <- renderDT({
    datatable(sample_data(), options = list(
      pageLength = 15
    ))
  })
  
  # function to plot samples violin plots
  sample_violin_plot <- function(data, var) {
    ggplot(data, aes(x = diagnosis, y = .data[[var]], fill = diagnosis)) +
      geom_violin(trim =FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
      labs(title = paste(tools::toTitleCase(gsub("_", " ", var)), "Grouped by Diagnosis"),
           x = "Diagnosis",
           y = gsub("_", " ", var)) +
      theme_minimal() +
      # plot labels were small so making larger
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(angle = 90, size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15)
      )
  }
  
  # samples: plot tab output
  output$sample_plot <- renderPlot({
    sample_data <- sample_data()
    sample_violin <- sample_violin_plot(sample_data, input$samples_variable)
    return(sample_violin)
  })
  
  
  
  # counts tab
  # read in counts csv
  counts_data <- reactive({
    req(input$counts)
    data <- read_csv(input$counts$datapath) %>%
      rename_with(.cols = 1, ~"genes")
    
    # error if input file is empty
    if (nrow(data) == 0) {
      showNotification("Error: The CSV file is empty.", type = "error")
      return(NULL)  
    }
    
    return(data)
  })
  
  # function to filter the counts
  filter_counts <- function(data, filter1, filter2) {
    gene_names <- data$genes
    data_numeric <-data[, sapply(data, is.numeric)]
    
    gene_var <- apply(data_numeric, 1, var, na.rm = TRUE)
    var_threshold <- quantile(gene_var, probs = filter1 / 100, na.rm = TRUE)
    filter_var <- gene_var >= var_threshold
    
    non_zero <- apply(data_numeric, 1, function(gene) sum(gene != 0, na.rm = TRUE))
    filter_zero <- non_zero >= filter2
    
    final_filtered <- filter_var & filter_zero
    filtered_data <- data[final_filtered, , drop = FALSE]
    
    filtered_data$genes <- gene_names[final_filtered]
    
    return(filtered_data)
    
  }
  
  # reactive function to to summarize filtering
  # initially did this when rendering table, but table wasn't updating when changing inputs
  summarize_filter <- reactive({
    
    counts_pre <- counts_data()
    counts_post <- filter_counts(counts_pre, input$filter1, input$filter2)
    
    num_samples <- ncol(counts_pre) - 1
    total_num_genes <- nrow(counts_pre)
    filter_num_genes <- nrow(counts_post)
    percent_pass <- (filter_num_genes / total_num_genes) * 100
    not_pass_genes <- total_num_genes - filter_num_genes
    percent_not_pass <- (not_pass_genes / total_num_genes) * 100
    
    sum_counts_table <- tibble(
      Measure = c("Number of Samples", 
                  "Total Number of Genes",
                  "Number of genes passing filter",
                  "% of genes passing filter",
                  "Number of genes not passing",
                  "% of genes not passing"
      ),
      Value = c(num_samples,
                total_num_genes,
                filter_num_genes,
                paste0(percent_pass, "%"),
                not_pass_genes,
                paste0(percent_not_pass, "%")
      )
    )
    return(sum_counts_table)
  })
  
  # counts: summary tab output
  output$counts_summary <- renderTable({
    summarize_filter()
  })
  
  
  # function to plot median count vs variance
  median_vs_variance <- function(counts_pre, counts_post) {
    data_numeric <- counts_pre[, sapply(counts_pre, is.numeric)]
    median_counts <- apply(data_numeric, 1, median, na.rm =TRUE)
    
    gene_var <- apply(data_numeric, 1, var, na.rm = TRUE)
    
    plot_data <- tibble::tibble(median = median_counts, variance = gene_var)
    plot_data$rank <- rank(plot_data$median)
    plot_data$filtered <- factor(counts_pre$genes %in% counts_post$genes,
                                 levels = c(FALSE, TRUE),
                                 labels = c("Fail", "Pass"))
    
    mv_plot <- ggplot(plot_data, aes(x = rank, y = variance)) +
      geom_point(aes(color = filtered), alpha = 0.5) +
      labs(x = "Median Rank",
           y = "Log Variance",
           title = "Median Rank vs Variance") +
      scale_y_log10() +
      scale_color_manual(name = "Gene Filter", values = c("Fail" = "lightgray", 
                                                          "Pass" = "deepskyblue4")) +
      geom_smooth(color = "black") +
      theme_minimal() +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16)) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    
    return(mv_plot)
  }
  
  # function to plot median count vs number of zeros
  median_vs_zeros <- function(counts_pre, counts_post) {
    data_numeric <- counts_pre[, sapply(counts_pre, is.numeric)]
    median_counts <- apply(data_numeric, 1, median, na.rm = TRUE)
    num_zeros <- apply(data_numeric, 1, function(gene) sum(gene == 0, na.rm = TRUE))
    
    plot_data <- tibble::tibble(median = median_counts, zeros = num_zeros)
    plot_data$rank <- rank(plot_data$median)
    plot_data$filtered <- factor(counts_pre$genes %in% counts_post$genes,
                                 levels = c(FALSE, TRUE),
                                 labels = c("Fail", "Pass"))
    
    mz_plot <- ggplot(plot_data, aes(x = rank, y = zeros)) +
      geom_point(aes(color = filtered), alpha = 0.5) +
      labs(x = "Median Rank",
           y = "Number of Zeros",
           title = "Median Rank vs Number of Zeros") +
      scale_color_manual(name = "Gene Filter", values = c("Fail" = "lightgray",
                                                          "Pass" = "deepskyblue4")) +
      theme_minimal() +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16)) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    return(mz_plot)
  }
  
  # counts: scatterplot tab output
  output$counts_scatter <- renderPlot({
    counts_pre <- counts_data()
    counts_post <- filter_counts(counts_pre, input$filter1, input$filter2)
    
    mv <- median_vs_variance(counts_pre, counts_post)
    mz <- median_vs_zeros(counts_pre, counts_post)
    
    grid.arrange(mv, mz, nrow = 2)
  })
  
  
  # function to plot heatmap
  heatmap_counts <- function(filtered_counts, log_transform) {
    data <- filtered_counts[, -1]
    
    if (log_transform) {
      data <- log1p(data)
    }
    
    data <- as.matrix(data)
    
    heat <- pheatmap(
      mat = data,
      color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      legend = TRUE,
      main = "Clustered Heatmap",
      fontsize = 10
    )
    
    return(heat)
    
  }
  
  # counts tab: heatmap output
  output$counts_heat <- renderPlot({
    counts_pre <- counts_data()
    counts_post <- filter_counts(counts_pre, input$filter1, input$filter2)
    
    heatmap_counts(counts_post, input$log_transform)
  }, width = 800, height = 1000)
  
  # function to compute and plot PCA
  pca_plot <- function(counts, pc_x, pc_y) {
    counts_mat <- t(as.matrix(counts[, -1]))
    pca <- prcomp(counts_mat, center = TRUE, scale = TRUE)
    
    pca_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100
    pca_df <- as.data.frame(pca$x)
    
    # pc are strings so need to convert into integers to access index of pca_var
    num_x <- as.integer(pc_x)
    num_y <- as.integer(pc_y)
    var_x <- pca_var[num_x]
    var_y <- pca_var[num_y]
    
    ggplot(pca_df, 
           aes_string(x = paste0("PC", pc_x), 
                      y = paste0("PC", pc_y))) +
      geom_point(size = 3, color = "deepskyblue3") +
      labs(title = paste("PC", num_x, "vs PC", num_y),
           x = paste("PC", num_x, "(", round(var_x, 2), "% Variance)"),
           y = paste("PC", num_y, "(", round(var_y, 2), "% Variance)")
      ) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15))
  }
  
  # counts tab: PCA output
  output$counts_PCA <- renderPlot({
    counts <- counts_data()
    pca_plot(counts, input$pc_x, input$pc_y)
  })
  
  
  # differential expression tab
  # read in DESeq results
  DESeq_res <- reactive({
    req(input$deseq)
    data <- read_csv(input$deseq$datapath) %>%
      rename_with(.cols = 1, ~"genes")
    
    # error if input file is empty
    if (nrow(data) == 0) {
      showNotification("Error: The CSV file is empty.", type = "error")
      return(NULL)  
    }
    
    return(data)
  })
  
  # function to create filterable table for DESEq results
  get_deseq_table <- function(deseq_data, slider) {
    filtered <- deseq_data %>%
      filter(padj < 10^(slider))
    filtered$padj <- formatC(filtered$padj, format = "e", digits = 5)
    filtered$pvalue <- formatC(filtered$pvalue, format = "e", digits = 5)
    return(filtered)
  }
  
  # differential expression: table output
  output$deseq_table <- renderDT({
    deseq_res <- DESeq_res()
    deseq_tab <- get_deseq_table(deseq_res, input$obs)
    
    datatable(deseq_tab, options = list(
      pageLength = 15
    ))
  })
  
  # function to plot volcano plot of user input columns in DESeq results
  get_deseq_plot <- function(deseq_data, x_name, y_name, slider, color1, color2) {
    deseq_data <- deseq_data %>%
      mutate(
        # sets status if less than user input threshold to be true, otherwise false
        volc_plot_status = case_when(
          !!sym(y_name) < 10^(slider) ~ "TRUE",
          !!sym(y_name) > 10^(slider) ~ "FALSE"
        ),
        # calculating -log10 of y variable selected
        log_y = -log10(!!sym(y_name))
      )
    plot <- ggplot(deseq_data, aes(x = !!sym(x_name), y = log_y, color = volc_plot_status)) +
      geom_point() +
      labs(
        title = paste(x_name, "vs - log10(", y_name, ")"),
        x = x_name,
        y = paste("-log10(", y_name, ")")
      ) +
      theme_minimal() +
      scale_color_manual(values = c("TRUE" = color2, "FALSE" = color1)) +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16)) +
      guides(color = guide_legend(override.aes = list(size=5)))
    
    return(plot)
  }
  
  # differential expression: plot output
  output$deseq_plot <- renderPlot({
    deseq_res <- DESeq_res()
    get_deseq_plot(deseq_res, 
                   input$dist_x, 
                   input$dist_y, 
                   input$obs, 
                   input$col_1, 
                   input$col_2)
  })
  
  
  # individual gene expression tab
  # input csvs already uploaded in samples and counts tabs
  
  # updates SelectizeInput so it just has a subset of genes in drop down list.
  # when user types in input box, genes with those characters will pop up
  observe({
    deseq_res <- DESeq_res()
    updateSelectizeInput(session, "select_gene", choices = deseq_res$symbol, server = TRUE)
  })
  
  # reactive function to prepare data when button is clicked
  gene_exp_data <- eventReactive(input$submitbutton2, {
    req(input$samples, input$counts, input$deseq)
    
    deseq_res <- DESeq_res()
    counts <- counts_data()
    samples <- sample_data()
    
    # adding gene symbols to counts data (to be able to search by gene instead of ensembl id)
    counts_symbol <- counts %>%
      left_join(deseq_res %>%
                  select(genes, symbol),
                by = "genes") %>%
      relocate(symbol, .after = genes)
    
    # filter for selected gene row
    counts_filtered <- counts_symbol %>%
      filter(symbol == input$select_gene)
    
    # pivot longer for plotting
    counts_long <- counts_filtered %>%
      pivot_longer(cols = -c(genes, symbol), 
                   names_to = "title",
                   values_to = "count")
    
    # combining the long format counts with the samples matrix
    merged_data <- counts_long %>%
      inner_join(samples, by = "title")
    
    return(merged_data)
    
  })
  
  # function to plot bar plot
  gene_exp_barplot <- function(data, category, gene) {
    ggplot(data, aes(x = !!sym(category), y = count, fill = !!sym(category))) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = paste("Individual Expression of", gene, "Grouped by Diagnosis"), 
           x = "Diagnosis",
           y = "Counts") +
      scale_fill_manual(values = c("Neurologically normal" = "cadetblue3", 
                                   "Huntington's Disease" = "maroon")) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15))
  }
  
  # function to plot boxplot
  gene_exp_boxplot <- function(data, category, gene) {
    ggplot(data, aes(x = !!sym(category), y = count, fill = !!sym(category))) +
      geom_boxplot() +
      labs(title = paste("Individual Expression of", gene, "Grouped by Diagnosis"),
           x = "Diagnosis",
           y = "Counts") +
      scale_fill_manual(values = c("Neurologically normal" = "cadetblue3", 
                                   "Huntington's Disease" = "maroon")) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15))
  }
  
  # function to plot violin plot
  gene_exp_violin <- function(data, category, gene) {
    ggplot(data, aes(x = !!sym(category), y = count, fill = !!sym(category))) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      labs(title = paste("Individual Expression of", gene, "Grouped by Diagnosis"),
           x = "Diagnosis",
           y = "Counts") +
      scale_fill_manual(values = c("Neurologically normal" = "cadetblue3", 
                                   "Huntington's Disease" = "maroon")) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15))
  }
  
  # function to plot beeswarm plot
  gene_exp_beeswarm <- function(data, category, gene) {
    ggplot(data, aes(x = !!sym(category), y = count, color = !!sym(category))) +
      geom_beeswarm() +
      labs(title = paste("Individual Expression of", gene, "Grouped by Diagnosis"),
           x = "Diagnosis",
           y = "Counts") +
      scale_color_manual(values = c("Neurologically normal" = "cadetblue3", 
                                    "Huntington's Disease" = "maroon")) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15))
  }
  # individual gene expression: plot output
  output$exp_plot <- renderPlot({
    gene_exp <- gene_exp_data()
    
    # selecting which plot function to run based on user input
    plot <- switch(input$type_plot,
                   "barplot" = gene_exp_barplot(gene_exp, input$gene_exp, input$select_gene),
                   "boxplot" = gene_exp_boxplot(gene_exp, input$gene_exp, input$select_gene),
                   "violin" = gene_exp_violin(gene_exp, input$gene_exp, input$select_gene),
                   "beeswarm" = gene_exp_beeswarm(gene_exp, input$gene_exp, input$select_gene))
    return(plot)
  })
  
  
  
  # stopping session after closing app to prevent RStudio from crashing
  session$onSessionEnded(function() {
    stopApp()
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)
