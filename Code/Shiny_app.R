library(shiny)
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(rjson)
library(Seurat)
library(readr)
library(RColorBrewer)
#specification csv:
#coords must be called x and y
#cluster by  expression column must be named cluster_expression
#cluster by Pixel column must be named cluster_pixel

#image must be in PNG format

#must have obj seurat with cvs file in it, with the same specifications
load("~/obj_biohackathon16um.rda")
obj_seurat <- obj_biohackathon16um
obj_pixel <- subset(obj_seurat, 
                    subset = cluster_pixel != is.na(cluster_pixel))

options(shiny.maxRequestSize=3000*1024^2) 
# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("Visium Spatial Coordinates Plotter with User-Controlled Transformations"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("csvFile", "Upload Coordinate CSV (Background and Clusters)", accept = c(".csv")),
      fileInput("pngFile", "Upload Reference PNG Image", accept = c(".png")),
      uiOutput("clusterSelectUI"),
      checkboxInput("PlotAllClusterExp", "Plot All Expression Cluster", value = FALSE),
      numericInput("xScale", "X Scale Factor", value = 1, step = 0.01),
      numericInput("yScale", "Y Scale Factor", value = 1, step = 0.01),
      numericInput("xOffset", "X Offset", value = 0, step = 10),
      numericInput("yOffset", "Y Offset", value = 0, step = 10),
      checkboxInput("rotateCoords", "Rotate Coordinates 90 Degrees Counter Clockwise", value = FALSE),
      checkboxInput("flipHorizontal", "Flip Coordinates Horizontally (Mirror X)", value = FALSE),
      checkboxInput("flipVertical", "Flip Coordinates Vertically (Mirror Y)", value = FALSE),
      uiOutput("clusterSelectpixel"),
      checkboxInput("PlotAllClusterPixel", "Plot All Pixel Cluster", value = FALSE),
      #uiOutput("SeuratClusterSelectExp"),
      #checkboxInput("SeuratClusterAllExp", "Plot All Expression Cluster", value = FALSE),
      uiOutput("SeuratClusterSelectPixel"), 
      checkboxInput("SeuratClusterAllPixel", "Plot All Pixel Cluster", value = FALSE),
      uiOutput("SeuratClusterSelectExpPcells"), 
      checkboxInput("SeuratClusterAllExpPcells", "Plot All Expression Cluster (with Pixel Cells)", value = FALSE),
      textInput("geneSelected", "Gene Expression"),
      
      actionButton("resetBtn", "Reset Adjustments")
    ),
    
    mainPanel(
      fluidRow(
        column(6, plotOutput("scatterPlot")),
        column(6, plotOutput("scatterPlotAllClusterExp")),
        column(10, plotOutput("pngOverlay")),
        column(6, plotOutput("scatterPlot_pixel")),
        column(6, plotOutput("scatterPlotAllClusterPixel"))
      ),
      
      verbatimTextOutput("SeuratWarning"),
      
      fluidRow(
        #column(6, plotOutput("Seurat_cluster_exp")),
        #column(6, plotOutput("Seurat_cluster_exp_all")),
        column(6, plotOutput("Seurat_cluster_pixel")),
        column(6, plotOutput("Seurat_cluster_pixel_all")),
        column(6, plotOutput("Seurat_cluster_exp_Pcells")),
        column(6, plotOutput("Seurat_cluster_exp_all_Pcells")),
        column(10, plotOutput("expression_plot"))
      ),
      
      #verbatimTextOutput("debugInfo")
    )
    
  )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
  
  # Reactive value to hold uploaded coordinates
  allCoords <- reactiveVal(NULL)
  
  # Load the coordinate CSV file
  observeEvent(input$csvFile, {
    coord_data <- read.csv(input$csvFile$datapath)
    allCoords(coord_data)
  })
  
  # Reactive value to hold uploaded PNG image
  pngImage <- reactiveVal(NULL)
  pngDims <- reactiveVal(NULL)
  
  # Load the PNG file and get its dimensions
  observeEvent(input$pngFile, {
    img <- readPNG(input$pngFile$datapath)
    pngImage(img)
    pngDims(dim(img))
  })
  
  # Dynamically generate the cluster selection(expression)  UI based on the CSV file
  output$clusterSelectUI <- renderUI({
    req(allCoords())
    cluster_ids <- na.omit(unique(allCoords()$cluster_expression))
    selectInput("cluster_expression", "Select Expression Cluster:", choices = cluster_ids)
  })
  
  #cluster pixel selection
  output$clusterSelectpixel <- renderUI({
    req(allCoords())
    cluster_ids <- unique(allCoords()$cluster_pixel)
    selectInput("cluster_pixel", "Select Pixel Cluster:", choices = cluster_ids)
  })
  
  #cluster by expression, seurat object
  output$SeuratClusterSelectExp <- renderUI({
    req(obj_seurat)
    cluster_ids <- na.omit(unique(obj_seurat@meta.data$cluster_expression))
    selectInput("cluster_expression_seurat", "Select Seurat Expression Cluster:", choices = cluster_ids)
  })
  
  #cluster by pixel, seurat object
  output$SeuratClusterSelectPixel <- renderUI({
    req(obj_pixel)
    cluster_ids <- unique(obj_pixel@meta.data$cluster_pixel)
    selectInput("cluster_pixel_seurat", "Select Seurat Pixel Cluster:", choices = cluster_ids)
  })
  
  #cluster by expression, seurat object with pixel cells
  output$SeuratClusterSelectExpPcells <- renderUI({
    req(obj_pixel)
    cluster_ids <- unique(obj_pixel@meta.data$cluster_expression)
    selectInput("cluster_expression_seurat_Pcells", "Select Seurat Expression Cluster (Pixel Cells):", choices = cluster_ids)
  })
  
  
  # Reset adjustments
  observeEvent(input$resetBtn, {
    updateNumericInput(session, "xScale", value = 1)
    updateNumericInput(session, "yScale", value = 1)
    updateNumericInput(session, "xOffset", value = 0)
    updateNumericInput(session, "yOffset", value = 0)
    updateCheckboxInput(session, "rotateCoords", value = FALSE)
    updateCheckboxInput(session, "flipHorizontal", value = FALSE)
    updateCheckboxInput(session, "flipVertical", value = FALSE)
    updateCheckboxInput(session, "PlotAllClusterExp", value = FALSE)
    updateCheckboxInput(session, "PlotAllClusterPixel", value = FALSE)
    updateCheckboxInput(session, "geneSelected", value = FALSE)
    updateCheckboxInput(session, "SeuratClusterAllExp", value = FALSE) 
    updateCheckboxInput(session, "SeuratClusterAllPixel", value = FALSE) 
    updateCheckboxInput(session, "SeuratClusterAllExpPcells", value = FALSE) 
    
    
  })
  
  # Render the plot with background and overlayed cluster
  output$scatterPlot <- renderPlot({
    req(allCoords())
    coord <- allCoords()
    
    # Apply rotation if selected
    if (input$rotateCoords) {
      coord <- coord %>%
        mutate(
          x_rot = y,
          y_rot = -x
        )
    } else {
      coord <- coord %>%
        mutate(
          x_rot = x,
          y_rot = y
        )
    }
    
    # Apply horizontal flip if selected
    if (input$flipHorizontal) {
      coord <- coord %>%
        mutate(
          x_rot = -x_rot
        )
    }
    
    # Apply vertical flip if selected
    if (input$flipVertical) {
      coord <- coord %>%
        mutate(
          y_rot = -y_rot
        )
    }
    
    # Plot the spatial coordinates as the background
    p <- ggplot(coord, aes(x = x_rot, y = y_rot)) +
      geom_point(color = "grey", size = 0.01) +
      scale_y_reverse() +
      labs(title = paste0("Spatial Coordinates - Expression Cluster: ", input$cluster_expression),
           x = "X Coordinate",
           y = "Y Coordinate") +
      theme_classic() +
      coord_fixed(ratio = 1)
    
    
    # Overlay the selected cluster if available
    if (!is.null(input$cluster_expression)) {
      cl_coords_filtered <- coord %>% filter(cluster_expression == input$cluster_expression)
      p <- p + geom_point(data = cl_coords_filtered, aes(x = x_rot, y = y_rot), color = "red", size = 0.5)
    }
    
    print(p)
  })
  
  
  #render pixel cluster
  output$scatterPlot_pixel <- renderPlot({
    req(allCoords())
    coord <- allCoords() 
    
    # Apply rotation if selected
    if (input$rotateCoords) {
      coord <- coord %>%
        mutate(
          x_rot = y,
          y_rot = -x
        )
    } else {
      coord <- coord %>%
        mutate(
          x_rot = x,
          y_rot = y
        )
    }
    
    # Apply horizontal flip if selected
    if (input$flipHorizontal) {
      coord <- coord %>%
        mutate(
          x_rot = -x_rot
        )
    }
    
    # Apply vertical flip if selected
    if (input$flipVertical) {
      coord <- coord %>%
        mutate(
          y_rot = -y_rot
        )
    }
    
    # Plot the spatial coordinates as the background
    p <- ggplot(coord, aes(x = x_rot, y = y_rot)) +
      geom_point(color = "grey", size = 0.01) +
      scale_y_reverse() +
      labs(title = paste0("Spatial Coordinates - Pixel Cluster: ", input$cluster_pixel),
           x = "X Coordinate",
           y = "Y Coordinate") +
      theme_classic() +
      coord_fixed(ratio = 1)
    
    
    # Overlay the selected cluster_pixel if available
    if (!is.null(input$cluster_pixel)) {
      cl_coords_filtered <- coord %>% filter(cluster_pixel == input$cluster_pixel)
      p <- p + geom_point(data = cl_coords_filtered, aes(x = x_rot, y = y_rot), color = "red", size = 0.5)
    }
    
    print(p)
  })
  
  #render gene expression plot
  output$expression_plot <- renderPlot({
    req(obj_seurat)
    req(input$geneSelected)
    
    # Apply vertical flip if selected
    if (!is.null(input$geneSelected)) {
      gene <- input$geneSelected
      exp_plot <- SpatialFeaturePlot(obj_seurat, features = gene,
                                     pt.size.factor = 5)
      print(exp_plot)
    }}
  )
  
  
  # Render the PNG image with overlayed coordinates
  output$pngOverlay <- renderPlot({
    req(pngImage(), allCoords(), pngDims(), input$cluster_expression)
    
    # Get the PNG image and its dimensions
    img <- pngImage()
    img_dims <- pngDims()
    
    # Get the coordinates and filter by selected cluster
    coord <- allCoords()
    
    # Apply rotation if selected
    if (input$rotateCoords) {
      coord <- coord %>%
        mutate(
          x_rot = y,
          y_rot = -x
        )
    } else {
      coord <- coord %>%
        mutate(
          x_rot = x,
          y_rot = y
        )
    }
    
    # Apply horizontal flip if selected
    if (input$flipHorizontal) {
      coord <- coord %>%
        mutate(
          x_rot = -x_rot
        )
    }
    
    # Apply vertical flip if selected
    if (input$flipVertical) {
      coord <- coord %>%
        mutate(
          y_rot = -y_rot
        )
    }
    
    cl_coords_filtered <- coord %>% filter(cluster_expression == input$cluster_expression)
    
    # Apply scaling and offset
    cl_coords_filtered <- cl_coords_filtered %>%
      mutate(
        x_scaled = (x_rot + input$xOffset - min(x_rot) + 1)* input$xScale,
        y_scaled = (y_rot + input$yOffset - min(y_rot) + 1)* input$yScale
      )
    
    # Plot the PNG image with overlayed coordinates
    grid.raster(img)
    grid.points(cl_coords_filtered$x_scaled, #img_dims[2] - 
                cl_coords_filtered$y_scaled, 
                pch = 16, gp = gpar(col = "red", cex = 0.3))
  })
  
  
  
  #render all cluster expression
  output$scatterPlotAllClusterExp <- renderPlot({
    req(allCoords())
    coord <- allCoords() 

    # Apply rotation if selected
    if (input$rotateCoords) {
      coord <- coord %>%
        mutate(
          x_rot = y,
          y_rot = -x
        )
    } else {
      coord <- coord %>%
        mutate(
          x_rot = x,
          y_rot = y
        )
    }
    
    # Apply horizontal flip if selected
    if (input$flipHorizontal) {
      coord <- coord %>%
        mutate(
          x_rot = -x_rot
        )
    }
    
    # Apply vertical flip if selected
    if (input$flipVertical) {
      coord <- coord %>%
        mutate(
          y_rot = -y_rot
        )
    }
    
    # Plot the spatial coordinates 
    if (!is.null(input$cluster_expression)) {
      if(input$PlotAllClusterExp){
        coord$cluster_expression <- as.character(coord$cluster_expression)
        p <- ggplot(coord, aes(x = x_rot, y = y_rot, color = cluster_expression)) +
          geom_point(size = 1) +
          scale_y_reverse() +
          labs(title = "Spatial Coordinates - All Expression Clusters",
               x = "X Coordinate",
               y = "Y Coordinate") +
          theme_classic() +
          coord_fixed(ratio = 1)
        
      }}
    
    print(p)
  })
  
  
  #render all cluster pixel
  output$scatterPlotAllClusterPixel <- renderPlot({
    req(allCoords())
    coord <- allCoords() 
    #coord <- coord[rownames(coord) %in% rownames(na.omit(coord)),]
    
    # Apply rotation if selected
    if (input$rotateCoords) {
      coord <- coord %>%
        mutate(
          x_rot = y,
          y_rot = -x
        )
    } else {
      coord <- coord %>%
        mutate(
          x_rot = x,
          y_rot = y
        )
    }
    
    # Apply horizontal flip if selected
    if (input$flipHorizontal) {
      coord <- coord %>%
        mutate(
          x_rot = -x_rot
        )
    }
    
    # Apply vertical flip if selected
    if (input$flipVertical) {
      coord <- coord %>%
        mutate(
          y_rot = -y_rot
        )
    }
    
    # Plot the spatial coordinates as the background
   
    if (!is.null(input$cluster_pixel)) {
      if(input$PlotAllClusterPixel){
        coord$cluster_pixel <- as.character(coord$cluster_pixel)
        p <- ggplot(coord, aes(x = x_rot, y = y_rot, color = cluster_pixel)) +
          geom_point(color = "grey", size = 0.01) +
          geom_point(data = (coord %>% filter(cluster_pixel != is.na(cluster_pixel)))) +
          scale_y_reverse() +
          labs(title = "Spatial Coordinates - All Pixel Clusters",
               x = "X Coordinate",
               y = "Y Coordinate") +
          theme_classic() +
          coord_fixed(ratio = 1)
        
      }}
    
    print(p)
  })
  
  # 
  # # Debugging output to explore coordinate ranges
  # output$debugInfo <- renderPrint({
  #   req(allCoords())
  #   coord <- allCoords()
  #   
  #   # Print summaries for debugging
  #   cat("Coordinate Summary:\n")
  #   print(summary(coord))
  #   cat("\nPNG Dimensions:\n")
  #   print(pngDims())
  # })
  # 
  # output$SeuratWarning <- renderPrint({
  #   print("Seurat Object Part")
  # })
  
  
  
  # #rendering expression seurat cluster overlap
  # output$Seurat_cluster_exp <- renderPlot({
  #   req(obj_seurat)
  # 
  #   obj_aux <- obj_seurat
  #   obj_aux@meta.data$cluster_expression <- as.character(obj_aux@meta.data$cluster_expression)
  #   
  #   # Plot the spatial tissue image in the background
  #   if (!is.null(input$cluster_expression_seurat)) {
  #     obj_aux <- subset(obj_aux, 
  #                       subset = cluster_expression == input$cluster_expression_seurat)
  #     obj_aux@meta.data <- obj_aux@meta.data[names(obj_aux@active.ident),]
  #     p <-  SpatialDimPlot(obj_aux, pt.size.factor = 5, group.by = "cluster_expression")
  #   }
  #   
  #   print(p)
  # })
  # 
  # 
  # #render all cluster seurat expression
  # output$Seurat_cluster_exp_all <- renderPlot({
  #   req(obj_seurat)
  #   obj_aux <- obj_seurat
  #   obj_aux@meta.data$cluster_expression <- as.character(obj_aux@meta.data$cluster_expression)
  #   
  #   # Plot the spatial coordinates 
  #   if(input$SeuratClusterAllExp){
  #     
  #     p <- SpatialDimPlot(obj_aux, 
  #                         group.by = "cluster_expression", 
  #                         pt.size.factor = 5)
  #   }
  #   
  #   print(p)
  # })
  # 
  
  #rendering pixel seurat cluster overlap
  output$Seurat_cluster_pixel <- renderPlot({
    req(obj_pixel)
    
    obj_aux <- obj_pixel
    obj_aux@meta.data$cluster_pixel <- as.character(obj_aux@meta.data$cluster_pixel)
    
    # Plot the spatial tissue image in the background
    if (!is.null(input$cluster_pixel_seurat)) {
      obj_aux <- subset(obj_aux, 
                        subset = cluster_pixel == input$cluster_pixel_seurat)
      obj_aux@meta.data <- obj_aux@meta.data[names(obj_aux@active.ident),]
      p <-  SpatialDimPlot(obj_aux, pt.size.factor = 15, group.by = "cluster_pixel") + 
        labs(title = paste0("Seurat Pixel Cluster: ", input$cluster_pixel_seurat),
             x = "X Coordinate",
             y = "Y Coordinate",
        ) +
        theme_classic() & NoLegend()
    }
    print(p)
  })
  
  
  #render all cluster seurat Pixel
  output$Seurat_cluster_pixel_all <- renderPlot({
    req(obj_pixel)
    obj_aux <- obj_pixel
    obj_aux@meta.data$cluster_pixel <- as.character(obj_aux@meta.data$cluster_pixel)
    
    # Plot the spatial coordinates 
    if(input$SeuratClusterAllPixel){
      
      p <- SpatialDimPlot(obj_aux, 
                          group.by = "cluster_pixel", 
                          pt.size.factor = 15) + 
        labs(title = "Seurat All Pixel Cluster",
             x = "X Coordinate",
             y = "Y Coordinate",
        ) +
        theme_classic() & NoLegend()
    }
    
    print(p)
  })
  
  
  #rendering expression seurat cluster overlap with pixel cells
  output$Seurat_cluster_exp_Pcells <- renderPlot({
    req(obj_pixel)
    
    obj_aux <- obj_pixel
    obj_aux@meta.data$cluster_expression <- as.character(obj_aux@meta.data$cluster_expression)
    
    # Plot the spatial tissue image in the background
    if (!is.null(input$cluster_expression_seurat_Pcells)) {
      obj_aux <- subset(obj_aux, 
                        subset = cluster_expression == input$cluster_expression_seurat_Pcells)
      obj_aux@meta.data <- obj_aux@meta.data[names(obj_aux@active.ident),]
      p <-  SpatialDimPlot(obj_aux, pt.size.factor = 15, group.by = "cluster_expression") + 
        labs(title = paste0("Seurat Expression Cluster: ", input$cluster_expression_seurat_Pcells),
             subtitle = "Cells With Pixel information",
             x = "X Coordinate",
             y = "Y Coordinate",
        ) +
        theme_classic() & NoLegend()
    }
    
    print(p)
  })
  
  
  #render all cluster seurat expression with Pixel cells
  output$Seurat_cluster_exp_all_Pcells <- renderPlot({
    req(obj_pixel)
    obj_aux <- obj_pixel
    obj_aux@meta.data$cluster_expression <- as.character(obj_aux@meta.data$cluster_expression)
    
    obj_aux@meta.data$cluster_expression <- factor(obj_aux@meta.data$cluster_expression, 
                                                   levels = c(unique(obj_aux@meta.data$cluster_expression)))

    #Idents(obj_aux) <- "cluster_expression"
    # Plot the spatial coordinates 
    if(input$SeuratClusterAllExpPcells){
      
      p <- SpatialDimPlot(obj_aux, 
                          group.by = "cluster_expression", 
                          pt.size.factor = 15) + 
        labs(title = "Seurat All Expression Cluster",
             subtitle = "Cells With Pixel information",
             x = "X Coordinate",
             y = "Y Coordinate",
        ) +
        theme_classic() & NoLegend()
    }
    
    print(p)
  })
  
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

