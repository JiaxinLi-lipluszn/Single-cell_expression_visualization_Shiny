#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Where the data come in
#source("helpers.R")
#counties <- readRDS("counties.rds")
#library(maps)
#library(mapproj)

# The bmcite dataset
library(Seurat)

#InstallData("bmcite")
load("./Data/Sestan.fetalHuman.Psychencode.Rdata")
r = row.names(meta2)
for(i in 1:length(r)){
    r[i] = gsub("_", "-", r[i])
}
row.names(meta2) = r

r = row.names(cpm2)
for(i in 1:length(r)){
    r[i] = gsub("\\|", "-", r[i])
}
rownames(cpm2) <- r


c = colnames(cpm2)
for(i in 1:length(c)){
    c[i] = gsub("_", "-", c[i])
}

colnames(cpm2) = c


bm <- CreateSeuratObject(counts = cpm2, project = "Yale", assay = "RNA", meta.data = meta2)


meta <- bm@meta.data

available_var <- unlist(colnames(meta))

#available_var <- c("donor", "celltype.l1", "celltype.l2")
# Define the UI
ui <- fluidPage(
    titlePanel("Violion plot of gene expression "),
    
    sidebarLayout(
        sidebarPanel(
            helpText("Expression of genes in single cell data"),
            
            textInput("Gene_name", label = "Please enter the gene you want to visualize here:", 
                      value = row.names(cpm2)[1], width = NULL, placeholder = NULL),
            
            selectInput("x_axis", 
                        label = "Choose a variable to be the x axis of your violion plot",
                        choices = available_var,
                        selected = colnames(meta2)[1]),
            selectInput("split_conditions", 
                        label = "Choose a variable to be the condition to split on",
                        choices = available_var,
                        selected = colnames(meta2)[2])
        ),
        mainPanel(plotOutput("map"))
    )
)


# Server logic ----
server <- function(input, output) {
    output$map <- renderPlot({
        VlnPlot(bm, features = input$Gene_name, group.by = input$x_axis, split.by = input$split_conditions)
        
        #percent_map(var = data, color = color, legend.title = "Map", max = 100, min = 0)
    })
}
# Run app ----
shinyApp(ui, server)
