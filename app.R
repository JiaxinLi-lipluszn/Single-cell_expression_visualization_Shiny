# Author: Jiaxin Li

library(shiny)
library(Seurat)
library(readr)

load("/Users/jiaxinli/Box/Li_Jiaxin/scRNA_visualize/Data/Sestan.fetalHuman.Psychencode.Rdata")
expr_path = "/Users/jiaxinli/desktop/David/bmcite/rna_no_scale.csv"
meta_path = "/Users/jiaxinli/desktop/David/bmcite/meta.csv"

cpm2 = read_csv(expr_path)
gene_list = cpm2[,1]
cpm2 = cpm2[,-1]
rownames(cpm2) <- unlist(gene_list)


meta2 = read_csv(meta_path)
cell_list = meta2[,1]
meta2 = meta2[,-1]
rownames(meta2) <- unlist(cell_list)

r = rownames(meta2)
for(i in 1:length(r)){
    r[i] = gsub("_", "-", r[i])
}
rownames(meta2) = r

r = rownames(cpm2)
for(i in 1:length(r)){
    r[i] = gsub("\\|", "-", r[i])
}
rownames(cpm2) <- r


c = colnames(cpm2)
for(i in 1:length(c)){
    c[i] = gsub("_", "-", c[i])
}

colnames(cpm2) = c

c = colnames(cpm2)
for(i in 1:length(c)){
    c[i] = gsub("\\.", "-", c[i])
}

colnames(cpm2) = c

meta2 = data.frame(meta2)
rownames(meta2) <- unlist(colnames(cpm2))

bm <- CreateSeuratObject(counts = cpm2, project = "Yale", assay = "RNA", meta.data = meta2)

cpm <- cpm2
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
