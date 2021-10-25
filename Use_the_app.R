library(shiny)

# Please change the path to the folder scRNA_visulization
runApp("/Users/jiaxinli/Box/Li_Jiaxin/scRNA_visualize")

library(rsconnect)
rsconnect::setAccountInfo(name='single-cell-expression-visualization',
                          token='3C8680F00D1F4701457B38C121553B2D',
                          secret='Jned1xmId2778Pnz/1VfnKA+qLGK09jBYSAeQ7De')
rsconnect::deployApp('/Users/jiaxinli/Box/Li_Jiaxin/scRNA_visualize')


