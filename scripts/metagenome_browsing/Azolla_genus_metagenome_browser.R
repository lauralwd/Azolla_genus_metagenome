# load libraries
library(shiny)
library(ggplot2)
library(data.table)

# import data only once
header <-  c('assembly', 'sample', 'scaffolded', 'node',    'length',  'coverage', 'ORFs',    'ORFs_classified', 'superkingdom', 'phylum', 'class',  'order',  'family', 'genus', 'species')
classes <- c('factor',   'factor', 'factor',     'numeric', 'numeric', 'numeric',  'numeric', 'numeric',         'factor',       'factor', 'factor', 'factor', 'factor', 'factor', 'factor')
metrics_shiny <- fread(input = "assembly_stats_and_taxonomy_2500.tab", # -hybrid_stats_and_taxonomy.tab", #  _stats_and_taxonomy_2500.tab
                 sep = '\t',
                 fill = TRUE,
                 header = F,
                 stringsAsFactors = T,
                 data.table = T,
                 col.names = header,
                 colClasses = classes
)
levels(metrics_shiny$assembly) <- c('doublefiltered','doublefiltered','singlefiltered')
levels(metrics_shiny$superkingdom)



singlelibs <- c("Azfil_lab_250", "Azfil_lab_500", "Azfil_lab_800", "Azfil_minuscyano_170", "Azfil_minuscyano_350", "Azfil_wild_galgw_E_1", "Azfil_wild_galgw_E_2", "Azfil_wild_galgw_E_3", "Azfil_wild_galgw_P_2", "Azfil_wild_galgw_P_3", "Azfil_wild_galgw_P_4") 
hybridlibs <- c('Azfil_wild', 'Azfil_minuscyano', 'Azfil_lab')
metrics_shiny$assemblytype <- factor(x = 'single library',levels = c('single library','hybrid library','partial library'))
metrics_shiny[sample %in% singlelibs]$assemblytype <- 'partial library'
metrics_shiny[sample %in% hybridlibs]$assemblytype <- 'hybrid library'

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Azolla genus metagenome assembly browser"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("filterstage","Choose a filter stage of the assembly data",
                        choices = list("host filtered" = "doublefiltered",
                                       "double filtered" = "singlefiltered",
                                       "both" = " "
                        ),
                        selected = "singlefiltered"
            ),
            selectInput("assemblylevel","Choose an assembly strategy",
                        choices = list("Single library assemblies" = "hybrid library",
                                       "Hybrid assemblies" = "partial library",
                                       "All assemblies" = " "
                        ),
                        selected = "partial library"
            ),
            selectInput("taxonomy","Choose a taxonomy level to colour the dotplot",
                        choices = list("superkingdom",
                                       "phylum",
                                       "class",
                                       "order",
                                       "family",
                                       "genus",
                                       "species" 
                                       ),
                        selected = "order"
                        ),
            selectInput("format","Choose contigs or scaffolds as assembly output",
                        choices = as.list(levels(metrics_shiny$scaffolded)),
                        selected = "scaffolds"
            ),            
                        selectInput("dotsize","Choose what is represented by the size of the dots",
                        choices = list("Open reading frames" = "ORFs",
                                       "Open reading frames classified" = "ORFs_classified",
                                       "Contig/scaffold length" = "length"
                        )
            ),
            checkboxGroupInput(inputId = 'filter',label = "Remove taxa on superkingdom level",
                               choices = as.list(levels(metrics_shiny$superkingdom)),
                               selected = list('NA','','not classified')
            ),
#            sliderInput(inputId = 'transparency',
#                        label = "Dot transparency",
#                        min = 1,
#                        max = 1000,
#                        value = 1,
#                        round = T,
#                        step = 50
#                        ),
            sliderInput(inputId = 'minORFs_classified',
                        label = "Minimum ORFs classified by CAT",
                        min = 0,
                        max = 50,
                        value = 5,
                        round = T,
                        step = 1
                        ),
            sliderInput(inputId = 'minORFs',
                        label = "Minimum ORFs called",
                        min = 0,
                        max = 50,
                        value = 5,
                        round = T,
                        step = 1
                        )

        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plot",click = "plot_click",width = '100%',height = '1000px'),
           verbatimTextOutput("info")
          
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  `%notin%` <- Negate(`%in%`)
    output$plot <- renderPlot({
      length_dist <- ggplot(metrics_shiny[scaffolded == input$format &
                                          superkingdom %notin% input$filter & 
                                          ORFs_classified >= input$minORFs_classified &
                                          ORFs >= input$minORFs &
                                          assembly != input$filterstage &
                                          assemblytype != input$assemblylevel  ],
                            aes_string(x='length',
                                       y='coverage',
                                       size=input$dotsize,
                                       alpha=0.001,
                                       col=input$taxonomy
                                       )
                            )
      length_dist <- length_dist + facet_grid(assembly ~ sample)
      length_dist <- length_dist + geom_point()
      length_dist <- length_dist + scale_x_log10()
      length_dist <- length_dist + scale_y_log10()
      length_dist <- length_dist + labs(x='contig length',
                                        y='contig coverage'
                                        )
      length_dist <- length_dist + theme_classic()
      length_dist <- length_dist + theme(legend.position = "bottom",
      #                                   strip.text.x = element_text(angle = 80),
                                         axis.text.x = element_text(angle = 80,hjust = 1))
      length_dist
    })
    output$info <- renderText({paste0(input$plot_click)})
}




# Run the application 
shinyApp(ui = ui, server = server)
