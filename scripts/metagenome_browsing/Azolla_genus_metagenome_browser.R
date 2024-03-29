# load libraries
library(shiny)        # basic shiny 
library(ggplot2)      # my favourite graphing library
library(data.table)   # for efficient handling of very big tables
library(ggtext)       # for markdown parsing inside ggplot

# import data
header <-  c('assembly', 'sample', 'scaffolded', 'node',    'length',  'coverage', 'ORFs',    'ORFs_classified', 'superkingdom', 'phylum', 'class',  'order',  'family', 'genus', 'species')
classes <- c('factor',   'factor', 'factor',     'numeric', 'numeric', 'numeric',  'numeric', 'numeric',         'factor',       'factor', 'factor', 'factor', 'factor', 'factor', 'factor')
metrics_shiny <- fread(input = "assembly_stats_and_taxonomy_2500.tab",
                 sep = '\t',
                 fill = TRUE,
                 header = F,
                 stringsAsFactors = T,
                 data.table = T,
                 col.names = header,
                 colClasses = classes
)
levels(metrics_shiny$assembly) <- c('doublefiltered',
                                    'doublefiltered',
                                    'singlefiltered'
)

# polish metadata

## make a separate factor to distinct hybrid assemblies from single library assemblies:
singlelibs <- c("Azfil_lab_250", 
                "Azfil_lab_500", 
                "Azfil_lab_800", 
                "Azfil_minuscyano_170",
                "Azfil_minuscyano_350", 
                "Azfil_wild_galgw_E_1", 
                "Azfil_wild_galgw_E_2", 
                "Azfil_wild_galgw_E_3", 
                "Azfil_wild_galgw_P_2", 
                "Azfil_wild_galgw_P_3", 
                "Azfil_wild_galgw_P_4"
)

hybridlibs <- c('Azfil_wild', 
                'Azfil_minuscyano', 
                'Azfil_lab'
)

metrics_shiny$assemblytype <- factor(x = 'single library',
                                     levels = c('single library',
                                                'hybrid library',
                                                'partial library'
                                                )
)

metrics_shiny[sample %in% singlelibs]$assemblytype <- 'partial library'
metrics_shiny[sample %in% hybridlibs]$assemblytype <- 'hybrid library'
rm(singlelibs,hybridlibs)

## Make a separate factor to indicate species name rather than only a library code
metrics_shiny$hostspecies <- factor(x = '*A. filiculoides*',
                                levels = c('*A. caroliniana*',
                                           '*A. microphylla*',
                                           '*A. mexicana*',
                                           '*A. filiculoides*',
                                           '*A. rubra*',
                                           '*A. pinnata*',
                                           '*A. nilotica*'
                                           )
)

metrics_shiny[sample %in% c('Azmex_IRRI_486'   )]$hostspecies <- '*A. mexicana*'
metrics_shiny[sample %in% c('Azmic_IRRI_456'   )]$hostspecies <- '*A. microphylla*'
metrics_shiny[sample %in% c('Aznil_IRRI_479'   )]$hostspecies <- '*A. nilotica*'
metrics_shiny[sample %in% c('Azrub_IRRI_479'   )]$hostspecies <- '*A. rubra*'
metrics_shiny[sample %in% c('Azspnov_IRRI1_472')]$hostspecies <- '*A. caroliniana*'
metrics_shiny[sample %in% c('Azspnov_IRRI2_489')]$hostspecies <- '*A. caroliniana*'

# Define UI

ui <- fluidPage(
## Application title
    titlePanel("Azolla genus metagenome taxonomy browser"),
    markdown("This R Shiny application is associated with the [Azolla genus metagenome project](https://github.com/lauralwd/Azolla_genus_metagenome) by [Laura Dijkhuizen](https://www.uu.nl/medewerkers/lwdijkhuizen). 
             It aims to reveal similarities amongst metagenome assemblies in a single glance.
             Assemblies were derived from public sequencing data of *Azolla* plants, known for their endophytic microbes and systematic vertical transmission of those microbes across generations of plants.
             The assembled contigs and scaffolds are augmented with taxonomy estimations by [CAT](https://github.com/dutilh/CAT), then plotted in a dot-plot by contig/scaffold length (x) and depth in the assembly graph (y).
             The app user may manipulate the input assemblies of the plot and filter by taxonomy metadata.
             The [code for this app](https://github.com/lauralwd/Azolla_genus_metagenome/blob/master/scripts/metagenome_browsing/Azolla_genus_metagenome_browser.R) falls under the [license](https://github.com/lauralwd/Azolla_genus_metagenome/blob/master/LICENSE) of the Github repository in which it is made public.
             Please read the figure legend below for further information."),

## Sidebar 
    sidebarLayout(
        sidebarPanel(
          h3("First, choose input data to plot in the graph:"),
            selectInput("filterstage","Filter stage of the assembly data:",
                        choices = list("host filtered" = "doublefiltered",
                                       "double filtered" = "singlefiltered",
                                       "both" = " "
                        ),
                        selected = "singlefiltered"
            ),
            selectInput("assemblylevel","Assembly strategy:",
                        choices = list("Single library assemblies" = "hybrid library",
                                       "Hybrid assemblies" = "partial library",
                                       "All assemblies" = " "
                        ),
                        selected = "partial library"
            ),
            selectInput("format","Contigs or scaffolds:",
                        choices = as.list(levels(metrics_shiny$scaffolded)),
                        selected = "scaffolds"
            ),
            h3("Second, choose how to colour the the dotplot:"),
            selectInput("taxonomy","Taxonomy level:",
                        choices = list("superkingdom" = "superkingdom",
                                       "phylum" = "phylum",
                                       "class" = "class",
                                       "order" = "order",
                                       "family" = "family",
                                       "genus" = "genus",
                                       "species" = "species",
                                       "single colour" = "assembly"
                                       ),
                        selected = "order"
                        ),
            selectInput("dotsize","Dot size represents:",
                        choices = list("Open reading frames" = "ORFs",
                                       "Open reading frames classified" = "ORFs_classified",
                                       "Contig/scaffold length" = "length"
                        )
                      ),
            h3("Third, filter the data displayed:"),
            checkboxGroupInput(inputId = 'filter',label = "Remove taxa on superkingdom level",
                               choices = as.list(levels(metrics_shiny$superkingdom)),
                               selected = list('NA','','not classified')
            ),
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
                        ),
            sliderInput(inputId = 'minTaxSize',
                        label = 'Minimum size of a taxonomic group to be displayed in the plot (in Mbase)',
                        min = 0,
                        max = 100,
                        value = 2),
            uiOutput("filter_fine"),
          downloadButton('downloadPlot','Download plot svg')
          #width = 3
        ),

## The main panel
    mainPanel(
      #width = 7,
      div(style = "position:relative",
          plotOutput("plot",hover = "plot_hover",
                     click = "plot_click",
                     brush = "plot_brush",
                     width = '100%',
                     height = '1000px'
                     ),
          uiOutput("hover_info"),
          markdown("**Figure legend:** Metagenome assemblies of 6 species of the fern genus *Azolla* (horizontal panels).
                   Sequencing data were derived from three public projects:
                   first, the ['Azolla genome project'](https://doi.org/10.1038/s41477-018-0188-8) data: [PRJNA430527](https://www.ebi.ac.uk/ena/browser/view/PRJNA430527)
                   second, the ['foul play in the pocket paper'](https://doi.org/10.1111/nph.14843) data: [PRJEB19522](https://www.ebi.ac.uk/ena/browser/view/PRJEB19522)
                   third, the original sequencing data for [the *Nostoc azollae*](https://doi.org/10.1371/journal.pone.0011486) genome paper: [PRJNA30807](https://www.ebi.ac.uk/ena/browser/view/PRJNA30807)
                   .
                   Sequencing reads of the former two projects was rid of host plant DNA reads by mapping to the *Azolla filiculoides* genome version 1.1 available at [fernbase.org](https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/).
                   The remaining reads were assembled with SPAdes in metagenome mode, resulting contigs and scaffolds were then assigned an approximate taxonomy with the Contig Annotation Tool [(CAT)](https://github.com/dutilh/CAT).
                   Contigs assigned 'Eukaryote' were then used for a second filtering step of the previously filtered data, then assembled again with SPAdes and assigned taxonomy again with CAT.
                   Hence the datasets are termed 'hostfiltered' and 'doublefiltered' (vertical panels). 
                   Only of the doublefiltered data, hybrid assemblies were generated per plant accession (hybrid assemblies; assemblies using more than one sequencing library as input).
                   Using CAT output, a table can be generated with contig/scaffold length, depth, details on open reading frames (ORFs) and details on taxonomy.
                   To this end, [this particular script](https://github.com/lauralwd/Azolla_genus_metagenome/blob/master/scripts/make_assembly_stats_and_taxonomy.bash) was used available at the project [Github repository](https://github.com/lauralwd/Azolla_genus_metagenome).
                   The graph displays the metagenome assemblies as a dot-plot with contig/scaffold length on the x-axis and depth in the assembly graph on the y-axis; both are log10 transformed.
                   By default, dot colour represents the order, and the size represents ORF count but a user can manipulate both.
                   For clarity, contigs/scaffolds without ORFs found (classified 'blank'), and contigs scaffolds with ORFs found but without taxonomy (classified 'not classified') are filtered out before displaying the plot.
                   Second, noisy contigs/scaffolds are omitted by default for clarity. 
                   The default threshold for noisy is when contigs/scaffolds have fewer than 5 ORFs classified or when a taxonomic group amounts to less than 2Mbase in the entire figure.
                   "),
          h2('Selection details'),
          uiOutput("brush_info"),
          h2('Taxa table'),
          markdown("Taxa size in Mbase are split up by sample and taxon level in this table. Only groups that amount to more than 1Mbase are shown."),
          dataTableOutput(outputId = 'tableout')
      )
    )
)
)

# Define the server backbone
server <- function(input, output) {
  `%notin%` <- Negate(`%in%`)
  
## First, define a function containing the data set, filtered according to input settings.
  metrics_subset <- shiny::reactive({
    dt <- metrics_shiny[scaffolded == input$format &
                        superkingdom %notin% input$filter & 
                        ORFs_classified >= input$minORFs_classified &
                        ORFs >= input$minORFs &
                        assembly != input$filterstage &
                        assemblytype != input$assemblylevel] 
    dt[, taxonomy := get(input$taxonomy)]
    
    low_abundant <- dt[,
                       .(length_mb=sum(length)/1000000),
                       .(get(input$taxonomy))][length_mb <= input$minTaxSize][,1] 
    low_abundant <- as.matrix(low_abundant)
    dt[taxonomy %in% low_abundant, 'taxonomy'] <- 'low abundant'
    dt
    }) %>% debounce(2000)

  output$filter_fine <- renderUI({
    tax_list <- unique(metrics_subset()[,taxonomy])
    choices <- as.character(tax_list[order(tax_list)])
    names(choices) <-  tax_list[order(tax_list)]
    checkboxGroupInput(inputId = "fine_filter",
                       label = "Select taxa to omit from the plot",
                       choices = choices,
                       selected = list('low abundant','not classified')
                       )
  }) 
  
## Second, using the pre-filtered data, a plot is build with ggplot2
  plotcontainer <- reactiveValues() 
  output$plot <- renderPlot({
    `%notin%` <- Negate(`%in%`)
      textsize <- 12
      xlabel <- 'contig length'
      ylabel <- 'contig coverage'
      if(input$format == 'scaffolds'){xlabel <- 'scaffold length' 
                                       ylabel <- 'scaffold coverage'} 
      
      dotplot <- ggplot(metrics_subset()[taxonomy %notin% input$fine_filter],
                            aes_string(x='length',
                                       y='coverage',
                                       size=input$dotsize,
                                       alpha=0.001,
                                       col='taxonomy'
                                       )
                            )
      dotplot <- dotplot + facet_grid(assembly ~ hostspecies + sample)
      dotplot <- dotplot + geom_point()
      dotplot <- dotplot + scale_x_log10(limits = c(10000,1000000))
      dotplot <- dotplot + scale_y_log10(limits = c(0.05,10000))
#      dotplot <- dotplot + scale_color_brewer(type = 'qual',palette = 'Set1' ,input$taxonomy)
      dotplot <- dotplot + scale_alpha(guide = F)
      dotplot <- dotplot + labs(x=xlabel,
                                y=ylabel
                                        )
      dotplot <- dotplot + theme_classic()
      dotplot <- dotplot + theme(legend.position = "bottom",
                                 legend.justification = 'center',
                                 legend.direction = 'horizontal',
                                 legend.box = 'vertical',
                                 text = element_text(size = textsize),
                                 strip.text  = element_markdown(size = 10),
                                 legend.text = element_text(size = textsize),
                                 axis.text.x = element_text(angle = 40,hjust = 1,size = textsize),
                                 axis.text.y = element_text(textsize),
                                 axis.title  = element_text(textsize)
                                 )
      #dotplot <- dotplot + annotation_logticks()
      plotcontainer$plot <- dotplot
      dotplot
    })
  
## third, an interactive window is defined which appears when the plot is clicked.
    output$hover_info <- renderUI({
      hover <- input$plot_hover
      # if pointer is outside plot (no x cord) then return markdown element with hover instructions
      if (length(hover$x) == 0) return(markdown("Hover over the plot to see more information about any particular contig or scaffold:"))
      # subset data table according to hover pane, taxonomy fine filter, and hover coordinates
      point <- nearPoints(df = metrics_subset()[assembly == input$plot_hover[[8]] &
                                                sample   == input$plot_hover[[6]] &
                                                taxonomy %notin% input$fine_filter],
                          xvar = 'length',
                          yvar = 'coverage',
                          coordinfo = input$plot_hover,
                          addDist = TRUE)
      
      # bellow section should allow to show tooltip on the plot near hovering location, but this currently does not work well.
      #left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      #top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      #left_px <- log10(hover$range$left + left_pct * (hover$range$right - hover$range$left))
      #top_px <- log10(hover$range$top + top_pct * (hover$range$bottom - hover$range$top))
      
      #style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
      #                "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      wellPanel(#style = style,
                p(HTML(
                    paste0(as.character(
                      unique(point[,taxonomy])
                      ))
                  ))
      )
    })
    
## fourth, calculate a table based on brush info
    output$brush_info <- renderUI({  
      brush <- input$plot_brush
      # if no square/brush is drawn, return this instruction:
      if (is.null(brush)) return(markdown("Select a rectangle in the plot to see more information about a particular set of contigs or scaffolds"))
      # subset data table according to brush dimensions and panel and taxonomy fine filter
      square <- metrics_subset()[assembly == brush[[7]] &
                                 sample   == brush[[5]] &
                                 length   >= brush$xmin &
                                 length   <= brush$xmax &
                                 coverage >= brush$ymin &
                                 coverage <= brush$ymax &
                                 taxonomy %notin% input$fine_filter
                                 ]
      # render final table
      renderDataTable({
      square[,
             .(contig_count    = length(length), 
               length_mb       = round( sum(length)/1000000),
               ORF_count       = sum(ORFs),
               ORFs_classified = sum(ORFs_classified),
               mean_coverage   = round(mean(coverage),2),
               sd_coverage    = round(sd(coverage),2)
             ),
             by='taxonomy']
    },
    options = list(scrollX = TRUE,
                   paging = F,
                   searching = F)
    )
    })
    
## fifth, a table is rendered displaying the top 14 taxa at the given filter, their contig count and total size in Mbase
    output$tableout<- renderDataTable({
      df_long <- (metrics_subset()[taxonomy %notin% input$fine_filter,
                                 .(length_mb       = round( sum(length)/1000000)
                                   #ORF_count       = sum(ORFs),
                                   #ORFs_classified = sum(ORFs_classified)
                                   ),
                                 by=c('taxonomy', #by=c(eval(input$taxonomy),
                                      'assembly',
                                      'sample')]
                                [length_mb >= 1]
                                [order(-rank(length_mb))]
      )
      df_cast <- dcast.data.table(data = df_long
                                  ,formula = assembly + taxonomy ~ sample
                                  ,fun.aggregate = sum
                                  ,value.var = 'length_mb'
                                  )
      df_cast
    },
    options = list(scrollX = TRUE,
                   paging = F
                   )
    )
## sixth, render a png/svg for downloading (file format is hardcoded, sorry)
    output$downloadPlot <- downloadHandler(
      #filename = function(){ paste('Azolla_genus_metagenome','-',input$taxonomy, '.svg', sep='') },
      filename = function(){ paste('Azolla_genus_metagenome','-',input$taxonomy, '.png', sep='') },
      content  = function(file){
        ggsave(filename = file
             ,plot = plotcontainer$plot
             #,device = 'svg' 
             ,width = unit(x = 18,units = 'cm')
             ,height = unit(x = 14,units = 'cm')
             ,scale = .55 # somehow I don't get the right dimensions out and this is an approximate fix...
             # for landscape format and png:
             ,device = 'png'
             #,width = unit(x = 27,units = 'cm')
             #,height = unit(x = 18,units = 'cm')
             ,dpi = 400)
      }
    )
}

# Run the application 

shinyApp(ui = ui, server = server)
