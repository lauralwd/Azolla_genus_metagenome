# load libraries
library(shiny)
library(ggplot2)
library(data.table)
library(ggtext)

# import data
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
levels(metrics_shiny$assembly) <- c('doublefiltered',
                                    'doublefiltered',
                                    'singlefiltered'
                                    )
levels(metrics_shiny$superkingdom)

# polish metadata

# make a separate factor to distinct hybrid assemblies from single library assemblies:
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

# Make a separate factor to indicate species name rather than only a library code
metrics_shiny$species <- factor(x = '*A. filiculoides*',
                                levels = c('*A. caroliniana*',
                                           '*A. microphylla*',
                                           '*A. mexicana*',
                                           '*A. filiculoides*',
                                           '*A. rubra*',
                                           '*A. pinnata*',
                                           '*A. nilotica*'
                                           )
                                )
metrics_shiny[sample %in% c('Azmex_IRRI_486')]$species <- '*A. mexicana*'
metrics_shiny[sample %in% c('Azmic_IRRI_456')]$species <- '*A. microphylla*'
metrics_shiny[sample %in% c('Aznil_IRRI_479')]$species <- '*A. nilotica*'
metrics_shiny[sample %in% c('Azrub_IRRI_479')]$species <- '*A. rubra*'
metrics_shiny[sample %in% c('Azspnov_IRRI1_472')]$species <- '*A. caroliniana*'
metrics_shiny[sample %in% c('Azspnov_IRRI2_489')]$species <- '*A. caroliniana*'

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Azolla genus metagenome assembly browser"),
    markdown("This R Shiny application is associated with the [Azolla genus metagenome project](https://github.com/lauralwd/Azolla_genus_metagenome) by [Laura Dijkhuizen](https://www.uu.nl/medewerkers/lwdijkhuizen). 
             It aims to reveal similarities amongst metagenome assemblies in a single glance.
             Assemblies were derived from public sequencing data of *Azolla* plants, known for their endophytic microbes and systematic vertical transmission of those microbes across generations of plants.
             The assembled contigs and scaffolds are augmented with taxonomy estimations by [CAT](https://github.com/dutilh/CAT), then plotted in a dot-plot by contig/scaffold length (x) and depth in the assembly graph (y).
             The app user may manipulate the input assemblies of the plot and filter by taxonomy metadata.
             The code for this app falls under the [license](https://github.com/lauralwd/Azolla_genus_metagenome/blob/master/LICENSE) of the Github repository in which it is made public.
             Please read the info paragraph at the bottom of this page for further information."),

    # Sidebar with a slider input for number of bins 
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
            uiOutput("filter_fine")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("plot",click = "plot_click",width = '100%',height = '1000px'),
          markdown("**Figure legend:** Metagenome assemblies of 6 species of the fern genus *Azolla* (horizontal panels).
                   Sequencing data was derived from three public projects:
                   first, the ['azolla genome project'](https://doi.org/10.1038/s41477-018-0188-8) data: [PRJNA430527](https://www.ebi.ac.uk/ena/browser/view/PRJNA430527)
                   second, the ['foul play in the pocket paper'](https://doi.org/10.1111/nph.14843) data: [PRJEB19522](https://www.ebi.ac.uk/ena/browser/view/PRJEB19522)
                   third, the original sequencing data for [the *Nostoc azollae*](https://doi.org/10.1371/journal.pone.0011486) genome paper: [PRJNA30807](https://www.ebi.ac.uk/ena/browser/view/PRJNA30807)
                   .
                   Sequencing reads of the former two projects was rid of host plant DNA reads by mapping to the *Azolla filiculoides* genome version 1.1 available at [fernbase.org](https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/).
                   The remaining reads were assembled with SPAdes in metagenome mode, resulting contigs and scaffolds were then assigned an approximate taxonomy with the Contig Annotation Tool [(CAT)](https://github.com/dutilh/CAT).
                   Contigs assigned 'Eukaryote' were then used for a second filtering step of the previously filtered data, then assembled again with SPAdes, and assigned taxonomy again with CAT.
                   Hence the datasets are termed 'hostfiltered' and 'doublefiltered' (vertical panels). 
                   Only of the doublefiltered data, hybrid assemblies were generated per plant accession (hybrid assemblies; assemblies using more than one sequencing library as input).
                   Using CAT output, a table can be generated with contig/scaffold length, depth, details on open reading frames (ORFs) and details on taxonomy.
                   To this end, [this particular script](https://github.com/lauralwd/Azolla_genus_metagenome/blob/master/scripts/make_assembly_stats_and_taxonomy.bash) was used available at the project [github repository](https://github.com/lauralwd/Azolla_genus_metagenome).
                   The graph displays the metagenome assemblies as a dotplot with contig/scaffold length on the x-axis, and depth in the assembly graph on the y-axis, the latter is log10 transformed.
                   By default dot colour represents the order, and the size represents ORF count, but both can be manipulated by the user.
                   For clarity, contigs/scaffolds without ORFs found (classified 'blank') and contigs scaffolds with ORFs found, but without taxonomy (classified 'not classified') are filtered out before displaying the plot.
                   Second, noisy contigs/scaffolds are ommitted by default for clarity. 
                   The default threshold for noisy is when contigs/scaffolds have fewer than 5 ORFs classified, or when a taxonomic group amounts to less than 2Mbase in the entire figure.
                   "),
#          h2('more details'),
#          verbatimTextOutput("info"),
          h2('Taxa table'),
          markdown("Taxa present at the second stage of filtering are displayed in this table by contig count and total size in mbases."),
          tableOutput(outputId = 'tableout')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  `%notin%` <- Negate(`%in%`)
  
## First, define a function containing the dataset, filtered according to input settings.
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
    })
  output$filter_fine <- renderUI({
    checkboxGroupInput(inputId = "fine_filter",
                       label = "Select taxa to omit from the plot",
                       choices = as.list(levels(metrics_subset()[,'taxonomy'])),
                       
                       )
  })
  
## Second, using the pre-filtered data, a plot is defined with ggplot2
  output$plot <- renderPlot({
      length_dist <- ggplot(metrics_subset(),
                            aes_string(x='length',
                                       y='coverage',
                                       size=input$dotsize,
                                       alpha=0.001,
                                       col='taxonomy'
                                       )
                            )
      length_dist <- length_dist + facet_grid(assembly ~ species + sample)
      length_dist <- length_dist + geom_point()
      length_dist <- length_dist + scale_x_log10()
      length_dist <- length_dist + scale_y_log10()
      length_dist <- length_dist + labs(x='contig length',
                                        y='contig coverage'
                                        )
      length_dist <- length_dist + theme_classic()
      length_dist <- length_dist + theme(legend.position = "bottom",
                                         text = element_text(size = 12),
                                         strip.text = element_markdown(size = 12),
                                         legend.text = element_text(size = 12),
                                         axis.text.x = element_text(angle = 80,hjust = 1,size = 12)
                                         )
      length_dist
    })
## third, an interactive window is defined which appears when the plot is clicked.
## this plot wil display details on the clicked contigs
    output$info <- renderText({
      levels(metrics_subset()$taxonomy)
      paste0(input$plot_click)})
## fourth, a table is rendered displaying the top 14 taxa at the given filter, their contig count and total size in Mbase
    output$tableout<- renderTable({
      as.matrix(metrics_subset()[,
                                 .(contig_count=length(length), 
                                   length_mb=sum(length)/1000000,
                                   ORF_count=sum(ORFs),
                                   ORFs_classified=sum(ORFs_classified),
                                   mean_coverage=mean(coverage),
                                   var_coverage=var(coverage)
                                   ),
                                 by=c(eval(input$taxonomy),
                                      'assembly')]
                [length_mb >= 1]
                [order(-rank(length_mb))]
    )})
}




# Run the application 
shinyApp(ui = ui, server = server)
