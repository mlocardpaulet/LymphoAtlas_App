############################################################################
# Copyright CNRS 2018
# Contributor : Marie Locard-Paulet (20/09/2018) [marie.locard@ipbs.fr]
# This software is a computer program whose purpose is to visualize and inspect the LymphoAtlas phosphoproteomics data.
# This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can  use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software,that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
############################################################################

# Packages:
############################################################################
library(shiny)
library(ggplot2)
library(reshape2)
library(shinyBS)
library(plotly)

############################################################################
# Functions:
############################################################################

############################################################################

# Import data:
############################################################################
load("DataLymphoAtlas/12_PhosphoTableWithProteinWaring.RData")
coloursLines <- c("blue4", "deepskyblue4", "lightseagreen", "mediumspringgreen")
options(max.print=nrow(export))

# App:
############################################################################
## UI:
############################################################################
ui <- fluidPage(
  
  # Application title
  titlePanel("LymphoAtlas"),
  
  # Sidebar with a slider input for the phosphorylation site of interest: 
  sidebarLayout(
    sidebarPanel(
      # selectInput("psite",
      #             "Phosphorylation site of interest:",
      #             choices = sort(as.character(export$GeneID)),
      #             selected = NULL,
      #             multiple = FALSE,
      #             selectize = TRUE)#,
      # Select the protein of interest:
      #################################
      selectizeInput("protein",
                     "Protein of interest:",
                     choices = sort(as.character(export$Accession)),
                     selected = NA,
                     multiple = FALSE),
      # Once the protein of interest is selected, select one or all the sites of the protein:
      #################################
      checkboxInput("allSites", "Show all the sites of the selected protein", FALSE), 
      conditionalPanel(condition = "input.allSites==false",
                       selectizeInput("psite",
                                      "Phosphorylation site of interest:",
                                      choices = sort(as.character(export$GeneID)),
                                      selected = "Cd3e_Y170",
                                      multiple = FALSE)#,
                       # bsTooltip("psite", 
                       #           "Select which phosphorylation site to plot",
                       #           "right")
      ),
      # Buttons for download:
      downloadButton("Download", "Download .pdf"),
      downloadButton("Download1", "Download .png"),
      downloadButton("Download2", "Download .svg")
    ),
    
    # Show a plot
    mainPanel(
      plotlyOutput("psiteplot")
    )
  ),
  # Footer
  tabsetPanel(
    tabPanel(
      HTML('<footer><font size="0.8">copyright 2018 - CNRS - All rights reserved - LymphoAtlas V1.0</font></footer>')
    )
  )
)

############################################################################
## Server:
############################################################################

server <- function(input, output) {
  psite2 <- reactive({
    if (is.null(input$psite)) {
      return(NULL)
    } else {
      return(input$psite)
    }
  })
  output$psiteplot <- renderPlotly({
    # Table preparation:
    cols <- which(grepl("MeanLoops", names(export)))
    gtab <- export[export$GeneID %in% psite2(),cols]
    names(gtab) <- gsub("MeanLoops_", "", names(gtab), fixed = T)
    gtab <- melt(gtab)
    
    validate (
      need(length(gtab$value[!is.na(gtab$value)]) > 0, "Not enough quantification values for this phosphorylation site in the data set.")
    )
    gtab$variable <- as.character(gtab$variable)
    gtab$TimePoint <- sapply(gtab$variable, function(x) {
      strsplit(x, "_", fixed = T)[[1]][3]
    })
    gtab$TimePoint[grepl("NS.", gtab$TimePoint)] <- 0 
    gtab$TimePoint[grepl("S30.", gtab$TimePoint, fixed = T)] <- 30 
    gtab$TimePoint[grepl("S15.", gtab$TimePoint, fixed = T)] <- 15 
    gtab$TimePoint[grepl("S120.", gtab$TimePoint, fixed = T)] <- 120 
    gtab$TimePoint[grepl("S300.", gtab$TimePoint, fixed = T)] <- 300 
    gtab$TimePoint[grepl("S600", gtab$TimePoint, fixed = T)] <- 600 
    gtab$Replicate <- sapply(gtab$variable, function(x) {
      strsplit(x, "_", fixed = T)[[1]][1]
    })
    gtab$TimePoint <- as.numeric(as.character(gtab$TimePoint))
    gtab <- gtab[!is.na(gtab$value),]
    
    g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) + theme_minimal() + ggtitle(paste0(psite2(), " upon TCR activation")) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + xlab("Time after stimulation (in seconds)") # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
    p <- ggplotly(g, height = 600) %>%
      config(displayModeBar = F)
    p
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

