############################################################################
# Copyright CNRS 2018
# Contributor : Marie Locard-Paulet (20/09/2018) [marie.locard@ipbs.fr]
# This software is a computer program whose purpose is to visualize and inspect the LymphoAtlas phosphoproteomics data.
# This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can  use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software,that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
############################################################################

# Logs are in /var/log/shiny-server

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
export$GeneNames <- sapply(export$GeneID, function(x) {
  strsplit(x, "_", fixed = T)[[1]][1]
})

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
      # Select the protein of interest:
      #################################
      selectizeInput("protein",
                     "Protein of interest:",
                     choices = sort(unique(paste0(export$GeneNames, " / ", export$Accession))),
                     selected = "Zap70 / P43404",
                     multiple = FALSE),
      bsTooltip("protein", 
                "If the list of proteins to select does not include your protein of interest, enter it in the selection box.",
                "right"),
      # Once the protein of interest is selected, select one or all the sites of the protein:
      #################################
      checkboxInput("allSites", "Show all the sites of the selected protein", FALSE), 
      conditionalPanel(condition = "input.allSites==false",
                       selectizeInput("psite",
                                      "Phosphorylation site of interest:",
                                      choices = NULL,
                                      selected = NULL,
                                      multiple = FALSE)#,
                       # bsTooltip("psite", 
                       #           "Select which phosphorylation site to plot",
                       #           "right")
      ),
      # Buttons for download:
      downloadButton("Download", "Download .pdf"),
      bsTooltip("Download", 
                "Will plot one plot per page",
                "right"),
      downloadButton("Download1", "Download .png"),
      bsTooltip("Download1", 
                "Will plot all the plots on one page",
                "right"),
      downloadButton("Download2", "Download .svg"),
      bsTooltip("Download2", 
                "Will plot all the plots on one page",
                "right")
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

server <- function(session, input, output, clientData) {
  # Define protein:
  #################
  prot2 <- reactive({
    if (is.null(input$protein)) {
      return(NULL)
    } else {
      res <- input$protein
      res <- sapply(res, function(x) {
        strsplit(x, " / ", fixed = T)[[1]][2]
      })
      return(res)
    }
  })
  # Update the field of the psite selection:
  #################
  observe({
    if (!is.null(prot2())) {
      updateSelectizeInput(session, "psite",
                           "Phosphorylation site of interest:",
                           choices = sort(as.character(export$GeneID[export$Accession == prot2()])),
                           selected = NULL)
    }
  })
  # If the checkbox "allSites" is checked, plot all the sites of the protein, else, plot only the one selected in the field "psite":
  #################
  psite2 <- reactive({
    if (input$allSites) {
      res <- export$phosphoSites[export$Accession == prot2()]
      return(res)
    } else {
      if (is.null(input$psite)) {
        return(NULL)
      } else {
        return(input$psite)
      }
    }
  })
  # Make the table for the plot:
  ################
  plotTable <- function() {
    if (is.null(prot2())) {
      return(NULL)
    } else {
      # If allSites is checked: plot all the sites for the protein selected:
      #########################
      if (input$allSites) {
        # Table preparation:
        cols <- which(grepl("MeanLoops", names(export)))
        gtab <- export[export$Accession %in% prot2(),cols]
        names(gtab) <- gsub("MeanLoops_", "", names(gtab), fixed = T)
        gtab$phosphosite <- export$GeneID[export$Accession %in% prot2()]
        gtab <- melt(gtab)
        validate (
          need(length(gtab$value[!is.na(gtab$value)]) > 0, "Select phosphorylation sites of interest")
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
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          # Table preparation:
          cols <- which(grepl("MeanLoops", names(export)))
          gtab <- export[export$GeneID %in% psite2(),cols]
          names(gtab) <- gsub("MeanLoops_", "", names(gtab), fixed = T)
          gtab <- melt(gtab)
          validate (
            need(length(gtab$value[!is.na(gtab$value)]) > 0, "Select a phosphorylation site of interest")
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
        }
      }
      return(gtab)
    }
  }
  # Make the plot:
  ################
  plotInput <- function(){
    if (is.null(plotTable())) {
      return(NULL)
    } else {
      gtab <- plotTable()
      # If allSites is checked: plot all the sites for the protein selected:
      #########################
      if (input$allSites) {
        g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 4, alpha = 0.8) + theme_minimal() + ggtitle(paste0("Sites of ", input$protein, " upon TCR activation")) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + xlab("Time after stimulation (in seconds)") + facet_wrap(~phosphosite)# + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
        # If allSites is notchecked: plot all the sites for the protein selected:
        #########################
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) + theme_minimal() + ggtitle(paste0(psite2(), " upon TCR activation")) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + xlab("Time after stimulation (in seconds)") # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
        }
      }
      return(g)
    }
  }
  output$psiteplot <- renderPlotly({
    p <- plotInput()
    ggplotly(p, height = 600) %>%
      config(displayModeBar = F)
  })
  
  # For export:
  ############
  plotInputExport <- function(){
    if (is.null(plotTable())) {
      return(NULL)
    } else {
      gtab <- plotTable()
      # If allSites is checked: plot all the sites for the protein selected:
      #########################
      if (input$allSites) {
        lg <- list()
        ylimplot <- range(gtab$value)
        for (el in unique(gtab$phosphosite)) {
          g <- ggplot(gtab[gtab$phosphosite %in% el,], aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 4, alpha = 0.8) + theme_minimal() + ggtitle(paste0(el, " upon TCR activation")) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + xlab("Time after stimulation (in seconds)") + ylim(ylimplot) # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
          lg[[length(lg) + 1]] <- g
        }
        # If allSites is notchecked: plot all the sites for the protein selected:
        #########################
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) + theme_minimal() + ggtitle(paste0(psite2(), " upon TCR activation")) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + xlab("Time after stimulation (in seconds)") # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
          lg <- list(g)
        }
      }
      return(lg)
    }
  }
  # pdf output:
  output$Download <- downloadHandler(
    validate(
      need(!is.null(psite2()), "Select a site of interest to plot.")
    ),
    filename = function(){
      if (input$allSites) {
        na <- gsub(" / ", "-", input$protein, fixed = T)
      } else {
        na <- gsub("+", "And", psite2(), fixed = T)
      } 
      paste0("LymphoAtlas_", na, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 8.5, height = 6.5)
      for (plot in plotInputExport()) {
       print(plot) 
      }
      dev.off()
    })
  # png output:
  output$Download1 <- downloadHandler(
    validate(
      need(!is.null(psite2()), "Select a site of interest to plot.")
    ),
    filename = function(){
      if (input$allSites) {
        na <- gsub(" / ", "-", input$protein, fixed = T)
      } else {
        na <- gsub("+", "And", psite2(), fixed = T)
      } 
      paste0("LymphoAtlas_", na, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = 800, res = 120)
      }
      ggsave(file, plot = plotInput(), device = device)
    })
  # svg output:
  output$Download2 <- downloadHandler(
    validate(
      need(!is.null(psite2()), "Select a site of interest to plot.")
    ),
    filename = function(){
      if (input$allSites) {
        na <- gsub(" / ", "-", input$protein, fixed = T)
      } else {
        na <- gsub("+", "And", psite2(), fixed = T)
      } 
      paste0("LymphoAtlas_", na, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::svg(..., width = 8.5, height = 6.5)
      }
      ggsave(file, plot = plotInput(), device = device)
    })
  ############
}


# Run the application 
shinyApp(ui = ui, server = server)

