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

# Colour scales for the different clusters:
# coloursLines <- c("blue4", "deepskyblue4", "lightseagreen", "mediumspringgreen")
colClusters <- colorRampPalette(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "darkblue"))(length(unique(export$ClusterMerged)))
collight <- as.character(
  sapply(colClusters, function(x) {
    colorRampPalette(c("white", x))(10)[5]
  })
)
matriceColours <- cbind("dark" = colClusters, "light" = collight, "cluster" = 1:length(unique(export$ClusterMerged)))
matriceColours <- rbind(matriceColours, c("grey30", "grey70", "NA"))

colorcodes <- paste0(rep(matriceColours[,3], each = 4), rep(c(" R1", " R3", " R4", " R5"), nrow(matriceColours)))
colorcodes <- paste0("Cluster ", colorcodes)
colorcodes <- gsub("Cluster NA", "Not regulated", colorcodes, fixed = T)
i <- 1
j <- 1
colrep <- vector()
while (i <= nrow(matriceColours)) {
  colrep[j:(j+3)] <- colorRampPalette(c(matriceColours[i,2], matriceColours[i,1]))(4)
  i <- i + 1
  j <- j + 4
}
matriceColours2 <- data.frame("colorcodes" = colorcodes, "value" = colrep)


# App:
############################################################################
## UI:
############################################################################
ui <- fluidPage(
  fluidRow(
    column(3, titlePanel("LymphoAtlas")
    ),
    column(1, actionButton(inputId='ab1', label="?", 
                           onclick ="window.open('https://masstools.ipbs.fr/lymphoAtlasHelp.html', '_blank')",
                           style="color: #fff; background-color: #673a49; border-color: #000000")
    )
  ),
  tags$style(type='text/css', "#ab1 { width:80%; margin-top: 25px; font-family : Cursive; font-weight: 900; font-size: 160%;}"),
  
  # Sidebar with a slider input for the phosphorylation site of interest: 
  
  radioButtons("SearchMode", "Search mode:",
               c("by identification" = 'ID',
                 "by kinetics" = 'kinetics'),
               selected = 'ID',
               inline = TRUE
  ), 
  sidebarLayout(
    sidebarPanel(width = 4,
                 checkboxInput("Scalex", "Scaled x-axis (same space between time points)", FALSE), # To switch to factors on the x-scale (help visualising the early time points -15 and 30sec-)
                 # Conditional panels:
                 # Part 1:
                 conditionalPanel(condition="input.SearchMode== 'ID'",
                                  # Select the protein of interest: ---
                                  selectizeInput("protein",
                                                 "Protein of interest:",
                                                 choices = sort(unique(paste0(export$GeneNames, " / ", export$Accession))),
                                                 selected = "Zap70 / P43404",
                                                 multiple = FALSE),
                                  bsTooltip("protein", 
                                            "If the list of proteins to select does not include your protein of interest, enter it in the selection box.",
                                            "right"),
                                  # Once the protein of interest is selected, select one or all the sites of the protein: ---
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
                                  conditionalPanel(condition = "input.allSites==true",
                                                   checkboxInput("fixedAxis",
                                                                 "Free the y-axis",
                                                                 FALSE),
                                                   bsTooltip("fixedAxis",
                                                             "Check if you want the y-axis to be free",
                                                             "right")
                                  )
                 ),
                 conditionalPanel(condition="input.SearchMode== 'kinetics'",
                                  # Select the cluster of interest: ---
                                  selectizeInput("cluster",
                                                 "Cluster of interest:",
                                                 choices = as.character(sort(unique(export$ClusterMerged))),
                                                 selected = "1",
                                                 multiple = FALSE),
                                  selectizeInput("psite2",
                                                 "Phosphorylation site of interest:",
                                                 choices = NULL,
                                                 selected = NULL,
                                                 multiple = FALSE)
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
                           "right"
                 )
    ),
    
    # Show a plot
    mainPanel(width = 8,
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
  # 
  ChekAllSites <- reactive({
    if (input$SearchMode == "kinetics") {
      return(FALSE)
    } else {
      return(input$allSites)
    }
  })
  # Define protein: ---
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
  # Update the field of the psite selection: ---
  observe({
    if (input$SearchMode == "ID") {
      updateSelectizeInput(session, "psite",
                           "Phosphorylation site of interest:",
                           choices = sort(as.character(export$GeneID[export$Accession == prot2()])),
                           selected = NULL)
    } else {
      updateSelectizeInput(session, "psite2",
                           "Phosphorylation site of interest:",
                           choices = sort(as.character(export$GeneID[as.character(export$ClusterMerged) %in% input$cluster])),
                           selected = NULL)
    }
  })
  # If the checkbox "allSites" is checked, plot all the sites of the protein, else, plot only the one selected in the field "psite": ---
  psite2 <- reactive({
    if (input$SearchMode == "ID") {
      if (ChekAllSites()) {
        res <- export$phosphoSites[export$Accession == prot2()]
        return(res)
      } else {
        if (is.null(input$psite)) {
          return(NULL)
        } else {
          return(input$psite)
        }
      }
    } else {
      if (is.null(input$psite2)) {
        return(NULL)
      } else {
        return(input$psite2)
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
      if (ChekAllSites()) {
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
        gtab$Cluster <- as.character(export$ClusterMerged[match(gtab$phosphosite, export$GeneID)])
        gtab$colorgroup <- paste0("Cluster ", gtab$Cluster, " ", gtab$Replicate)
        gtab$colorgroup <- gsub("Cluster NA", "Not regulated", gtab$colorgroup, fixed = T)
        #gtab <- gtab[!is.na(gtab$value),]
        gtab <- gtab[order(gtab$variable),]
        gtab <- gtab[order(gtab$phosphosite),]
        # Remove chunk of NA data:
        i <- 1
        vec <- gtab$value
        vec[is.na(vec)] <- 0
        toremove <- vector()
        while (i <= nrow(gtab)) {
          val <- sum(vec[i:(i+5)])
          if (val == 0) {
            toremove <- c(toremove, c(i:(i+5)))
          }
          i <- i + 6
        }
        gtab <- gtab[-toremove,]
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          # Table preparation:
          cols <- which(grepl("MeanLoops", names(export)))
          gtab <- export[export$GeneID %in% psite2(),cols]
          if (nrow(gtab) > 1) {
            gtab <- gtab[order(rowSums(gtab, na.rm = T), decreasing = T),][1,]
            gtab$phosphosite <- unique(export$GeneID[export$GeneID %in% psite2()])
          } else {
            gtab$phosphosite <- export$GeneID[export$GeneID %in% psite2()]
          }
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
          gtab$Cluster <- as.character(export$ClusterMerged[match(gtab$phosphosite, export$GeneID)])
          gtab$colorgroup <-  paste0("Cluster ", gtab$Cluster, " ", gtab$Replicate)
          gtab$colorgroup <- gsub("Cluster NA", "Not regulated", gtab$colorgroup, fixed = T)
          #gtab <- gtab[!is.na(gtab$value),]
          gtab <- gtab[order(gtab$variable),]
          gtab <- gtab[order(gtab$phosphosite),]
          # Remove chunk of NA data:
          i <- 1
          vec <- gtab$value
          vec[is.na(vec)] <- 0
          toremove <- vector()
          while (i <= nrow(gtab)) {
            val <- sum(vec[i:(i+5)])
            if (val == 0) {
              toremove <- c(toremove, c(i:(i+5)))
            }
            i <- i + 6
          }
          gtab <- gtab[-toremove,]
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
      if (input$Scalex) {
        # If allSites is checked: plot all the sites for the protein selected:
        #########################
        if (ChekAllSites()) {
          g <- ggplot(gtab, aes(x = factor(TimePoint), y = value)) +
            geom_line(aes(x = factor(TimePoint), y = value, group = Replicate, col = colorgroup), size = 1.2) +
            geom_vline(xintercept = 0, size = 1.2) +
            geom_point(col = "black", shape = "+", size = 4, alpha = 0.8) +
            theme_minimal() + ggtitle(paste0("Sites of ", input$protein, " upon TCR activation")) +
            scale_color_manual(values = as.character(matriceColours2$value[match(sort(unique(gtab$colorgroup)), matriceColours2$colorcodes)]), name = "Cluster and\nbiological repeat") +
            ylab("log2-transformed normalised MS intensities") +
            xlab("Time after stimulation (in seconds)")
          if (input$fixedAxis == FALSE) {
            g <- g + facet_wrap(~phosphosite) 
          } else {
            g <- g + facet_wrap(~phosphosite, scales = "free_y", shrink = FALSE) 
          }
          
          # If allSites is notchecked: plot all the sites for the protein selected:
          #########################
        } else {
          if (is.null(psite2())) {
            return(NULL)
          } else {
            g <- ggplot(gtab, aes(x = factor(TimePoint), y = value)) +
              geom_line(aes(x = factor(TimePoint), y = value, group = Replicate, col = Replicate), size = 1.2) +
              geom_vline(xintercept = 0, size = 1.2) +
              geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) +
              theme_minimal() + ggtitle(paste0(psite2(), " upon TCR activation")) +
              scale_color_manual(values = as.character(unique(matriceColours2$value[match(gtab$colorgroup, matriceColours2$colorcodes)])), name = "Replicate") + 
              ylab("log2-transformed normalised MS intensities") +
              xlab("Time after stimulation (in seconds)") 
          }
        }
      } else {
        # If allSites is checked: plot all the sites for the protein selected:
        #########################
        if (ChekAllSites()) {
          g <- ggplot(gtab, aes(x = TimePoint, y = value)) +
            geom_line(aes(x = TimePoint, y = value, group = Replicate, col = colorgroup), size = 1.2) +
            geom_vline(xintercept = 0, size = 1.2) +
            geom_point(col = "black", shape = "+", size = 4, alpha = 0.8) +
            theme_minimal() + ggtitle(paste0("Sites of ", input$protein, " upon TCR activation")) +
            scale_color_manual(values = as.character(matriceColours2$value[match(sort(unique(gtab$colorgroup)), matriceColours2$colorcodes)]), name = "Cluster and\nbiological repeat") +
            ylab("log2-transformed normalised MS intensities") +
            xlab("Time after stimulation (in seconds)")
          if (input$fixedAxis == FALSE) {
            g <- g + facet_wrap(~phosphosite) 
          } else {
            g <- g + facet_wrap(~phosphosite, scales = "free_y", shrink = FALSE) 
          }
          
          # If allSites is notchecked: plot all the sites for the protein selected:
          #########################
        } else {
          if (is.null(psite2())) {
            return(NULL)
          } else {
            g <- ggplot(gtab, aes(x = TimePoint, y = value)) +
              geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) +
              geom_vline(xintercept = 0, size = 1.2) +
              geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) +
              theme_minimal() + ggtitle(paste0(psite2(), " upon TCR activation")) +
              scale_color_manual(values = as.character(unique(matriceColours2$value[match(gtab$colorgroup, matriceColours2$colorcodes)])), name = "Replicate") + 
              ylab("log2-transformed normalised MS intensities") +
              xlab("Time after stimulation (in seconds)") 
          }
        }
      }
      return(g)
    }
  }

  output$psiteplot <- renderPlotly({
    p <- plotInput()
    ggplotly(p, height = 600) %>%
      config(displayModeBar = F)%>%
      layout(margin = list(l = 110, b = 40, r = 40, t = 110, pad = -2))
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
      if (ChekAllSites()) {
        lg <- list()
        ylimplot <- range(gtab$value)
        for (el in unique(gtab$phosphosite)) {
          gtab2 <- gtab[gtab$phosphosite %in% el,]
          gtab2$Cluster[is.na(gtab2$Cluster)] <- "NA"
          
          g <- ggplot(gtab2, aes(x = TimePoint, y = value))  +
            geom_line(aes(x = TimePoint, y = value, group = colorgroup, col = colorgroup), size = 1.2) +
            geom_vline(xintercept = 0, size = 1.2) +
            geom_point(col = "black", shape = "+", size = 4, alpha = 0.8) +
            theme_minimal() + ggtitle(paste0(el, " upon TCR activation")) +
            ylab("log2-transformed normalised\nMS intensities") +
            xlab("Time after stimulation (in seconds)")
          if (input$fixedAxis == FALSE) {
            g <- g + ylim(ylimplot)
          } 
          # I choose the colours later.
          lg[[length(lg) + 1]] <- g
        }
        # If allSites is not checked: 
        #########################
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          g <- ggplot(gtab, aes(x = TimePoint, y = value)) +
            geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) +
            geom_vline(xintercept = 0, size = 1.2) +
            geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) +
            theme_minimal() +
            ggtitle(paste0(psite2(), " upon TCR activation")) +
            ylab("log2-transformed normalised\nMS intensities") +
            xlab("Time after stimulation (in seconds)")  # I choose the colours later.
          lg <- list(g)
        }
      }
      return(lg)
    }
  }
  # Colour scales for download:
  loopcolourvalue <- function(){
    if (is.null(plotTable())) {
      return(NULL)
    } else {
      gtab <- plotTable()
      # If allSites is checked: plot all the sites for the protein selected:
      #########################
      if (ChekAllSites()) {
        colourval <- list()
        for (el in unique(gtab$phosphosite)) {
          gtab2 <- gtab[gtab$phosphosite %in% el,]
          gtab2$Cluster[is.na(gtab2$Cluster)] <- "NA"
          colourval[[length(colourval) + 1]] <- as.character(unique(matriceColours2$value[match(sort(gtab2$colorgroup[gtab2$phosphosite %in% el]), matriceColours2$colorcodes)]))
        }
        # If allSites is not checked: 
        #########################
      } else {
        if (is.null(psite2())) {
          return(NULL)
        } else {
          colourval[[length(colourval) + 1]] <- as.character(unique(matriceColours2$value[match(sort(gtab$colorgroup), matriceColours2$colorcodes)]))
        }
      }
      return(colourval)
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
      pdf(file, width = 4.5, height = 3)
      if (ChekAllSites()) {
        for (iterplot in seq_along(plotInputExport())) {
          loopcolval <- loopcolourvalue()[[iterplot]]
          toplot <- plotInputExport()[[iterplot]]
          toplot <- toplot + scale_color_manual(values = loopcolval, name = "Cluster and\nbiological repeat") 
          print(toplot) 
        }
      } else {
        print(plotInput())
      }
      dev.off()
    })
  # png output:
  output$Download1 <- downloadHandler(
    validate(
      need(!is.null(psite2()), "Select a site of interest to plot.")
    ),
    filename = function(){
      if (ChekAllSites()) {
        na <- gsub(" / ", "-", input$protein, fixed = T)
      } else {
        na <- gsub("+", "And", psite2(), fixed = T)
      } 
      paste0("LymphoAtlas_", na, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = 800, res = 200)
      }
      ggsave(file, plot = plotInput(), device = device)
    })
  # svg output:
  output$Download2 <- downloadHandler(
    validate(
      need(!is.null(psite2()), "Select a site of interest to plot.")
    ),
    filename = function(){
      if (ChekAllSites()) {
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

