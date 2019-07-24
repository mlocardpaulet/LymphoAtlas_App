##############################################################################################################
# Copyright CNRS 2019
# Contributors : Marie Locard-Paulet [marie.locard@ipbs.fr] and Guillaume Voisinne [voisinne@ciml.univ-mrs.fr]
# This software is a computer program whose purpose is to visualize and inspect the LymphoAtlas 
# phosphoproteomics data.
# This software is governed by the CeCILL license under French law and abiding by the rules 
# of distribution of free software. You can  use, modify and/or redistribute the software under 
# the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the 
# following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute 
# granted by the license, users are provided only with a limited warranty and the software's author, 
# the holder of the economic rights, and the successive licensors have only limited liability. 
# In this respect, the user's attention is drawn to the risks associated with loading, using, 
# modifying and/or developing or reproducing the software by the user in light of its specific 
# status of free software,that may mean that it is complicated to manipulate, and that also 
# therefore means that it is reserved for developers and experienced professionals having in-depth 
# computer knowledge. Users are therefore encouraged to load and test the software's suitability as 
# regards their requirements in conditions enabling the security of their systems and/or data to be 
# ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license 
# and that you accept its terms.
##########################################################################################################3#

library(shiny)
library(shinyBS)
library(shinydashboard)

library(ggplot2)
library(ggrepel)
library(data.table)
library(DT)
library(plotly)
library(heatmaply)
library(RColorBrewer)
library(viridis)

#devtools::install_github("VoisinneG/pannot")
library(pannot)

load("./data/df_merge.rda")

df_merge$Gene <- gsub("|", ";", df_merge$Gene, fixed = TRUE)

names(df_merge)[names(df_merge)=="Kinase_reported_mouse_human"] <- "Kinase-substrate"
names(df_merge)[names(df_merge)=="Protein.families"] <- "Protein families"

levels(df_merge$Cluster) <- c( levels(df_merge$Cluster) , "not regulated" )
df_merge$Cluster[is.na(df_merge$Cluster)] <- "not regulated"
df_reg <- df_merge[df_merge$Cluster != "not regulated", ]


########################################################################################################
# Mapping conditions

idx_intensity <- grep("^MeanLoops_", names(df_merge))
names_intensity <- names(df_merge)[idx_intensity]
s<-strsplit(as.character(names_intensity), split="_")

time <- sapply(s, function(x){x[4]})
time <- factor(time, levels = c("NS.", "S15.", "S30.", "S120.", "S300.", "S600."))
levels(time) <- c("0s","15s","30s","120s","300s","600s")
replicate <- sapply(s, function(x){x[2]})
dataset <- sapply(s, function(x){x[3]})

df_cond <- data.frame(name = names_intensity, time = time, replicate = replicate, dataset = dataset)


########################################################################################################
# Preprocessing data (averaging over technical replicates)

df_melt <- melt(df_merge, id.vars = c("psiteID"), measure.vars = names_intensity)
df_melt$time <- df_cond$time[match(df_melt$variable, df_cond$name)]
df_melt$replicate <- df_cond$replicate[match(df_melt$variable, df_cond$name)]

df_cast <- dcast(df_melt, psiteID +  time ~ replicate, mean , na.rm = TRUE)

df_melt1 <- melt(df_cast, id.vars = c("psiteID", "time"))
names(df_melt1)[names(df_melt1)=="variable"] <- "replicate"
df_melt1$Cluster <- df_merge$Cluster[match(df_melt1$psiteID, df_merge$psiteID)]
df_melt1$label <- df_merge$GeneID[match(df_melt1$psiteID, df_merge$psiteID)]


########################################################################################################
# App parameters

colClusters <- colorRampPalette(c("#9E0142", "#D53E4F", "#F46D43", 
                                  "#FDAE61", "#ABDDA4", "#66C2A5", 
                                  "#3288BD", "#5E4FA2", "darkblue"))(length(levels(df_merge$Cluster)) -1)
colClusters <- c(colClusters, "#FFFFFF")
names(colClusters) <- levels(df_reg$Cluster)
delim <- ";"
var_choices <- c("psiteID", "GeneID", "Gene", "Accession",  "Cluster", "Residue",
                 "Keywords", "Protein families", "Kinase-substrate",
                 "GO", "GO(biological process)", "GO(molecular function)", "GO(cellular component)", "Sequence")

var_display <- c(var_choices, "Sequence")

########################################################################################################
#User interface ----

body2<- dashboardBody(
  fluidRow(
    
  
  column(width = 6, height = NULL,
         tabBox(title = "t-SNE plot",
             width = NULL,
             tabPanel("t-SNE",
                      plotlyOutput("tsne", height = "500px")
                      ),
             tabPanel("Plot options",
                      checkboxInput("show_bckg", "Show all TCR-regulated phospho-sites", value = TRUE)
                      )
         ), 
         box(
           title = "Selection", 
           width = NULL, 
           height = NULL,
           selectizeInput("col_data_selection",
                          label = "Select columns to display",
                          choices = var_choices,
                          selected = c("GeneID", "Gene", "Cluster"),
                          multiple = TRUE),
           div(style = 'overflow-x: scroll', DT::dataTableOutput("selection"))
         )
  ),
  column(width = 6,height = NULL,
         tabBox(title = "Heatmap plot",
             width = NULL,
             tabPanel("Heatmap",
                      plotlyOutput("heatmap", height = "500px")),
             tabPanel("Plot options",
                      checkboxInput("scale", "scale values by phospho-site", value = TRUE),
                      conditionalPanel(condition = "input.scale",
                                       numericInput("max_scale", label = "scale limit", 2),
                                       br()
                                       )
             
             )
             
         ),
         tabBox(
           title = "Focus (single phospho-site)", 
           id = "tabset1",
           width = NULL, height = NULL,
           tabPanel("Plot", 
                    plotlyOutput("plot_focus", height = "300px")
           ),
           tabPanel("Plot options", 
                    checkboxInput("scale_focus", label = "Scale values", value = FALSE),
                    checkboxInput("boxplot_focus", label = "Show box plots", value = FALSE),
                    checkboxInput("scale_x_axis", label = "Scaled x axis (equal space between time points)", value = TRUE)
                    
           ),
           tabPanel("Info", 
                    selectizeInput("col_selected",
                                   label = "Select columns to display",
                                   choices = var_choices,
                                   selected = c("GeneID", "Keywords", "Protein families"),
                                   multiple = TRUE),
                    div(style = 'overflow-x: scroll', DT::dataTableOutput("data_focus"))
           )
         ),
         box(title = "General Info", 
             solidHeader = TRUE,
             width = NULL, 
             height = NULL,
             HTML(
                '<div> 
                LymphoAtlas is an open-source project released under the 
                <a href="https://www.cecill.info/index.en.html"> 
                  CeCILL license 
                </a>.
                Phosphoproteomics analysis was conducted on effector CD4+ T cells following TCR activation
                using cross-linking of anti-CD3 and anti-CD4 antibodies.
                Please see 
                <a href="a great paper here"> 
                  Locard-Paulet et al.
                </a>
                for more details.
                </div>'
                
              )
           
         )
  )
  ),
  HTML('<footer>
       <font size="0.8">copyright 2019 - CNRS - All rights reserved - LymphoAtlas V1.0 -</font> 
       <a href =  "https://github.com/mlocardpaulet/LymphoAtlas_App"> <img src="GitHub-Mark-32px.png", height = 20></a>
       </footer>')
  
)



sidebar <- dashboardSidebar(
  sidebarMenu(id = "menu",
              menuItem("Select",
                       tabName = "dashboard", 
                       startExpanded = TRUE,
                       icon = icon("check-circle"),
                       selectInput("var",
                                   label = "Variable",
                                   choices = var_choices,
                                   selected = "Keywords"),
                       bsTooltip("var",
                                 "Choose a selection variable",
                                 placement = "right"),
                       selectizeInput("selected",
                                      label = "Selected value(s)",
                                      choices = list(),
                                      selected = NULL,
                                      multiple = TRUE),
                       bsTooltip("selected",
                                 "Choose values corresponding to the selection variable above. Multiple choices are allowed.",
                                 placement = "right"),
                       checkboxInput("reg_only", "TCR-regulated sites only", value = TRUE),
                       bsTooltip("reg_only",
                                 "Select only phospho-sites whose phosphorylation significantly varies during the first 10 min following TCR-stimulation",
                                 placement = "right"),
                       br()
              ),
              menuItem("Annotations",
                       tabName = "annot",
                       startExpanded = FALSE,
                       icon = icon("check-circle"),
                       br(),
                       actionButton("update_annot", label = "Update annotations"),
                       bsTooltip("update_annot",
                                 "Retrieve up-to-date protein annotations from UniProt",
                                 placement = "right"),
                       br(),
                       HTML('<div> <font size="2"> Warning : This may take a long time. </font> </div>'),
                       br(),
                       textOutput("status_update"),
                       br()
                      # bsTooltip(id = "update_annot", title = "Warning : This may take a\n long time or even fail.\n Please be patient!")
              )
             
              
  )
  
  
)

ui <- dashboardPage(
  dashboardHeader(title = "LymphoAtlas"),
  sidebar,
  body2
)


########################################################################################################
# Server logic ----
server <- function(session, input, output) {

  cdata <- session$clientData
  
  react_val <- reactiveValues(psiteID_selected = NULL,
                              psiteID_heatmap = NULL,
                              psiteID_focus = NULL,
                              hover_psiteID = NULL,
                              hover_label = NULL,
                              hover_focus_x = NULL,
                              hover_focus_y = NULL,
                              delim = NULL,
                              value_selected = NULL,
                              data_merge = df_merge,
                              status = "Please be patient!")
  
  #####################################################################################################
  # Reactive functions
  {

    df_melt_selected <- reactive({
      df <- df_melt1
      df <- df[(df$psiteID %in% react_val$psiteID_selected) | (df$Cluster != "not regulated"), ]
      df
    })
      
    gtab_selection <- reactive({
      
      validate(
        need( length(react_val$psiteID_heatmap)>0, "Empty selection. Please modify selection choices." )
      )
      
      df <- df_melt_selected()
      df <- df[df$psiteID %in% react_val$psiteID_heatmap, ]
      df <- df[!is.na(df$value), ]
      
      validate(
        need( dim(df)[1]>0, "No data to display. Please modify selection choices." )
      )
      
      df$time_replicate <- paste(as.character(df$time), as.character(df$replicate), sep  = "-")
      lev <- levels(df$time)
      new_lev <- NULL
      for(i in 1:length(lev)){
        new_lev <- c( new_lev, paste(lev[i], levels(df$replicate), sep  = "-"))
      }
      levels(df$time_replicate) <- new_lev
      
      df_cast <- dcast(df, psiteID  ~ time_replicate, mean , na.rm = TRUE)
      df_cast <- df_cast[ , c("psiteID", levels(df$time_replicate)[levels(df$time_replicate) %in% names(df_cast)] )]
      
      df_cast

    })
    
    gtab_selection_scaled <- reactive({
      
      df_cast <- gtab_selection()
      
      #####################################################################################################
      # Scale data
      
      if(input$scale){
        df_scale <- as.data.frame(t(scale(t(df_cast[, -1]))))
      }else{
        df_scale <- df_cast[, -1]
      }
      df_scale$psiteID <- df_cast$psiteID
      df_cast <- df_scale
      
      #####################################################################################################
      
      idx_match <- match(df_cast$psiteID, react_val$data_merge$psiteID)
      df_cast[["Cluster"]] <- react_val$data_merge[["Cluster"]][idx_match]
      df_cast[["GeneID"]] <- react_val$data_merge[["GeneID"]][idx_match]
      df_cast <- df_cast[order(df_cast$Cluster, df_cast$GeneID), ]
      
      df_cast
      
    })
    
    
    
    gtab_focus <- reactive({
      
      validate(
        need( length(react_val$psiteID_focus) > 0, "Empty selection. Please select a phospho-site." )
      )
      
      df <- df_melt_selected()
      df <- df[df$psiteID == react_val$psiteID_focus, ]
      df <- df[!is.na(df$value), ]
      
      
      validate(
        need( dim(df)[1] > 0, "No data to display for this phospho-site" )
      )
      
      if(input$scale_focus){
        df$value <- scale(df$value)
      }
      
      df_cast <- dcast(df, psiteID ~ time, mean, na.rm = TRUE)
      df_melt <- melt(df_cast, id.vars = c("psiteID"))
      names(df_melt)[which(names(df_melt)=="variable")] <- "time"
      df_melt$variable <- "mean"
      df_melt
      
    })
    
    data_selection <- reactive({
      validate(  
        need( length(react_val$psiteID_heatmap)>0, "No phospho-site selected" )
      )
      
      idx_match <- match(react_val$psiteID_heatmap, react_val$data_merge$psiteID)
      df <- react_val$data_merge[idx_match, ]
      
      df <- df[order(df$Cluster, df$GeneID), ]
      
      react_val$psiteID_focus <- df$psiteID[1]
      
      df
      
    })
  }

  
  #####################################################################################################
  # Observe environments
  {
    observe({
      
      react_val$delim <- switch(input$var,
                                "Protein families" = ", ",
                                "Keywords" =";",
                                "Kinase-substrate" = ";",
                                "Gene" = ";",
                                "; ")
    })
    
    observeEvent(c(input$var, input$reg_only), {
      if(input$reg_only){
        all_levels <- paste(as.character(react_val$data_merge[[input$var]][react_val$data_merge$Cluster != "not regulated"]), collapse=react_val$delim)
      }else{
        all_levels <- paste(as.character(react_val$data_merge[[input$var]]), collapse=react_val$delim)
      }
      
      
        
      unique_levels <- unique(strsplit(all_levels, split = react_val$delim)[[1]])
      if(input$var == "Gene"){
        
        unique_levels <- unique(strsplit(paste(unique_levels, collapse = " "), split=" ", fixed = TRUE)[[1]])
      }
      
      unique_levels <- unique_levels[order(unique_levels)]
      
      #print(unique_levels)
      
      default_selection <- switch(input$var,
                                  "psiteID" = "Q9Z0R6_Y554",
                                  "GeneID" = "Itsn2_Y554",
                                  "Accession" = "Q9Z0R6",
                                  "Gene" = "Itsn2",
                                  "Cluster" = 7,
                                  "Keywords" = "Actin-binding",
                                  "Protein families" = "Protein kinase superfamily",
                                  "Kinase-substrate" = "Akt1",
                                  "Residue" = "Y", 
                                  "GO" = "cytoplasmic vesicle [GO:0031410]",
                                  "GO(biological process)" = "endocytosis [GO:0006897]",
                                  "GO(molecular function)" = "calcium channel regulator activity [GO:0005246]",
                                  "GO(cellular component)" = "immunological synapse [GO:0001772]",
                                  NULL
      )
      
      updateSelectInput(session, "selected",
                          choices = unique_levels,
                          selected = default_selection)

      
    })
    
    observe({
      if( ! input$var %in% c("GO", 
                            "GO(biological process)", 
                            "GO(molecular function)", 
                            "GO(cellular component)",
                            "Gene")){
        idx_match <-unlist(
          lapply(input$selected, function(x){
            grep(paste("(",react_val$delim,"|^)", x,"($|", react_val$delim, ")", sep=""),
                 as.character(react_val$data_merge[[input$var]]), fixed=FALSE)
          })
        )
      }else if(input$var == "Gene"){
        idx_match <-unlist(
          lapply(input$selected, function(x){
            grep(paste("(",react_val$delim,"|\\s|^)", x,"($|\\s|", react_val$delim, ")", sep=""),
                 as.character(react_val$data_merge[[input$var]]), fixed=FALSE)
          })
        )
        print(idx_match)
      }else{
        idx_match<-unlist(
          lapply(input$selected, function(x){
            grep(x, as.character(react_val$data_merge[[input$var]]), fixed=TRUE)
          })
        )
      }
      
      df <- react_val$data_merge[idx_match, ]
      
      if(input$reg_only){
        react_val$psiteID_selected <- df$psiteID[df$Cluster != "not regulated"]
      }else{
        react_val$psiteID_selected <- df$psiteID
      }
      react_val$psiteID_heatmap <- react_val$psiteID_selected
      
    })
    
    observe({
      df <- gtab_selection_scaled()
      event.data <- event_data("plotly_click", source = "select_heatmap")
      idx <- dim(df)[1] - event.data$pointNumber[[1]][1]
      if(length(idx)>0){
        react_val$psiteID_focus <- df$psiteID[idx]
      }
    })
    
    observeEvent(event_data("plotly_selected", source = "select_tsne"), {
      
      event.data <- event_data("plotly_selected", source = "select_tsne")
      cluster_displayed <- unique(df_reg$Cluster[match(react_val$psiteID_selected, df_reg$psiteID)])
      cluster_displayed <- cluster_displayed[!is.na(cluster_displayed)]
      cluster_displayed <- cluster_displayed[order(cluster_displayed, decreasing = FALSE)]
      
      cat(cluster_displayed)
      
      if(!is.null(event.data) ){
        if(length(event.data$curveNumber)>0){
        
        if(!input$show_bckg){
          event.data$curveNumber <- event.data$curveNumber + 1
        }
        
        df_tsne <- df_reg
        df_tsne$Cluster <- as.character(df_tsne$Cluster)
        idx_bckg <- which(! as.character(df_tsne$psiteID) %in% react_val$psiteID_selected )
        
        if(length(idx_bckg)>0){
          df_tsne$Cluster[idx_bckg] <- 0
        }
        
        psiteID <- rep(NA, dim(event.data)[1])
        for(i in 1:dim(event.data)[1]){
          if(event.data$curveNumber[i] > 0){
            idx <- which(df_tsne$Cluster == cluster_displayed[event.data$curveNumber[i]] )
            psiteID[i] <- as.character(df_tsne$psiteID[idx[event.data$pointNumber[i] + 1]])
          }else{
            psiteID[i] <- as.character(df_tsne$psiteID[event.data$pointNumber[i] + 1])
          }
          
        }
        
        react_val$psiteID_heatmap <- psiteID
        }else{
          react_val$psiteID_heatmap <- react_val$psiteID_selected
        }
      }else{
        react_val$psiteID_heatmap <- react_val$psiteID_selected
      }
      
    })
    
    observe({
      event.data <- event_data("plotly_click", source = "select_tsne")
      if(!is.null(event.data)){
        idx_min <- which.min( (event.data$x - df_reg$X)^2 + (event.data$y - df_reg$Y)^2 )
        if(length(idx_min)>0){
          react_val$psiteID_focus <- df_reg$psiteID[idx_min]
        }else{
          react_val$psiteID_heatmap <- react_val$psiteID_selected
        }
      }else{
        react_val$psiteID_heatmap <- react_val$psiteID_selected
      }

    })

    observe({
      if(length(input$selection_rows_selected)>0){
        react_val$psiteID_focus <- data_selection()$psiteID[input$selection_rows_selected]
      }
    })
    
    
    
    observeEvent(input$update_annot,{
      
        # Create a Progress object
        react_val$status <- "Querying Uniprot..."
        
        progress <- shiny::Progress$new(min = 0, max = 100)
        progress$set(message = "Querying UniProt...", value = 0)
        on.exit(progress$close())
        updateProgress <- function(value = NULL, detail = NULL) {
          progress$set(value = value, detail = detail)
        }
        
        df_annot <- pannot::get_annotations_uniprot(id = react_val$data_merge$Entry,
                                                     max_keys = 400,
                                                     columns = c("id",
                                                                 "keywords",
                                                                 "families",
                                                                 "go",
                                                                 "go(biological_process)",
                                                                 "go(molecular_function)",
                                                                 "go(cellular_component)"),
                                                     updateProgress = updateProgress)
                                           
        if(!is.null(df_annot)){
          idx_match <- match(react_val$data_merge$Entry, df_annot$query_id)
          react_val$data_merge[["GO"]] <- df_annot[["Gene.ontology..GO."]][idx_match]
          react_val$data_merge[["GO(biological process)"]] <- df_annot[["Gene.ontology..biological.process."]][idx_match]
          react_val$data_merge[["GO(molecular function)"]] <- df_annot[["Gene.ontology..molecular.function."]][idx_match]
          react_val$data_merge[["GO(cellular component)"]] <- df_annot[["Gene.ontology..cellular.component."]][idx_match]
          react_val$data_merge[["Protein families"]] <- df_annot[["Protein.families"]][idx_match]
          react_val$data_merge[["Keywords"]] <- df_annot[["Keywords"]][idx_match]
          
          react_val$status <- "Done!"
          
        } else{
          react_val$status <- "Query failed..."
        }                                         

        
      
    })
    
  }
  
  #####################################################################################################
  # Output text
  
  output$status_update <- renderText({
    react_val$status
  })
  
  #####################################################################################################
  # Output data tables
  {
    output$data_focus <- DT::renderDataTable({
      validate(
        need( length(input$col_selected)>0, "Please select data columns to display" )
      )
      validate(
        need( length(react_val$psiteID_focus)>0, "No phospho-site selected" )
      )
      
      DT::datatable(
        react_val$data_merge[match(react_val$psiteID_focus, react_val$data_merge$psiteID), input$col_selected],
        rownames = FALSE
      )
    })
    
    output$selection <- DT::renderDataTable({
      
      validate(
        need( length(input$col_data_selection)>0, "Please select data columns to display" )
      )

      DT::datatable(
        data_selection()[ , input$col_data_selection],
        selection = list(mode = "single", target = "row"),
        rownames = FALSE
      )

    })
    
  }
  
  
  #####################################################################################################
  # Plots 
  {
    output$tsne <- renderPlotly({
      p <- ggplot(df_reg, aes(x=X, y=Y, label = GeneID))
      
      if(input$show_bckg){
        p <- p + geom_point(alpha = 0.1, size = 2, color = "black")
      }

      p <- p +
        scale_color_manual(values = colClusters) +
        theme(title = element_text(size = 7)) +
        ggtitle(paste(input$var, ": ",paste(input$selected, collapse = "+\n"), sep ="")) +
        xlab("t-SNE 1") + 
        ylab("t-SNE 2") +
        geom_point(df_reg[df_reg$psiteID %in% react_val$psiteID_selected,],
                   mapping = aes(x=X, y=Y, color=Cluster, label = GeneID),
                   size = 4,
                   alpha = 0.8,
                   inherit.aes = FALSE) +
        coord_cartesian(xlim = range(df_reg$X), ylim = range(df_reg$Y), expand = TRUE)
      
      
      s<-ggplotly(p, tooltip = "label") %>% layout(dragmode = "select")
      s$x$source <- "select_tsne"
      s
    })
    
    output$heatmap <- renderPlotly({
      
      df <- gtab_selection_scaled()
      df_plot <- df[ , - which(names(df) %in% c("psiteID", "GeneID", "Cluster")) ]
      
      row_labels <- df$GeneID
      psite_label_size <- 8
      if(dim(df)[1]<30){
        psite_label_size <- 12
      }else if(dim(df)[1]<40){
        psite_label_size <- 10
      }else if(dim(df)[1]>60){
        psite_label_size <- 1
      }
      
      if(input$scale){
        df_plot[df_plot > input$max_scale] <- input$max_scale
        df_plot[df_plot < -input$max_scale] <- -input$max_scale
        zlims <- c(-input$max_scale, input$max_scale)
        main_title <- "Scaled values (Z-score)"
      }else{
        zlims <- NULL
        main_title <- "log2-transformed Intensities"
      }
      
      rsc <- data.frame("Cluster" = df$Cluster)
      
      p <- heatmaply(df_plot, 
                     labRow = row_labels,
                     Rowv = FALSE, 
                     Colv = FALSE, 
                     colors = viridis,
                     row_side_colors = rsc,
                     row_side_palette = colClusters[df$Cluster],
                     margins = c(80, 80, 50, 0),
                     main = main_title,
                     plot_method="plotly",
                     column_text_angle = 90,
                     colorbar_yanchor = "middle",
                     fontsize_row = psite_label_size,
                     limits = zlims,
                     xlab=("time / replicate"),
                     ylab = "phospho-sites (GeneID)")
      
      p$x$data[[2]]$showscale <- FALSE
      p$x$data[[1]]$colorbar$yanchor <- "middle"
      p$x$data[[1]]$colorbar$y <- 0.5
      p$x$data[[1]]$colorbar$xanchor <- "middle"
      p$x$data[[1]]$colorbar$x <- 1
      p$x$source <- "select_heatmap"
      
      
      # initiate a line shape object
      line <- list(
        type = "line",
        line = list(color = "white", width = 5),
        xref = "x",
        yref = "y"
      )
      
      lines <- list()
      idx_rep <- sapply(levels(time), 
             function(x){ 
               idx <- grep(paste("^", x, "-",sep=""), colnames(df_plot))
               ifelse(length(idx)>0, idx[length(idx)], NULL)
               })

      for (i in idx_rep) {
        line[c("x0", "x1")] <- i + 0.5
        line[["y0"]] <- 0
        line[["y1"]] <- length(df$psiteID) + 0.5
        lines <- c(lines, list(line))
      }
      
      p <- layout(p, shapes = lines)
      
      p
      
    })
    
    output$plot_focus <- renderPlotly({
      
      validate(
        need( length(react_val$psiteID_focus) > 0, "Empty selection. Please select a phospho-site." )
      )
      
      df <- df_melt_selected()
      df <- df[df$psiteID == react_val$psiteID_focus, ]
      
      if(input$scale_focus){
        df$value <- scale(df$value)
        ylabel <- "scaled values (Z score)"
      }else{
        ylabel <- "Intensity (log2)"
      }
      
      df_mean <- gtab_focus()
      
      xlabel <- "time"
      if(!input$scale_x_axis){
        df$time <- as.numeric( do.call(rbind, strsplit(as.character(df$time), split="s")) )
        df_mean$time <- as.numeric( do.call(rbind, strsplit(as.character(df_mean$time), split="s")) )
        xlabel <- "time (s)"
      }
      
      p <- ggplot(df, aes(x=time, y=value, color = replicate)) + 
        #theme(legend.title=element_blank()) +
        ggtitle(react_val$data_merge$GeneID[ match(react_val$psiteID_focus, react_val$data_merge$psiteID) ])
        
      if(input$boxplot_focus){
        p <- p + geom_boxplot(mapping=aes(x=time, y=value), inherit.aes = FALSE, color = rgb(0.75, 0.75, 0.75) )
      }
      
      p <- p + ylab( ylabel ) +
        xlab(xlabel) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(mapping = aes(group = replicate), alpha = 0.5) +
        geom_line(data = df_mean,
                  mapping = aes(x=time, y=value, color = variable, group = 1),
                  inherit.aes = FALSE,
                  alpha = 0.5) +
        scale_color_discrete(name = "replicate")
      
      
      ggplotly(p)
      
    })
  }
  
}

shinyApp(ui, server)


