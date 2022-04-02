#
# This is a Shiny application. You can run the application by clicking
# the 'Run App' button above.
#
# Author: David Garcia Seisdeods
# Web: https://github.com/davco6
# 
#

#Libraries
library(shiny)
library(shinythemes)
library(dplyr)
library(readr)
library(stringi)
#source("https://neuroconductor.org/neurocLite.R")
# Default Install
#neuro_install('EveTemplate')

library("EveTemplate")
library(neurohcp)
library("data.table")
library(fslr)
setwd(".")
links <- fread("./dataset_links.txt")
dict_Onotology_terms <- fread("./Onology_terms_to_app_dic.txt")
labels <- getEveMapLabels(type="I")
map <- readEveMap(type="I")
eve_t1 <- readEve()

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                titlePanel("Protein expression in the brain"),
                sidebarLayout(
                  sidebarPanel(
                    
                    # Select type of trend to plot
                    selectInput(inputId = "experiment", label = strong("Dataset (Authorship)"),
                                choices = links$Authorship,
                                selected = "Mendonca CF et al."),
                    
                    uiOutput("selectgene"),
                    
                    uiOutput("selectcondition"),
                    
                    conditionalPanel(condition = "input.Disease == 'Alzheimer disease'",
                                     selectInput(inputId = "Stage",
                                                 label = strong("Select Braak stage"),
                                                 choices = c("Early braak stage", "Late braak stage"))),
                    
                    sliderInput(inputId = "x.axis", label = "x axis:",
                                min = 1, max = 181, value = 90, step = 1),
                    sliderInput(inputId = "y.axis", label = "y axis:",
                                min = 1, max = 217, value = 109, step = 1),
                    sliderInput(inputId = "z.axis", label = "z axis:",
                                min = 1, max = 181, value = 93, step = 1)
                  ),
                  # Output: Description, lineplot, and reference
                  mainPanel(
                    textOutput(outputId = "desc"),
                    textOutput(outputId = "desc2"),
                    textOutput(outputId = "desc3"),
                    plotOutput(outputId = "plotBrainRegions", height = "300px"),
                    tableOutput(outputId = "table"),
                    
                    plotOutput(outputId = "plotBrainExpression", height = "300px")
                    
                    
                    
                  )
                )
)

# Define server function
server <- function(input, output) {
  
  
  selected_exp_sdrf <- reactive({
    req(input$experiment)
    sdrf <- fread(links$sdrf_link[which(links$Authorship==input$experiment)], stringsAsFactors = TRUE)
    sdrf 
  }) 
  
  get_factor_names <- reactive({
    sdrf <- selected_exp_sdrf()
    factor_table <- colnames(sdrf)[grep("Factor Value\\[", colnames(sdrf))]
    factor_names <- sort(gsub("\\]", gsub("Factor Value\\[",factor_table, replacement = ""), replacement = ""))
    factor_names
  })
  
  output$selectgene <- renderUI({
    proteinGroups <- selected_exp_table()
    selectInput(inputId = "Gene",
                label = strong("Select gene"),
                selectize = TRUE,
                choices = proteinGroups$`Gene Name`)
  })
  
  selected_exp_table <- reactive({
    
    req(input$experiment)
    factor_names <- get_factor_names()
    
    proteinGroups <- as.data.frame(fread(links$ProteinGroups_link[which(links$Authorship==input$experiment)],
                                         stringsAsFactors = FALSE))
    proteinGroups[is.na(proteinGroups)] <- 0  
    
    if ("individual" %in% factor_names){
      index_to_remove <- which( factor_names == "individual")
      new_col_names <- lapply(strsplit(colnames(proteinGroups), ", " ), function(x) paste0(x[-index_to_remove], collapse = ", "))
      new_proteinGroups <- data.frame("Gene ID" = proteinGroups$`Gene ID`,
                                      "Gene Name"= proteinGroups$`Gene Name`,
                                      check.names = FALSE)
      for (i in unique(new_col_names)[c(-1,-2)]){
        index_individuals <- which(i == new_col_names)
        if (length(index_individuals)>1){
          new_proteinGroups[i]<-apply(proteinGroups[,index_individuals],1, mean)}
        else(
          new_proteinGroups[i]<-proteinGroups[,index_individuals]
        )
      }
      proteinGroups <- new_proteinGroups
    }
    proteinGroups
  })
  
  
  
  
  output$selectcondition <- renderUI({
    if ("disease" %in% get_factor_names()){
      if ("disease staging" %in% get_factor_names()){
        selectInput(inputId = "Disease", label = strong("Select condition"),
                    choices = c("Control", "Alzheimer disease"),
                    selected = "Control")
        
        
      }
      else {
        sdrf <- selected_exp_sdrf()
        if (dim(unique(sdrf[,'Factor Value[disease]']))[1]==1){
          selectInput(inputId = "Disease", label = strong("Select condition"),
                      choices = c("Control"),
                      selected = "Control")
        }
        else {
          selectInput(inputId = "Disease", label = strong("Select condition"),
                      choices = c("Control", "Alzheimer"),
                      selected = "Control")
        }
        
      }
    }
    else {
      selectInput(inputId = "Disease", label = strong("Select condition"),
                  choices = c("Control"),
                  selected = "Control")
    }
    
  })
  
  
  selected_query_gene_table <- reactive({
    req(input$Gene)
    req(input$experiment)
    validate(need(toupper(input$Gene) %in% selected_exp_table()$`Gene Name`, paste0("Error: Gene not found in ", input$experiment)))
    Row_out <- which(selected_exp_table()$`Gene Name`==toupper(input$Gene))
    max_value<-max(na.exclude(selected_exp_table()[Row_out, c(-1,-2)]))
    list(Row_out, max_value)
  })
  
  selected_query_conditions <- reactive({
    sdrf <- selected_exp_sdrf()
    factor_names <- get_factor_names()
    organism_parts <- levels(sdrf$`Sample Characteristic[organism part]`)
    tmp_organism_parts <- c()
    for ( i in organism_parts){
      if (i %in% dict_Onotology_terms$`Ontology term`){
        tmp_organism_parts<-c(tmp_organism_parts, i)
      }
    }
    organism_parts <- tmp_organism_parts
    rm(tmp_organism_parts)
    queries <- c()
    if ("individual" %in% factor_names){
      
      index_to_remove <- which( factor_names == "individual")
      get_factor_names <- factor_names[-index_to_remove]
    } else {
      get_factor_names <- factor_names
    }
    for (i in seq_along(organism_parts)){
      tmp <- ""
      for (j in seq_along(get_factor_names)){
        if (tmp!=""){tmp<-paste0(tmp, ", ")}
        if (get_factor_names[j]=="disease"){
          tmp <- paste0(tmp,dict_Onotology_terms$`Ontology term`[which(dict_Onotology_terms$`App term`==input$Disease)])
          
        }
        
        if (get_factor_names[j]=="disease staging"){
          if (input$Disease=="Control"){tmp <- paste0(tmp,"normal")}
          else {tmp <- paste0(tmp,dict_Onotology_terms$`Ontology term`[which(dict_Onotology_terms$`App term`==input$Stage)])}
        }
        if (get_factor_names[j]=="organism part"){
          tmp <- paste0(tmp,organism_parts[i])
        }
      }
      queries[i]<-tmp
    }
    list(queries, organism_parts)
  })
  
  generateTable <- reactive({
    proteinGroups <- selected_exp_table()
    Gene_row <- selected_query_gene_table()[[1]]
    Condition_col <- selected_query_conditions()[[1]]
    
    Organism_parts <- selected_query_conditions()[[2]]
    values_from_table =list()
    for (i in seq_along(Organism_parts)){
      
      values_from_table[[Organism_parts[i]]] <-   proteinGroups[Gene_row,Condition_col[i]]
      
    }
    values_from_table
  })
  
  
  
  preplots <- reactive({
    map_roi <- c(map)
    organism_parts <- names(generateTable())
    
    organism_parts_eve_template <- c()
    organism_parts_eve_template_index <- c()
    organism_parts_eve_template_index_rep <-c()
    for (i in seq_along(organism_parts)){
      organism_parts_eve_template[i]<- dict_Onotology_terms$EveTemplate[which(dict_Onotology_terms$`Ontology term`==organism_parts[i])]
      grep_tmp <- grep(organism_parts_eve_template[i], labels$text_label)
      organism_parts_eve_template_index<-c(organism_parts_eve_template_index, grep_tmp-1)
      organism_parts_eve_template_index_rep <- c(organism_parts_eve_template_index_rep,  rep(i, length(grep_tmp)))
    }
    
    
    map_roi [!map_roi %in% organism_parts_eve_template_index] <-0
    map_roi = plyr::mapvalues(x=map_roi,
                              from = organism_parts_eve_template_index,
                              to = organism_parts_eve_template_index_rep)
    map_roi = niftiarr(map, map_roi)
    colorConditions <- c()
    for (i in seq_along(generateTable())){
      
      colorConditions[i]<-as.integer(255*generateTable()[[i]]/selected_query_gene_table()[[2]])
      
      if (colorConditions[i]==0){
        colorConditions[i] <- 1
      }
    }
    
    list(map_roi, organism_parts, organism_parts_eve_template, colorConditions)
  })
  
  output$desc <- renderText({
    
    paste("Experiment: ", input$experiment)
    
  })
  
  output$desc2 <- renderText({
    
    paste("Publication: ", links$Publication[which(links$Authorship==input$experiment)])
    
  })
  
  output$desc3 <- renderText({
    paste("Condition: ", input$Disease, if(input$Disease == 'Alzheimer disease'){paste("-",input$Stage)})
    
  })
  
  
  
  
  
  output$table <- renderTable(generateTable())
  
  output$plotBrainRegions <- renderPlot({
    
    
    ortho2(eve_t1,
           y=preplots()[[1]],
           xyz = c(input$x.axis,input$y.axis,input$z.axis),
           mfrow = c(1,4),
           col.y=rainbow(length(preplots()[[2]])),
           legend = preplots()[[3]],
           leg.col=rainbow(length(preplots()[[2]])),
           addlegend = TRUE,
           leg.cex=2,
           leg.x=0,
           leg.y=50)
    
  })
  output$plotBrainExpression <- renderPlot({
    
    
    mycolors <- colorRampPalette(c("blue", "red"))(255)
    ortho2(eve_t1,
           y=preplots()[[1]],
           xyz = c(input$x.axis,input$y.axis,input$z.axis),
           mfrow = c(1,4),
           col.y = mycolors[preplots()[[4]]])
    
    
    text(0.82,.91,"Scale", col = "white")
    text(0.85,.2,"0", col = "white", cex = 0.8, pos = 4)
    text(0.85,.85,selected_query_gene_table()[[2]], col = "white", pos = 4, cex = .8)
    colfunc <- colorRampPalette(c("red", "blue"))
    legend_image <- as.raster(matrix(colfunc(20), ncol=1))
    rasterImage(legend_image, .8, .2, .85,.85)
  })
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)
