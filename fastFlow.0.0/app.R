#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(DT)
library(shiny)
library(dplyr) 
library(readxl)
library(tidyr)
library(ggplot2)
library(plotly)
library(rhandsontable)
library(stringr)
library(tidyverse)
# eventually we will put the custom functions in their own file and use ">source()" to reference them, but for now I'm just putting them in-line
cleanFlow_0.2 = function(mytable_in,control_ID = "TIA_rep1"){
    
    #get just the FL2-H data from plot 3 rows using an index of rows with the string "Plot 3"
    plot3_index = which(str_detect(mytable_in$...1,"Plot 3") == T)
    plot3_index = c(plot3_index,plot3_index+1)
    plot3_index = sort(plot3_index)
    mytable = mytable_in[plot3_index,]
    #copy JUST the sample metadata into the same row as the analysis data
    mytable[which(as.integer(rownames(mytable)) %% 2 == 0),1] = mytable[which(as.integer(rownames(mytable)) %% 2 == 0)-1,1]
    #subset out just the analysis data, which now has the necessary metadata appended
    clean_cytoDat = mytable[which(as.integer(rownames(mytable)) %% 2 == 0),]
    #append the column names back and correct the first column name
    colnames(clean_cytoDat) = mytable[1,]
    colnames(clean_cytoDat)[1] = "metadata"
    colnames(clean_cytoDat)[3] = "events_per_ul"
    #clean up memory
    rm(mytable)
    
    #parse the metadata
    #TODO the function works well for standard names that carp samples would have e.g., 18-057-05 
    #         but not great for how a random sample from another species might be e.g., LST 23 (lake sturgeon 23)
    #         do we generalize in the function, adjust the naming scheme for non-carp, or make the app only for carp?
    #         right now the non-standardized names still plot out but have "NA-rep_n" as their ID.
    #isolate the metadata
    metaDat_parse = clean_cytoDat[,1]
    #split apart all the elements
    metaDat_parse = metaDat_parse %>% separate(metadata, c("NA","plot_ID","well_ID","fiscal_year","case_ID","individual_ID"))
    
    #reconstruct the elements of the case ID
    metaDat_parse$case_ID = paste(metaDat_parse$fiscal_year,metaDat_parse$case_ID, sep = "-")
    metaDat_parse = metaDat_parse %>% select("case_ID","individual_ID","well_ID","plot_ID")
    #loop through the individual_ID values and append replicate number to present replicates
    metaDat_parse$replicate = NA
    metaDat_parse$individual_ID = as.integer(metaDat_parse$individual_ID)
    for (i in 1:nrow(metaDat_parse)) {
        metaDat_parse$replicate[i] = length(which(as.vector(metaDat_parse$individual_ID[1:i]) %in% metaDat_parse$individual_ID[i])) 
    }
    #re-order the columns to more logical orientation
    metaDat_parse = metaDat_parse %>% select("case_ID","individual_ID","replicate","well_ID","plot_ID")
    
    #append the formatted metadata to the analysis data
    clean_cytoDat = cbind.data.frame(metaDat_parse,clean_cytoDat)
    #drop the un-formatted metadata column, as this has now been parsed
    clean_cytoDat = clean_cytoDat %>% select(-"metadata")
    #sort by individual ID, as this is more logical than well
    clean_cytoDat = clean_cytoDat %>% arrange(individual_ID,Count)
    #format the control data individual ID
    clean_cytoDat$individual_ID[which(clean_cytoDat$case_ID == "TIA-NA")] = "TIA"
    #format CV and mean FL2-H to doubles for calculations
    clean_cytoDat$`Mean FL2-H` = as.double(clean_cytoDat$`Mean FL2-H`)
    clean_cytoDat$`CV FL2-H` = as.double(clean_cytoDat$`CV FL2-H`)
    
    #construct a normal distribution for each sample based on mean and CV of FL2-H
    #first calculate standard deviation based on the mean and CV
    clean_cytoDat$standardDev = (clean_cytoDat$`CV FL2-H`/100)*clean_cytoDat$`Mean FL2-H`
    #build a random normal distribution based on the mean and calculated standard deviation for each well, stored as list
    for (i in 1:nrow(clean_cytoDat)) {
        clean_cytoDat$distNorm[i] = I(list(rnorm(n = 10000,mean = clean_cytoDat$`Mean FL2-H`[i],sd = clean_cytoDat$standardDev[i])))
    }
    
    #insert any useful metadata that we want to show with the plot
    clean_cytoDat$sample_rep = paste(clean_cytoDat$individual_ID,"_rep",clean_cytoDat$replicate, sep = "")
    
    #we get the estimated DNA mass in here, this is  a calculation based on the TIA control we select as our standard
    clean_cytoDat$DNA_mass_pg = round((2.4/clean_cytoDat[which(clean_cytoDat$sample_rep == control_ID),"Mean FL2-H"])*clean_cytoDat[,"Mean FL2-H"],3)
    
    #convert the random normal distribution to a kernel density estimate for each rep for plotting
    density1 <- density(clean_cytoDat$distNorm[[1]])
    df_density_i = data.frame(x=density1$x, y=density1$y)
    #append the replicate metadata for splitting in the plot
    df_density_i$sample_rep = paste(clean_cytoDat$individual_ID[1],clean_cytoDat$replicate[1],sep = "_rep")
    
    for (i in 2:nrow(clean_cytoDat)) {
        density_i = density(clean_cytoDat$distNorm[[i]])
        density_i = data.frame(x=density_i$x, y=density_i$y)
        density_i$sample_rep = paste(clean_cytoDat$individual_ID[i],clean_cytoDat$replicate[i],sep = "_rep")
        df_density_i = rbind.data.frame(df_density_i,density_i)
    }
    
    #TODO: can i code the sample rep info as a double by concatenating with a "."??? then I can sort the legend by this... or keep as is for readability and sort based on other value will work fine too...
    
    
    
    if (is.na(control_ID) == F) {
        #formatting the output table to have all the valuable control data associated with the samples it was associated with
        clean_cytoDat$control_ID = control_ID
        clean_cytoDat$`control Median FL2-H` = clean_cytoDat$`Median FL2-H`[which(clean_cytoDat$sample_rep == control_ID)]
        clean_cytoDat$`control Mean FL2-H` = clean_cytoDat$`Mean FL2-H`[which(clean_cytoDat$sample_rep == control_ID)]
        clean_cytoDat$`control CV FL2-H` = clean_cytoDat$`CV FL2-H`[which(clean_cytoDat$sample_rep == control_ID)]
        clean_cytoDat$control_Count = clean_cytoDat$Count[which(clean_cytoDat$sample_rep == control_ID)]
        clean_cytoDat$control_events_per_ul = clean_cytoDat$events_per_ul[which(clean_cytoDat$sample_rep == control_ID)]
        mta_dat4Plot = clean_cytoDat[,c("case_ID","individual_ID","sample_rep","replicate","DNA_mass_pg","Median FL2-H","Mean FL2-H","CV FL2-H","Count","events_per_ul","control_ID","control Median FL2-H","control Mean FL2-H","control CV FL2-H","control_Count","control_events_per_ul")]
    }
    if (is.na(control_ID)) {
        mta_dat4Plot = clean_cytoDat[,c("case_ID","individual_ID","sample_rep","replicate","DNA_mass_pg","Median FL2-H","Mean FL2-H","CV FL2-H","Count","events_per_ul")]
    }
    #TODO: need to add ploidy, analysis date, and lab notes columns
    # [,c("case_ID","individual_ID","sample_rep","replicate","DNA_mass_pg",PLOIDY,"Median FL2-H","Mean FL2-H","CV FL2-H","Count","events_per_ul","control_ID","control Median FL2-H","control Mean FL2-H","control CV FL2-H","control_Count","control_events_per_ul",ANALYSIS DATE, LAB NOTES]
    
    #append the count and event per UL data based on sample_rep ID
    df_density_i$count = mta_dat4Plot$Count[match(df_density_i$sample_rep, mta_dat4Plot$sample_rep)]
    df_density_i$events_per_ul = mta_dat4Plot$events_per_ul[match(df_density_i$sample_rep, mta_dat4Plot$sample_rep)]
    df_density_i$DNA_mass_pg = mta_dat4Plot$DNA_mass_pg[match(df_density_i$sample_rep, mta_dat4Plot$sample_rep)]
    df_density_i$individual_ID = as.integer(mta_dat4Plot$individual_ID[match(df_density_i$sample_rep, mta_dat4Plot$sample_rep)])
    df_density_i$replicate = mta_dat4Plot$replicate[match(df_density_i$sample_rep, mta_dat4Plot$sample_rep)]
    df_density_i$sample_rep = factor(df_density_i$sample_rep, levels = unique(df_density_i[order(df_density_i$individual_ID,df_density_i$replicate),"sample_rep"]), ordered = T)
    #TODO need to reorder the levels so that they plot in numeric sequential despite being factors
    #   i.e. 2_rep1 should be before 10_rep1
    
    #generate hexidecimal colors for each sample ID
    df_density_i$rgb = "#000000"
    colorSwatch = NULL
    df_density_i[which(is.na(df_density_i$individual_ID)),"individual_ID"] = "TIA"
    # i = "11_rep1"
    i = "TIA"
    # i = 1
    #set all the black color for the control samples, identified with TIA beacuse we only use tilapia
    for (i in unique(df_density_i$individual_ID)) {
        index_matchI_rows = which(df_density_i$individual_ID == i)
        if ("TRUE" %in% grepl("TIA",unique(df_density_i[index_matchI_rows,"sample_rep"]))) {
            #check if there are multiple replicates for the sample 
            n_replicates = n_distinct(df_density_i[index_matchI_rows,"sample_rep"])
            colorSwatch = c(colorSwatch,rep("#000000",n_replicates))
            df_density_i$rgb[index_matchI_rows] = rgb(sample(seq(0,1,1),1), sample(seq(0,1,1),1), sample(seq(0,1,1),1), maxColorValue=255)
            next
        }
        
        else{
            #check if there are multiple replicates for the sample 
            n_replicates = n_distinct(df_density_i[index_matchI_rows,"sample_rep"])
            
            #assign the color to the swatch n times, where n = the number of replicates present
            colorValue = rgb(sample(seq(50,255,1),1), sample(seq(50,255,1),1), sample(seq(50,255,1),1), maxColorValue=255)
            colorSwatch = c(colorSwatch,rep(colorValue,n_replicates))
            df_density_i$rgb[index_matchI_rows] = colorValue
            next
        }
    }
    
    #plot the kernel density estimates by replicate with random colors
    x_ax = list(title = "FL2-H")
    y_ax = list(title = "frequency")
    
    fig <- plot_ly() %>% 
        add_trace(x = ~df_density_i$x, 
                  y = ~df_density_i$y, 
                  split = ~df_density_i$sample_rep,
                  name = ~df_density_i$sample_rep,
                  color = ~df_density_i$sample_rep,
                  colors = colorSwatch,
                  hoverinfo = 'text',
                  text = ~paste(
                      '</br> Sample_replicate: ', df_density_i$sample_rep,
                      '</br> DNA mass (pg): ', df_density_i$DNA_mass_pg,
                      '</br> Event Count: ', df_density_i$count,
                      '</br> events/uL: ', df_density_i$events_per_ul
                  )
        )%>% 
        layout(xaxis = x_ax,
               yaxis = y_ax)
    
    #return a list holding the clean data table, and the plotted input data
    cleanFlow_out = list(mta_dat4Plot, fig)
    return(cleanFlow_out)
}




# Define UI for application that draws a histogram
ui <- fluidPage(
    ############    
    sidebarPanel(
        #Allow user to read in .xlsx file with data to process
        fileInput('file1', label = 'Select file', accept = c(".xlsx")),
        # #Provide a button to choose the control sample for DNA mass calculation
        selectInput('control', label = 'Select control', "Select file above"),
        # download processed data
        downloadButton("downloadData", "Download")
    ),
    
    mainPanel(
        fluidRow(
            column(width = 12, offset = 0, plotlyOutput(outputId = "ploidyPlot"))
        ),
        fluidRow(
            column(width = 12, rHandsontableOutput("database_out"))
        ),
    )
)


server <- function(input, output, session) {
    ############
    
    #This function is responsible for producing the database table and plot, they are stored in a list in that specific order
    cleaned_flow <- reactive({infile = input$file1
    if(is.null(infile)) {
        # User has not uploaded a file yet
        return(NULL)
    }
    if(input$control == "Select file above"){
        fileData = read_excel(infile$datapath, col_names = F)
        cleanFlow_0.2(fileData)    
    }
    else{
        fileData = read_excel(infile$datapath, col_names = F)
        cleanFlow_0.2(fileData,input$control)
    }
    
    })
    # once an input file is present, this allows user to select the sample ID for the control, 
    #which will then trigger calculation of the DNA mass estimate for all samples
    #this assumes that the control is talapia, which is the standard workflow at the La Crosse FHC
    observe({
        updateSelectInput(session, "control",
                          label = "control",
                          choices = cleaned_flow()[[1]][3],
                          selected = input$control
        )
    })
    
    #database table is captured as output
    output$database_out = renderRHandsontable({
        rhandsontable(cleaned_flow()[[1]]) 
    })
    
    # capture density plot output
    output$ploidyPlot = renderPlotly({cleaned_flow()[[2]]})
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("results_",gsub(".xlsx","",input$file1), ".csv", sep = "")
        },
        content = function(file) {
            write.csv(cleaned_flow()[[1]], file, row.names = FALSE)
            # (database_out, file, row.names = FALSE)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)