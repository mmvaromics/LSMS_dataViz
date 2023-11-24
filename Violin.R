################################################################################
#######################Script for automated violin plots########################
################################################################################

###Packages###

{
#install.packages(tidyverse)
library(tidyverse)
#install.packages(tibble)
library(tibble)
#install.packages(ggforce)
library(ggforce)
#install.packages(ggdist)
library(ggdist)
#install.packages(gghalves)
library(gghalves)
}

#Info on gghalves utility:
#https://cran.r-project.org/web/packages/gghalves/vignettes/gghalves.html

###Data###

#Load the data (according to the templates) in a .txt file (tab delimited)

#This table contains quantification data for each variable and sample grouping information
#Column G is group
#Column subG is for sub-group
#Column S is for sample ID
#The remaining coumns are for the quantification values of your variables of interest
#Any number of columns may be added with other categorical variables

Q_values <- readr::read_delim("Q_values_R.txt", delim = "\t")

#This table contains other important information about each variable
#This can be used to annotate the final image with labels, names, fold-change, p-value, and others
#The variable IDs in the first column must match the ones on each column of the Q_values table

V_list <- readr::read_delim("V_list_R.txt", delim = "\t")

#Run the following commands to see the top rows of each table

tibble::tibble(Q_values)
tibble::tibble(V_list)

###Loop###

#Use this line to define for which variables the respective values should be plotted
#In the case of the template data there are 6 variables
#Their respective values are in columns 4 through 9 of the Q_Values table
#This explains the numbers in the following command and should be suited to your input data

loop.cols <- names(Q_values)[4:9]

###Theme###

#The theme_set function here is used to standardize several non-data aspects of all plots in this script
#In this case theme_classic is used with a base letter size of 16
#You can check out other available themes and their respective arguments in the following link
#https://ggplot2.tidyverse.org/reference/ggtheme.html

theme_set(theme_classic(base_size = 16))

################################################################################

###Violin plot###

#Standard violin with median and IQR. Colored groups and log scale option.
#Each line can be altered or omitted with (#) to best suit your aim


for (col in loop.cols) {
  p <-
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) + #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each variable)
    geom_violin(                                     #Violin plot
      aes(color = G),                                #Uses the dataset defined in ggplot() and separates the data by group (column G) with different colors
      fill = "NA",                                   #Removes the fill
      adjust = 1.5,                                  #Adjusts the "smoothness" of the plots
      linewidth = 1                                  #Controls the line thickness of the plots
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },    #Plots a line from the 0.25 quantile
      fun.max = function(z) { quantile(z,0.75) },    #to the 0.75 quantile (also known as IQR)
      fun = median,                                  #Draws a point at the median
      size = 0.8,                                    #Controls the size of the median point
      color = c("firebrick", "blue"),                #Changes the color of the plotted geoms
      linewidth = 0.7                                #Controls the thickness of the IQR line
    ) +
    scale_y_log10() +                                #Log10 scales the y-axis (omit with # to go back to the normal scale)
    scale_x_discrete(labels= c("Control","Disease"))+#X-axis tick labels (group names, in this case)
    scale_color_manual(                              #Overwrites the violin plot color
      values = c("firebrick", "blue") 
    ) +
    labs(
      x = "Group",                                   #Defines the x-axis label
      y = "Relative Abundance",                      #Defines the y-axis label
      title = "Control vs Disease",                  #Defines the title
      subtitle = paste(                              #Defines the subtitle
        as.character(
          c(
            filter(V_list,Protein == col) %>% select(Protein_names), #Gets the protein name from the V_list table
            " ",
            "(",
            col,                                                     #Gets the Accession Number from the Q_values table (it is also the column name for the respective variable)
            ")"
          )
        )
        , collapse = ""),
      caption = paste(                               #Defines the caption (lower right corner text)
        as.character(
          c(
            "p =",
            filter(V_list,Protein == col) %>% select(p_value),       #Gets the p-value from the V_list table.  Use round(,x) to round to x decimal places.
            "/",
            "FC =",
            filter(V_list,Protein == col) %>% select(FC)             #Gets the fold-change from the V_list table.  Use round(,x) to round to x decimal places.
          )
        ), collapse = " ") ,
      tag = filter(V_list,Protein == col) %>% select(Figure),        #Defines the tag (upper left corner text). Gets it from the figure column of the V_list table
    ) +
    theme(plot.title = element_text(hjust = 0.5),                    #Centers the title
          plot.subtitle = element_text(hjust = 0.5),                 #and the subtitle
          legend.position = "none"                                   #Removes the legend
    )
  
  print(p)
  
}

###Violin + Scatter plots###

#Standard violin + scatter plot, with median and IQR. Colored groups and log scale option.
#Each line can be altered or omitted with (#) to best suit your aim

for (col in loop.cols) {
  p <-
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) + #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each variable)
    geom_point(                                      #Scatter plot
      data = dplyr::filter(Q_values, G=="D"),        #Filters the data to include only samples from group D
      shape = 21,                                    #The shape of the plotted points
      color = "transparent",                         #Defines the outline color of each point
      fill = "black",                                #Defines the fill color of each point
      position = position_jitter(                    #Scrambles the points in the x-axis to avoid overplotting.
        width = .15,                                 #Controls the width of the scrambling
        height = 0,                                  #Stops scrambling on the y-axis
        seed = 2021),                                
      size = 2,                                      #Defines the point size
      alpha = .3                                     #Defines the point transparency
    ) +
    geom_point(                                      #Another scatter plot, plotted over the last one, to highlight the point outline
      data = dplyr::filter(Q_values, G=="D"),
      shape = 21,
      color = "black",                               #Notice that for this one the fill is transparent and not the color
      fill = "transparent",
      stroke = 1,                                    #Defines the thickness of the point outline
      position = position_jitter(
        width = .15, 
        height = 0,
        seed = 2021),
      size = 2,
    ) +
    geom_point(                                      #Scatter plot but now for the CT group
      data = dplyr::filter(Q_values, G=="CT"),
      shape = 21,
      color = "transparent",
      fill = "black",
      position = position_jitter(
        width = .15, 
        height = 0,
        seed = 27),
      size = 2,
      alpha = .3
    ) +
    geom_point(                                      #Another scatter plot, plotted over the last one, to highlight the point outline
      data = dplyr::filter(Q_values, G=="CT"),
      shape = 21,
      color = "black",
      fill = "transparent",
      stroke = 1,
      position = position_jitter(
        width = .15, 
        height = 0,
        seed = 27),
      size = 2,
    ) +
    geom_violin(                                     #Violin plot
      aes(color = G),                                #Uses the dataset defined in ggplot() and separates the data by group (column G) with different colors
      fill = "NA",                                   #Removes the fill
      adjust = 1.5,                                  #Adjusts the "smoothness" of the plots
      linewidth = 1                                  #Controls the line thickness of the plots
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },    #Plots a line from the 0.25 quantile
      fun.max = function(z) { quantile(z,0.75) },    #to the 0.75 quantile (also known as IQR)
      fun = median,                                  #Draws a point at the median
      size = 0.8,                                    #Controls the size of the median point
      color = c("firebrick", "blue"),                #Changes the color of the plotted geoms
      linewidth = 0.7                                #Controls the thickness of the IQR line
    ) +
    scale_y_log10() +                                #Log10 scales the y-axis (omit with # to go back to the normal scale)
    scale_x_discrete(labels= c("Control","Disease"))+#X-axis tick labels (group names, in this case)
    scale_color_manual(                              #Overwrites the violin plot color
      values = c("firebrick", "blue") 
    ) +
    labs(
      x = "Group",                                   #Defines the x-axis label
      y = "Relative Abundance",                      #Defines the y-axis label
      title = "Control vs Disease",                  #Defines the title
      subtitle = paste(                              #Defines the subtitle
        as.character(
          c(
            filter(V_list,Protein == col) %>% select(Protein_names), #Gets the protein name from the V_list table
            " ",
            "(",
            col,                                                     #Gets the Accession Number from the Q_values table (it is also the column name for the respective variable)
            ")"
          )
        )
        , collapse = ""),
      caption = paste(                               #Defines the caption (lower right corner text)
        as.character(
          c(
            "p =",
            filter(V_list,Protein == col) %>% select(p_value),       #Gets the p-value from the V_list table.  Use round(,x) to round to x decimal places.
            "/",
            "FC =",
            filter(V_list,Protein == col) %>% select(FC)             #Gets the fold-change from the V_list table.  Use round(,x) to round to x decimal places.
          )
        ), collapse = " ") ,
      tag = filter(V_list,Protein == col) %>% select(Figure),        #Defines the tag (upper left corner text). Gets it from the figure column of the V_list table
    ) +
    theme(plot.title = element_text(hjust = 0.5),                    #Centers the title
          plot.subtitle = element_text(hjust = 0.5),                 #and the subtitle
          legend.position = "none"                                   #Removes the legend
    )
  
  print(p)
  
}

#Standard violin + scatter plot, with median and IQR. Colored sub-groups and log scale option.
#Each line can be altered or omitted with (#) to best suit your aim

for (col in loop.cols) {
  
  p <- 
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) +   #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each variable)
    geom_violin(                                       #Violin plot
      fill = "NA",                                     #Removes the fill
      adjust = 1.5,                                    #Adjusts the "smoothness" of the plots
      linewidth = 1                                    #Controls the line thickness of the plots
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },      #Plots a line from the 0.25 quantile
      fun.max = function(z) { quantile(z,0.75) },      #to the 0.75 quantile (also known as IQR)
      fun = median,                                    #Draws a point at the median
      size = 0.8,                                      #Controls the size of the median point
      linewidth = 0.7                                  #Controls the thickness of the IQR line
    ) +
    scale_y_log10() +                                  #Log10 scales the y-axis (omit with # to go back to the normal scale)
    scale_x_discrete(labels= c("Control","Disease"))+  #X-axis tick labels (group names, in this case)
    gghalves::geom_half_point(                         #Scatter plot for group D
      data = dplyr::filter(Q_values, G=="D"),
      aes(fill = subG),                                #Sorts the points by sub-group (column subG) with different fill colors
      side = "l", 
      range_scale = 0.5,                               #Controls the dispersion of points, in each sub-group, on the x-axis
      transformation = position_jitter(                #Scrambles the points in the x-axis to avoid overplotting.
        height = 0,                                    #Stops scrambling on the y-axis
        seed = 1
      ),
      shape = 21,                                      #The shape of the plotted points
      color = "black",                                 #Defines the outline color of each point
      size = 2,                                        #Defines the point size
      stroke = 1.2                                     #Defines the thickness of the point outline
    ) +
    geom_point(                                        #Scatter plot for group CT
      data = dplyr::filter(Q_values, G=="CT"),
      aes(fill = subG),
      shape = 21,
      color = "black",
      stroke = 1.2,
      position = position_jitter(
        width = .1, 
        height = 0,
        seed = 27), 
      size = 2
    ) +
    labs(
      x = "Group",                                     #Defines the x-axis label
      y = "Relative Abundance",                        #Defines the y-axis label
      title = "Control vs Disease",                    #Defines the title
      subtitle = paste(                                #Defines the subtitle
        as.character(
          c(
            filter(V_list,Protein == col) %>% select(Protein_names), #Gets the protein name from the V_list table
            " ",
            "(",
            col,                                                     #Gets the Accession Number from the Q_values table (it is also the column name for the respective variable)
            ")"
          )
        )
        , collapse = ""),
      caption = paste(                                 #Defines the caption (lower right corner text)
        as.character(
          c(
            "p =",
            filter(V_list,Protein == col) %>% select(p_value),       #Gets the p-value from the V_list table.  Use round(,x) to round to x decimal places.
            "/",
            "FC =",
            filter(V_list,Protein == col) %>% select(FC)             #Gets the fold-change from the V_list table.  Use round(,x) to round to x decimal places.
          )
        ),
        collapse = " "),
      tag = filter(V_list,Protein == col) %>% select(Figure),        #Defines the tag (upper left corner text). Gets it from the figure column of the V_list table
      fill = "Sub-Groups"                                            #Overwrites the scatter plot legend title
    ) +
    theme(plot.title = element_text(hjust = 0.5),                    #Centers the title
          plot.subtitle = element_text(hjust = 0.5)                  #and the subtitle
    )
  
  print(p)
}

################################################################################

###Raincloud plot###

#based on: https://z3tt.github.io/Rainclouds/

set.seed(27)                 #Run this before ensure points are persistently plotted between runs.
for (col in loop.cols) {                                   #If you change the seed number also change it at the end of the for-loop.
  p <-
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) +       #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each variable)
    stat_halfeye(                                          #Draws a half violin plot
      width = .3,                                          #Controls the width of the plot
      .width = 0,                                          #Draws a median point, IQR line and full range of the data line (keep at 0 to hide full range, median and IQR can be added afterwards)
      adjust = 1                                           #Controls the "smothness" of the plot
    ) + 
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },          #Plots a line from the 0.25 quantile
      fun.max = function(z) { quantile(z,0.75) },          #to the 0.75 quantile (also known as IQR)
      fun = median                                         #Draws a point at the median
    ) +
    geom_half_point(                                       #Scatter plot
      side = "l",                                          #On the left side
      transformation = position_jitter(                    #Scrambles the points in the x-axis to avoid overplotting.
        height = 0,                                        #Stops scrambling on the y-axis
        seed = 1
      ),
      range_scale = 0.7,                                   #Point dispersion
      shape = 21,                                          #The shape of the plotted points (See Help folder for other shapes)
      color = "black",                                     #Defines the outline color of each point
      stroke = .8,                                         #Defines the thickness of the point outline
      size = 1.5                                           #Defines the point size
    ) +
    labs(
      x = "Group",                                   #Defines the x-axis label
      y = "Relative Abundance",                      #Defines the y-axis label
      title = "Control vs Disease",                  #Defines the title
      subtitle = paste(                              #Defines the subtitle
        as.character(
          c(
            filter(V_list,Protein == col) %>% select(Protein_names), #Gets the protein name from the V_list table
            " ",
            "(",
            col,                                                     #Gets the Accession Number from the Q_values table (it is also the column name for the respective variable)
            ")"
          )
        )
        , collapse = ""),
      caption = paste(                                               #Defines the caption (lower right corner text)
        as.character(
          c(
            "p =",
            filter(V_list,Protein == col) %>% select(p_value),       #Gets the p-value from the V_list table.  Use round(,x) to round to x decimal places.
            "/",
            "FC =",
            filter(V_list,Protein == col) %>% select(FC)             #Gets the fold-change from the V_list table.  Use round(,x) to round to x decimal places.
          )
        ), collapse = " ") ,
      tag = filter(V_list,Protein == col) %>% select(Figure),        #Defines the tag (upper left corner text). Gets it from the figure column of the V_list table
    ) +
    theme(plot.title = element_text(hjust = 0.5),                    #Centers the title
          plot.subtitle = element_text(hjust = 0.5),                 #and the subtitle
          legend.position = "none"                                   #Removes the legend
    ) +
    scale_y_log10() +                                                #Log10 scales the y-axis (omit with # to go back to the normal scale)
    scale_x_discrete(labels= c("Control","Disease"))+                #X-axis tick labels (group names, in this case)
    coord_flip()                                                     #Flips the coordinates of plotted elements (here used to flip the axis so that rainclouds are horizontal)
  
  print(p)
  
  set.seed(27)
}