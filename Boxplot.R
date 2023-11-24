################################################################################
#########################Script for automated box plots#########################
################################################################################

###Packages###

{
  #install.packages(tidyverse)
  library(tidyverse)
  #install.packages(tibble)
  library(tibble)
}

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

#Standard boxplot + scatter plot without outliers.
#Each line can be altered or omitted with (#) to best suit your aim

for (col in loop.cols) {

  b <- Q_values[,c(1,which( colnames(Q_values)==col ))]                         #Creates a dataframe (b) with the groups column and the column of the variable being plotted
                                                                                #Note that this is in a for-loop so it will cycle through more than one variable, depending on your input data
    list_quantiles <- tapply(pull(b, col), b$G, quantile)                       #List of quantiles (0%,25%,50%,75%,100%) for each group.

  Q1s <- sapply(1:2, function(i) list_quantiles[[i]][2])                        #List of the 25% quantiles in each group (its the lower limit of the IQR for that group/lower limit of the box in the boxplot)
  Q3s <- sapply(1:2, function(i) list_quantiles[[i]][4])                        #The same list but for the 75% quantiles. Note that in the template data there are only 2 groups.
                                                                                #If you have more than 2 groups just alter sapply(1:2 ...) to sapply(1:x ...) in both Q1s and Q3s, x being your number of groups.
  IQRs <- tapply(pull(b, col), b$G, IQR)                                        #IQR for each group

  Lowers <- Q1s - 1.5*IQRs                                                      #The lower limit of the whiskers in each group
  Uppers <- Q3s + 1.5*IQRs                                                      #The upper limit of the whiskers in each group

  bs <- split(b, b$G)                                                           #Dataframe b split into the various groups

  n <- NULL
  for (i in 1:2){                      #This loop creates a dataframe (n) with outliers removed.
    out <-                             #If you have more than 2 groups just alter for (i in 1:2) to for (i in 1:x), x being your number of groups.
      subset(
        bs[[i]],
        pull(
          bs[[i]],
          col
        ) 
        >= Lowers[i] 
        & 
        pull(
          bs[[i]], 
          col
        ) 
        <= Uppers[i]
      )
    n <- rbind(n, out)                
  }
  
  M <- max(pull(n, col))              #The maximum value for this variable, excluding outliers. Will help us rearrange the y-axis limit later on.
  
  m <- min(pull(n, col))              #The minimum value for this variable, excluding outliers. Will help us rearrange the y-axis limit later on.
  
  p <-
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) + #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each sample)
    geom_boxplot(                                    #Boxplot
      width = .3,                                    #Controls the width of the boxplots
      outlier.shape = NA,                            #Hides (but does not remove) outliers. To plot without outliers see other code chunk bellow this one.
      size = .75,                                    #Controls boxplot line thickness                                   
      coef = 1.5
    ) +
    geom_point(                                      #Scatter plot
      data = n,                                      #To plot only data without outliers (dataframe n)
      aes(
        as.factor(G),
        !!sym(col)
      ),
      position = position_jitter(                    #Scrambles the points in the x-axis to avoid overplotting.
        width = .1,                                  #Controls the width of the scrambling
        height = 0,                                  #Stops scrambling on the y-axis
        seed = 1
      ),
      alpha = .3                                     #Controls transparency of the points
    ) +
    coord_cartesian(                                 #Used to zoom in/out on specific coordinates.Here we use it to set the axis limits.
      ylim = c(m - m/20,                             #Defines the upper limit of the y-axis
               M + M/20                              #Defines the upper limit of the y-axis
      )
    ) +
    scale_y_continuous(labels = scales::scientific) +#Y-axis numbers in scientific format (add # before line to remove)
    scale_x_discrete(labels= c("Control","Disease"))+#X-axis tick labels (group names, in this case)
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
            filter(V_list,Protein == col) %>% select(p_value),       #Gets the p-value from the V_list table. Use round(,x) to round to x decimal places.
            "/",
            "FC =",
            filter(V_list,Protein == col) %>% select(FC)             #Gets the fold-change from the V_list table. Use round(,x) to round to x decimal places.
          )
        ), collapse = " ") ,
      tag = filter(V_list,Protein == col) %>% select(Figure),        #Defines the tag (upper left corner text). Gets it from the figure column of the V_list table
    ) +
    theme(plot.title = element_text(hjust = 0.5),                    #Centers the title
          plot.subtitle = element_text(hjust = 0.5),                 #and the subtitle
          legend.position = "none",                                  #Removes the legend
    )
  
  print(p)
  
}

################################################################################

#Standard boxplot + scatter plot with outliers.
#Each line can be altered or omitted with (#) to best suit your aim

for (col in loop.cols) {
  
  b <- Q_values[,c(1,which( colnames(Q_values)==col ))]                         #Creates a dataframe (b) with the groups column and the column of the variable being plotted
  
  list_quantiles <- tapply(pull(b, col), b$G, quantile)                         #List of quantiles (0%,25%,50%,75%,100%) for each group.
  
  Q1s <- sapply(1:2, function(i) list_quantiles[[i]][2])                        #List of the 25% quantiles in each group (its the lower limit of the IQR for that group/lower limit of the box in the boxplot)
  Q3s <- sapply(1:2, function(i) list_quantiles[[i]][4])                        #The same list but for the 75% quantiles. Note that in the template data there are only 2 groups.
                                                                                #If you have more than 2 groups just alter sapply(1:2 ...) to sapply(1:x ...) in both Q1s and Q3s, x being your number of groups.
  IQRs <- tapply(pull(b, col), b$G, IQR)                                        #IQR for each group
  
  Lowers <- Q1s - 1.5*IQRs                                                      #The lower limit of the whiskers in each group
  Uppers <- Q3s + 1.5*IQRs                                                      #The upper limit of the whiskers in each group
  
  bs <- split(b, b$G)                                                           #Dataframe b split into the various groups
  
  n <- NULL
  for (i in 1:2){                      #This loop creates a dataframe (n) with outliers removed.
    out <-                             #If you have more than 2 groups just alter for (i in 1:2) to for (i in 1:x), x being your number of groups.
      subset(
        bs[[i]],
        pull(
          bs[[i]],
          col
        ) 
        >= Lowers[i] 
        & 
          pull(
            bs[[i]], 
            col
          ) 
        <= Uppers[i]
      )
    n <- rbind(n, out)    
  } 
    o <- NULL
    for (i in 1:2){                      #This loop creates a dataframe (o) with only outliers.
      outo <-                            #If you have more than 2 groups just alter for (i in 1:2) to for (i in 1:x), x being your number of groups.
        subset(
          bs[[i]],
          pull(
            bs[[i]],
            col
          ) 
          < Lowers[i] 
          | 
            pull(
              bs[[i]], 
              col
            ) 
          > Uppers[i]
        )
      o <- rbind(o, outo) 
  }
  
  M <- max(pull(Q_values, col))              #The maximum value for this variable, including outliers. Will help us rearrange the y-axis limit later on.
  
  m <- min(pull(Q_values, col))              #The minimum value for this variable, including outliers. Will help us rearrange the y-axis limit later on.
  
  p <-
    ggplot(Q_values, aes(as.factor(G),!!sym(col))) + #Defines the dataset (Q_values), x values (column G, which stands for group), and y values (quantification values for each variable)
    geom_boxplot(                                    #Boxplot
      width = .3,                                    #Controls the width of the boxplots
      outlier.shape = NA,                            #Hides outliers. We will plot them latter with "jitter" to avoid overplotting.
      size = .75,                                    #Controls boxplot line thickness                                   
      coef = 1.5
    ) +
    geom_point(                                      #Scatter plot
      data = n,                                      #To plot only data without outliers (dataframe n)
      aes(
        as.factor(G),
        !!sym(col)
      ),
      position = position_jitter(                    #Scrambles the points in the x-axis to avoid overplotting.
        width = .1,                                  #Controls the width of the scrambling
        height = 0,                                  #Stops scrambling on the y-axis
        seed = 1
      ),
      alpha = .3                                     #Controls transparency of the points
    ) +
    geom_point(                                      #Scatter plot
      data = o,                                      #To plot only outliers (dataframe o)
      aes(
        as.factor(G),
        !!sym(col)
      ),
      position = position_jitter(                    #Scrambles the points in the x-axis to avoid overplotting.
        width = .1,                                  #Controls the width of the scrambling
        height = 0,                                  #Stops scrambling on the y-axis
        seed = 1
      ),
    ) +
    coord_cartesian(                                 #Used to zoom in/out on specific coordinates.Here we use it to set the axis limits.
      ylim = c(m - m/20,                             #Defines the lower limit of the y-axis
               M + M/20                              #Defines the upper limit of the y-axis
      )
    ) +
    scale_y_continuous(labels = scales::scientific) +#Y-axis numbers in scientific format (add # before line to remove)
    scale_x_discrete(labels= c("Control","Disease"))+#X-axis tick labels (group names, in this case)
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
          legend.position = "none",                                  #Removes the legend
    )
  
  print(p)
  
}