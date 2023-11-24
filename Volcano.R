################################################################################
############################Script for Volcano plots############################
################################################################################

###Packages###

{
  #install.packages(tidyverse)
  library(tidyverse)
  #install.packages(tibble)
  library(tibble)
}

#Info on gghalves utility:
#https://cran.r-project.org/web/packages/gghalves/vignettes/gghalves.html

###Data###

#Load the data (according to the templates) in a .txt file (tab delimited)

#This table contains other important information about each variable
#This can be used to annotate the final image with labels, names, fold-change, p-value, and others
#The variable IDs in the first column must match the ones on each column of the Q_values table

V_list <- readr::read_delim("V_list_R.txt", delim = "\t")

#Run the following commands to see the top rows of the table

tibble::tibble(V_list)

###Theme###

#The theme_set function here is used to standardize several non-data aspects of all plots in this script
#In this case theme_classic is used with a base letter size of 16
#You can check out other available themes and their respective arguments in the following link
#https://ggplot2.tidyverse.org/reference/ggtheme.html

theme_set(theme_classic(base_size = 16))

################################################################################

###Volcano plot###

{
a <- -log10(0.05)                                                               #Reference values used to draw cut-off lines in the plot
b <- log2(1.5)                                                                  #
c <- log2(0.66)                                                                 #

lp <- -log10(V_list$p_value)                                                    #Log transformation of p-values
lfc <- log2(V_list$FC)                                                          #and FC values

l1 <- data.frame(                                                               #Stores the log transformed values in a dataframe (l1)
  lp,
  lfc
)

l2 <-                                                                           #A dataframe (l2) only for up and down regulated variables
  filter(l1,(lp >= a & lfc >= b)|(lp >= a & lfc <= c)) %>%
  mutate(note = if_else(lfc >= b, "Up", "Down"))

l3 <-                                                                           #A dataframe (l3) for the remaining variables
  filter(l1,!((lp >= a & lfc >= b)|(lp >= a & lfc <= c)))

ggplot(l2, aes(lfc, lp, fill = note)) +             #Defines the dataset (l2), x values (log transformed FC values), and y values (log transformed p-values)
  geom_hline(yintercept = a) +                      #Draws a line at p-value = 0.05                            
  geom_vline(xintercept = b) +                      #Draws a line at FC = 1.5  
  geom_vline(xintercept = c) +                      #Draws a line at FC = 0.66
  geom_point(                                       #Scatter plot for up/down regulated variables
    shape = 21,                                     #Defines the shape of the points
    color = "black",                                #Defines the point outline color
    alpha = .75) +                                  #Defines the point transparency
  scale_fill_manual(
    values = c(
      "blue",                                       #Fill color for points where FC below/equal to 0.66
      "firebrick"                                   #Fill color for points where FC above/equal to 1.5
    ),
    labels =
      c("FC \U2264 0.66",                           #Legend text for points where FC below/equal to 0.66                   
        "FC \U2265 1.5"                             #Legend text for points where FC above/equal to 1.5
      )                                             #The "\UNICODE" standard for the (more than / less than) or (equal to symbols)
  ) +
  geom_point(                                       #Scatter plot for remaining variables
    data = l3,                                      #Overrides the dataframe being used (in this case to be l3 containg only the remaining variables)
    shape = 21,                                     #Defines the shape of the points
    color = "black",                                #Defines the point outline color
    fill = "black",                                 #Defines the point fill color
    alpha = .3                                      #Defines the point transparency
  ) +
  labs(
    x = "log2 FC",                                  #Defines the x-axis label
    y = "-log10 p-value",                           #Defines the y-axis label
    title = "Altered Proteins",                     #Defines the title
    subtitle = "Control vs Disease",                #Defines the subtitle
    tag = "A",                                      #Defines the tag (upper left corner text).
    fill = NULL                                     #Overwrites the legend title (in this case removes it)
  ) +
  theme(plot.title = element_text(hjust = 0.7),     #Centers the title
        plot.subtitle = element_text(hjust = 0.7)   #and the subtitle
  )
}