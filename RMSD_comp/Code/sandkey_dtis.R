# done in local laptop
library(networkD3)
library(dplyr)



database_name <- 'BindingDB'
print(database_name)



file_path <-  paste0('data_grouped_', database_name, '_w_other.csv')
#file_path <- 'data_grouped_NR_subclass'

links <- read.table(file_path, sep=";", header = TRUE,  quote = "")

# From these flows we need to create a node data frame: 
# it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
  as.character(links$target)) %>% unique()
)


# With networkD3, connection must be provided using id, 
# not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
 

 # Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "count", NodeID = "name", 
              sinksRight=FALSE)
p


# save the widget
library(htmlwidgets)
widget_name <- paste0("sankey", database_name, "_w_other.html")
saveWidget(p, file=widget_name)

