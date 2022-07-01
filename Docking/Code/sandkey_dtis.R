library(networkD3)




database_name <- 'NR'
print(database_name)
folder_path <- '../Results/data_per_dataset'
folder_db_path <-  file.path(folder_path, database_name)


file_dtis <- file.path(folder_path, database_name, paste0('dtis_', database_name, '.csv'))

dtis <- read.table(file_dtis, header = TRUE, sep = ";", quote = "")



dtis$value = 1

# Load energy projection data
URL <- "https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json"
Energy <- jsonlite::fromJSON(URL)

 
# Now we have 2 data frames: a 'links' data frame with 3 columns (from, to, value), and a 'nodes' data frame that gives the name of each node.
head( Energy$links )
head( Energy$nodes )
 

p <- sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30)
p

library(htmlwidgets)
saveWidget(p, file=paste0("sankeyEnergy.html"), selfcontained=FALSE)
