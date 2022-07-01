


library("rjson")
library('ggplot2')

data <- fromJSON(file="test_all_datasets_info.json")

names(data)

list_dbs <- c("NR", "IC")


pdf(file = "testR.pdf")

library(devtools)
library(circlize)
library(RColorBrewer)


list_names <- names(data)

colors <- brewer.pal(length((list_names)), "Set2")

list_degrees <- c()

pdf(file = "testR.pdf")
xl = c(1,1e1)
for (i in 1:length(list_names))
{
      print(list_names[i]);
     p <- data[[list_names[i]]][['list_number_degrees']] #[['list_subgraphs_size']]
      list_degrees <- append(list_degrees, p)
     plot(tabulate(p), type='l', xlim=xl, ylim=c(1,7000),col=colors[i], 
            ylab="Count", xlab="Degree", log='x', lwd = 2)
      par(new=TRUE)
}
legend("topright", legend=list_names,
       col=colors, lty=1, cex=0.8)
dev.off()

data.frame(data[[list_names[i]]][['list_number_degrees']])



df = data.frame()

# Defining a for loop with 30 iterations
df$new <- c(3, 3, 6, 7, 8, 12)

df_new <- cbind(df, new)





# Creating variable
x1= tabulate(data[[list_names[1]]][['list_number_degrees']])
x2= tabulate(data[[list_names[2]]][['list_number_degrees']])
x3= tabulate(data[[list_names[3]]][['list_number_degrees']])
x4= tabulate(data[[list_names[4]]][['list_number_degrees']])
x5= tabulate(data[[list_names[5]]][['list_number_degrees']])
x6= tabulate(data[[list_names[6]]][['list_number_degrees']])
x7= tabulate(data[[list_names[7]]][['list_number_degrees']])
x8= tabulate(data[[list_names[8]]][['list_number_degrees']])


# Finding maximum length
max_ln <- max(c(length(x1), length(x2), length(x3),length(x4), 
                  length(x5), length(x6), length(x7), length(x8)))

gfg_data<- data.frame(index=1:max_ln,
                      col1 = c(x1,rep(NA, max_ln - length(x1))),
                      col2 = c(x2,rep(NA, max_ln - length(x2))),
                      col2 = c(x3,rep(NA, max_ln - length(x3))),
                      col2 = c(x4,rep(NA, max_ln - length(x4))),
                      col2 = c(x5,rep(NA, max_ln - length(x5))),
                      col2 = c(x6,rep(NA, max_ln - length(x6))),
                      col2 = c(x7,rep(NA, max_ln - length(x7))),
                      col2 = c(x8,rep(NA, max_ln - length(x8))))


columns <- c('index')
columns <- append(columns, list_names)
names(gfg_data) <- columns

#melt data frame into long format
df <- melt(gfg_data ,  id.vars = 'index', variable.name = 'datasets')

p2 <- ggplot(df, aes(index, value)) +
    geom_line(aes(colour = datasets))+
    geom_point(alpha = 0.3) +
    geom_smooth(se = FALSE) +
    scale_x_continuous(limits = c(0, 10)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw()



pdf(file = "testR_3.pdf", height=20, width=20)
#create line plot for each column in data frame
ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets))+
            scale_x_continuous(trans='log10')
dev.off()


sub.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets))+
            xlim (0, 10) +
            ylim (0, 1000) 


inset.plot <- main.plot + theme(legend.position = "none")

library(cowplot)

main.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets), size=0.8)+
            #geom_point(size=0.1)+
            scale_x_log10(limits= c(1, 1e2))+
            xlab('Degree') + 
            ylab('Count')+ 
            theme(legend.title = element_text(face = "bold"))

sub.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets), size=0.6)+
            geom_point()+
            #xlim (1, 10) +
            scale_x_continuous(breaks = seq(1, 10, by = 2), limits = c(1, 10))+
            ylim (0, 700) +
            theme(legend.position = "none",
                  axis.title = element_blank()
                  )

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(sub.plot, x = 0.3, y = .5, width = .5, height = .4)


ggsave(filename = "../Results/degree_plot_dtis.pdf", 
       plot = plot.with.inset,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 330)


######

for(i in 1:3) 
{                                   # Head of for-loop
  new <- rep(i, nrow(data))                       # Create new column
  data[ , ncol(data) + 1] <- new                  # Append new column
  colnames(data)[ncol(data)] <- paste0("new", i)  # Rename column name
}

df<-data.frame()
df
col1 = c('A', 'B', 'C', 'D', 'E'))
names(df) <- "DrugBank"
,
               col2 = c(1, 2, 3, 4, 5))
  
vec2 = c(id = list_names[1], numobs = 10)


data.frame(str(list_names[1]): c(1,2,3,4))

tabulate(data$DrugBank$list_number_degrees)[1]
tabulate(data$BIOSNAP$list_number_degrees)[1]


pdf(file = "testR.pdf")
p <- data[[name]][['list_number_degrees']] #[['list_subgraphs_size']]
plot(tabulate(p), type='b', xlim=xl, col="#f47823", ylab="Count",yaxt="n", xlab="Degree", log='x')
dev.off()




data[[NR]]$list_subgraphs_size

degrees <- data$list_number_degrees

unique_degrees <- unique(degrees)

dataf <- data.frame(degrees)
#head(data)

subgraphs <- data$list_subgraphs_size

pdf(file = "testR.pdf")
ggplot(dataf, aes(degrees)) +
  theme(plot.title = element_text(hjust = 0.5)) +     
  ggtitle("Histogram of Degrees for X Dataset")+
  labs(x = "Degrees") +             
  geom_histogram(bins= 10, colour="black", pad= TRUE, fill="#25b5d5") +
  #geom_density(aes(y = ..count.. * 0.8))+
  geom_freqpoly(bins = 10) 
  
  #xlim(1, max(dataf))
dev.off()



pdf(file= "test2.pdf")
plot(tabulate(degrees), type='o')
#hist(degree)
dev.off()


# Loading all data in a json

all_data <- fromJSON('test_all_datasets_info.json')
all_data <- rjson::fromJSON(file = 'test_all_datasets_info.json')


for (db in ls(all_data))
{
  #print(all_data[db])
  all_data[db]
}

library(tibble) # nicer dataframes
weather <- tibble(place = all_data)

degrees_nr <- all_data$NR$list_number_degrees
degrees_e <- all_data$E$list_number_degrees
degrees_gpcr <- all_data$GPCR$list_number_degrees
degrees_ic <- all_data$IC$list_number_degrees

pdf(file= "test2.pdf")
plot(tabulate(degrees_nr), type='o') 
#hist(degree)
dev.off()



dtinet <- rjson::fromJSON(file = 'test_dtinet_info.json')

xl = c(1,1e3)
pdf(file= "test_degree.pdf")
plot(tabulate(dtinet$DrugBank$list_number_degrees), 
      type='l', xlim=xl, col="#f47823", ylab="Count",yaxt="n", xlab="Degree", log='x') 
par(new=TRUE)
plot(tabulate(dtinet$IC$list_number_degrees), 
      type='l', xlim=xl, col="#237af4", ylab="Count",yaxt="n",  xlab="Degree", log='x') 
par(new=TRUE)
plot(tabulate(dtinet$NR$list_number_degrees), 
      type='l', xlim=xl, col="#d123f4", ylab="Count",yaxt="n",  xlab="Degree", log='x') 
par(new=TRUE)
plot(tabulate(dtinet$Davis_et_al$list_number_degrees), 
      type='l', xlim=xl, col="#42124c", ylab="Count",yaxt="n",  xlab="Degree",log='x') 
legend("topright", legend=c("DrugBank", "IC", "NR"),
       col=c("#f47823","#237af4","#d123f4", "#42124c"), lty=1, cex=0.8)
dev.off()

plot(tabulate(dtinet$DrugBank$list_number_degrees), 
      type='o', col="#f47823", ylab="",yaxt="n", xlab="Degree", log='x') 




