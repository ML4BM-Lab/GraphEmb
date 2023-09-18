
library(rjson)
library(ggplot2)
library(devtools)
library(circlize)
library(RColorBrewer)
library(reshape2)





data <- fromJSON(file='../Results/statistics_extended_eegdti.json')

list_names <- names(data)

colors <- brewer.pal(length((list_names)), "Paired")



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

df <- melt(gfg_data ,  id.vars = 'index', variable.name = 'datasets')


library(cowplot)


main.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets), size=0.8)+
            geom_point(aes(colour = datasets), alpha = 0.5)+
            scale_color_manual(values = colors) +
            scale_x_log10(limits= c(1, 1e2))+
            xlab('Degree (log10 scale)') + 
            ylab('Count')+ 
            theme(legend.title = element_text(face = "bold"))+
            theme_classic()

sub.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets), size=0.6, show.legend = FALSE)+
            geom_point(aes(colour = datasets), show.legend = FALSE, alpha = 0.5)+
            scale_color_manual(values = colors) +
            #xlim (1, 10) +
            scale_x_continuous(breaks = seq(1, 10, by = 2), limits = c(1, 10))+
            ylim (0, 700) +
            theme_classic()+
            theme(axis.title.y = element_blank())+
            xlab("Degree")

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(sub.plot, x = 0.3, y = .5, width = .5, height = .4)


ggsave(filename = "../Results/degree_plot_ext_eegdti.pdf", 
       plot = plot.with.inset,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 330)




#########################
# subgraphs

list_names <- names(data)
colors <- brewer.pal(length((list_names)), "Paired")


# Creating variable
name_var <- 'list_subgraphs_size'
x1= tabulate(data[[list_names[1]]][['list_subgraphs_size']])
x2= tabulate(data[[list_names[2]]][['list_subgraphs_size']])
x3= tabulate(data[[list_names[3]]][['list_subgraphs_size']])
x4= tabulate(data[[list_names[4]]][['list_subgraphs_size']])
x5= tabulate(data[[list_names[5]]][['list_subgraphs_size']])
x6= tabulate(data[[list_names[6]]][['list_subgraphs_size']])
x7= tabulate(data[[list_names[7]]][['list_subgraphs_size']])
x8= tabulate(data[[list_names[8]]][['list_subgraphs_size']])


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

df <- melt(gfg_data ,  id.vars = 'index', variable.name = 'datasets')


xlim_ = 200

main.plot <- ggplot(df, aes(index, value)) +
            geom_line(aes(colour = datasets), size=0.8)+
            geom_point(aes(colour = datasets), alpha = 0.5)+
            scale_color_manual(values = colors) +
            scale_x_log10(limits= c(1, 1e4))+
            #scale_x_continuous(breaks = seq(1, 10, by = 1), limits = c(1, 10))+
            # scale_x_continuous(limits =c(0, 2000),
            #          breaks = seq(0, 2000, by = 200),
            #          #labels = c(0, seq(1600, 2000, by = 200)),
            #          expand = c(0,0,0.05,0)) +
            xlab('# nodes per connected component') + 
            ylab('Count')+ 
            theme_classic()+
            theme(legend.title = element_text(face = "bold"))


ggsave(filename = "../Results/subgraph_plot_ext_eegdti.pdf", 
       plot = main.plot,
       width = 20, 
       height = 20,
       units = "cm",
       dpi = 330)



#