rm(list=ls())
setwd("/home/hyunseok/Desktop/YS/Reassortment")
raw <- read.csv("./data/ultimate_test/winning_prob.csv",header=T)
library(ggplot2)
library(reshape2)
mu_values = c(0.00067,0.0008,0.0010,0.0013,0.0033,0.0067,0.01)
mat_p1 = subset(raw, mu == mu_values[3])
mat_p2 = subset(raw, mu == mu_values[4])

mat = mat_p1
p1 <- ggplot(mat,aes(c,y=1)) +
  geom_tile(aes(fill = prob), color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  ylab("") +
  xlab("c") +
  ggtitle(sprintf("mu=%.5f",unique(mat$mu))) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Probability") +
  facet_grid(n1r~s)

mat = mat_p2
p2 <- ggplot(mat,aes(c,y=1)) +
  geom_tile(aes(fill = prob), color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  ylab("") +
  xlab("c") +
  ggtitle(sprintf("mu=%.5f",unique(mat$mu))) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Probability") +
  facet_grid(n1r~s)





# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(p1,p2)
