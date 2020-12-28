#Helper file to plot bar charts

library(ggplot2)

path="/Users/carolinebao/Documents/Feng UROP (2019-)/Gene Therapy/shank3/data/npy/mouse/UCSC_mm10/results/intersections/gene_intersect_frequencies_h_m_up_20000_down_0_window_10_bps.txt"

v <- read_table_helper(path)
print(v["V2"])

ggplot(data = v, mapping = aes(x=V2)) + 
  geom_histogram(aes(y=V6),fill="bisque",color="white",alpha=0.7) + 
  geom_density() +
  geom_rug() +
  labs(x='mean education per house') +
  theme_minimal() + stat=identity



function windows_plotting(folder)
data <- read_table_helper(path, fns, graph_names) {

  # Create data
  data <- data.frame(
    name=c("A","B","C","D","E") ,  
    value=c(3,12,5,18,45)
  )
  
  # Barplot
  ggplot(aes(x=data)
}

# Use datn2 from above
ggplot(data=v, aes(x=factor(V2), y=V6)) +
  geom_bar(stat="identity", position=position_dodge())

# Use the original data frame, but put factor() directly in the plot specification
ggplot(data=datn, aes(x=factor(V2), y=length, fill=supp)) +
  geom_bar(stat="identity", position=position_dodge())