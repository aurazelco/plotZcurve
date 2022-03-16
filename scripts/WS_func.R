# defines the function, which takes as inputs a dataframe containing the coordinates for the 3 axes,
# the output filename and a list containing the formats of the output plots
plotWS = function(coord_input, ws_outputname, format_list, plot_title) {
  # creates a vector with an integer step-count of the genome sequence -> used for x-axis
  step=seq(1,nrow(coord_input))
  # plots the W/S disparity, corresponding to the GC content and saves it to an object
  # x-axis is the seq length, on the y-axis is the W/S disparity, from the Z column of the coord dataframe
  WS_plot <- ggplot(data=NULL, aes(step, coord_input[,'Z'], color=coord_input[,'Z'])) +
    # line
    geom_line() +
    # labels for x and y axes, legend and main plot title
    xlab('Sequence length') +
    ylab('W/S disparity') +
    labs(colour='W/S values') +
    ggtitle(plot_title) + 
    # empty background
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),    
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
  # prints the plot, otherwise it cannot be saved
  print(WS_plot)
  # for each format present in the input list
  for (file_format in format_list) {
    # creates a new file name
    file_output = paste(ws_outputname, file_format, sep=".")
    # saves the plot in the correct format
    ggsave(file_output, plot = WS_plot, dpi=300)
  }
}