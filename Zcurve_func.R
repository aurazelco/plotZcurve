# defines the function, which takes as inputs a dataframe containing the coordinates for the 3 axes,
# the output filename and a list containing the formats of the output plots
plotZcurve = function(coord_input, outputname, format_list, plot_title) {
  # creates a vector with an integer step-count of the genome sequence -> used for plot legend
  step=seq(1,nrow(coord_input))
  # creates a 3D plot, with the points represented as a line, from the 3 columns of the dataframe
  lines3D(coord_input[,'X'], coord_input[,'Y'], coord_input[,'Z'], 
          # colors the line: here I chose to color it according to the step, so we can know the direction
          colvar=step,
          # legend title
          clab = c("Sequence length"),
          # main title
          main = plot_title, 
          # axes titles, see README for more details
          xlab = "R/Y disparity",
          ylab ="M/K disparity", 
          zlab = "W/S disparity")
  # for each format present in the input list
  for (file_format in format_list) {
    # because of the different formats, we have to adjusts the size, since using units or res raises an error
    if (file_format == 'pdf') {
        # pdf needs a smaller scale to be visible, otherwise the plot created is too big
        param_plot = c(20,10)
      # if not a pdf, all other vectorial image should be fine
      } else {
        param_plot = c(700,350)
      }
    # copies the plot without having to re-enter the commands
    dev.copy(
      width = param_plot[1],
      height = param_plot[2],
      # takes file_format as a string
      eval(parse(text = file_format)),
      # creates the new filename with the file format after the dot -> it is the new outputname, including the file format
      paste(outputname, file_format, sep = ".")
    )
    # closes the plots and automatically saves it
    dev.off()
  }
}
