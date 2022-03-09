
plotZcurve = function(coord_input, outputname, format_list) {
  step=seq(1,nrow(coord_input))
  print('Plotting the Z-curve...\n')
  #png(file=outputname, width=600, height=350)
  lines3D(coord_input[,'X'], coord_input[,'Y'], coord_input[,'Z'], 
          theta = 15, phi = 20,
          colvar=step,
          main = "Z-curve", 
          xlab = "R/Y disparity",
          ylab ="M/K disparity", 
          zlab = "W/S disparity")
  for (file_format in format_list) {
    dev.copy(
      eval(parse(text = file_format)),
      paste(outputname, file_format, sep = ".")
    )
    dev.off()
  }
}

