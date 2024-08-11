###

# Extract image data from the MSImagingExperiment object
image_data <- data_matrix

# Calculate the mean metabolic content across all pixels for each mz value
average_metabolic_content <- apply(image_data, 1, mean)

# Create a data frame for plotting
average_metabolic_df <- data.frame(
  mz = fdata@mz,
  average_content = average_metabolic_content
)

# Plotting with ggplot2
ggplot(average_metabolic_df, aes(x = mz, y = average_content)) +
  geom_line() +
  labs(title = "Average Metabolic Content by mz", x = "mz", y = "Average Metabolic Content") +
  theme_minimal()

# Alternatively, plotting with base R
plot(average_metabolic_df$mz, average_metabolic_df$average_content, type = "l",
     main = "Average Metabolic Content by mz",
     xlab = "mz", ylab = "Average Metabolic Content")


###############
dim(image_data)

# Calculate the mean metabolic content across all pixels for each mz value
average_spot_content <- apply(image_data, 2, mean)
# Create a data frame for plotting
average_spot_df <- data.frame(
  spot = colnames(data_matrix),
  average_content = average_spot_content
)
######################
# Add mean metabolic content as a new variable to the PositionDataFrame
pdata$average_spot_content <- as.numeric(average_spot_df$average_content)
pdata
# Update the MSImagingExperiment object with the new PositionDataFrame
msi_data <- MSImagingExperiment(
  imageData = idata,
  featureData = fdata,
  pixelData = pdata
)

# Print updated object information
print(msi_data)

# Plot the image with the new variable mean_metabolic_content
image(msi_data, 
      #feature =1,
      groups="average_spot_conten",
      asp=1.5,
      #superpose = TRUE,
      col = '#d62728',
      main = "Mean Metabolic Content per Spot", xlab = "x", ylab = "y")

image(msi_data, 
      #i=3,
      #feature =182:210,
      #groups="average_spot_conten",
      asp=1.5,
      #superpose = TRUE,
      col = '#d62728',
      #main = "Mean Metabolic Content per Spot", 
      xlab = "x", ylab = "y")



