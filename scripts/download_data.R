
# Load libraries
library(SomaDataIO)
library(tidyverse)

# Download data
download.file(
  url = "https://raw.githubusercontent.com/SomaLogic/SomaLogic-Data/main/example_data.adat",
  destfile = "data/example_data.adat"
)

# Read in data
my_adat <- read_adat("data/example_data.adat")

# Check the content
is.soma_adat(my_adat)
methods(class = "soma_adat")
my_adat |> View() # By default sample column
my_adat # Explore object

# Number of somamers (analytes) -> 5K
my_adat |>
  getAnalytes() |>
  length()

