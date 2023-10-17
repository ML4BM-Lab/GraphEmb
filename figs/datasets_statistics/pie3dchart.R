library(rAmCharts4)
library(ggplot2)
library(webshot)


dat <- data.frame(
  group = c("A", "B", "0", "AB"),
  FR = c(20, 32, 32, 16)
)

amPieChart(
  data = dat,
  category = "group",
  value    = "FR",
  threeD = TRUE,
  variableDepth = TRUE
)

