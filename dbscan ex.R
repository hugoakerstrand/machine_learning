set.seed(123)
data <- tibble(
  FSC_A = rnorm(1000, mean = 100, sd = 20),
  FSC_H = rnorm(1000, mean = 100, sd = 15),
  SSC_A = rnorm(1000, mean = 80, sd = 15),
  SSC_H = SSC_A + rnorm(1000, mean = 0, sd = 5),
  Marker_Intensity = c(rnorm(800, mean = 3, sd = 0.5), rnorm(200, mean = 6, sd = 0.7)) # On-target vs off-target
)

data_recipe <- recipe(~ FSC_A + SSC_A, data = data) %>%
  step_normalize(all_numeric_predictors())

prepped_data <- prep(data_recipe) %>%
  bake(new_data = NULL)

dbscan_result <- dbscan::dbscan(prepped_data, eps = 0.5, minPts = 10)

tidy(dbscan_result)
augment(dbscan_result, data) |> 
  ggplot(aes(x = FSC_A, y = SSC_A, colour = `.cluster`)) +
  geom_point()

dbscan_result |> tail()
?map2
