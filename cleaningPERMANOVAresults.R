results <- read_csv("permanovaResults.csv")

sevYear <- c("U.2020_vs_L.2020", "U.2020_vs_H.2020", "L.2020_vs_H.2020",
             "U.2021_vs_L.2021", "U.2021_vs_H.2021", "L.2021_vs_H.2021",
             "U.2022_vs_L.2022", "U.2022_vs_H.2022", "L.2022_vs_H.2022",
             "U.2023_vs_L.2023", "U.2023_vs_H.2023", "L.2023_vs_H.2023",
             "U.2024_vs_L.2024", "U.2024_vs_H.2024", "L.2024_vs_H.2024")

permResult <- results %>%
  separate(element, into = c("sev.year", "element"), sep = ":") %>%
  pivot_wider(names_from = element, values_from = x) %>%
  filter(sev.year %in% sevYear)

write.csv(permResult, "PERMresults_repeatMeasure.csv")


envfit <- read_csv("envfit_results.txt")
write.csv(envfit, "envfit_results.csv")
