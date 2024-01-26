library("tidyverse")

the_data <- read_delim("summary.csv",delim=";")

    ggplot(data=the_data
           ,mapping = aes(x=baseline_mortality, y=Bfec)) +
        geom_tile(mapping=aes(fill=z_large))
