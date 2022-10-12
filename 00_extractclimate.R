library(terra)
library(ggplot2)
library(dplyr)

t2m <- rast("./data/b30.00_4kaDVTd.cam2.ncrcat.jja.TS.nc")
t2m <- rotate(t2m)
t2m.ts <- extract(t2m, cbind(-110, 40))

pre <- rast("./data/b30.00_4kaDVTd.cam2.ncrcat.jja.PRECT.nc")
pre <- rotate(pre)
pre <- pre[[1:2204]]
pre.ts <- extract(pre, cbind(-110, 40))

dates <- names(t2m.ts)
dates <- unlist(strsplit(dates, "_"))
dates_id <- seq(2, length(dates), by = 2)
dates <- as.numeric(dates[dates_id]) * 10


clim.df <- data.frame(age = dates,
                      t2m = t(t2m.ts),
                      pre = t(pre.ts))

clim.df <- clim.df %>%
  mutate(age = 22000 - age,
         t2m = t2m - 273.15,
         pre = pre * (60 * 60 * 24 * 30 * 12 * 100))

ggplot(clim.df, aes(x = age, y = t2m)) +
  geom_line()

ggplot(clim.df, aes(x = age, y = pre)) +
  geom_line()

write.csv(clim.df, "./data/clim.df.csv", row.names = FALSE)
