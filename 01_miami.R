library(ggplot2)
library(dplyr)

clim.df <- read.csv("./data/clim.df.csv")
clim.df <- clim.df %>% 
  mutate(t2m = t2m,
         pre = pre)

k1 = 3000
k2 = 1.315
k3 = 0.119
k4 = 3000
k5 = 0.000664
tmp <- seq(-20,30)

## 0.45 is the carbon:dry matter ratio (King et al 1997)
f1t = 0.45 * (k1 / (1 + exp(k2 - k3 * tmp)))
plot.df = data.frame(tmp = tmp, f1t = f1t)
ggplot(plot.df, aes(x = tmp, y = f1t)) +
  geom_line()
pre <- seq(0,4500, length.out = 200)
f2p = 0.45 * k4 * (1 - exp(-(4000 * k5 * pre)))
f2p = 0.45 * k4 * (1 - exp(-k5 * pre))
plot.df = data.frame(pre = pre, f2p = f2p)
ggplot(plot.df, aes(x = pre, y = f2p)) +
  geom_line()




miami <- function(tmp, pre, 
                  k1 = 3000, k2 = 1.315, k3 = 0.119, 
                  k4 = 3000, k5 = 0.0006644) {
  
  f1t = 0.45 * (k1 / (1 + exp(k2 - k3 * tmp)))
  f2p = 0.45 * k4 * (1 - exp(-k5 * pre))
  # f2p = 0.45 * k4 * (1 - exp(-4000 * k5 * pre))
  # 
  npp = (1 / 1000) * apply(cbind(f1t, f2p), 1, min) 
  return(npp)
}

clim.df$npp <- miami(clim.df$t2m, clim.df$pre)
ggplot(clim.df, aes(x = age, y = t2m)) +
  geom_line()
ggplot(clim.df, aes(x = age, y = pre)) +
  geom_line()
ggplot(clim.df, aes(x = age, y = npp)) +
  geom_line()

write.csv(clim.df, "data/clim.df.csv", row.names = FALSE)
