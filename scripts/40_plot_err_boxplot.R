library(ggplot2)
dat = read.table("analysis/4-admix/all.err", header = T)
print(dat)

p = ggplot(dat, aes( x = k, y = err)) +
geom_point() +
theme_classic()


options(bitmapType='cairo')
png("analysis/4-admix/all.err.png")
plot(p)
dev.off()