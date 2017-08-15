
# AY271541.1  4   14373

input_file = "data/2_read_depth_chloro.tsv"

library(ggplot2)

data = read.table(input_file, sep="\t")
names(data) = c("reference", "base", "depth")

p = ggplot(data, aes(base, depth)) + ylim(0, 500) + geom_line(aes(colour=depth))

print(p)
