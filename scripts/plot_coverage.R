
# AY271541.1  4   14373

#input_file = "assemblies/bwa/ITS/33/33_read_depth.tsv"
input_file = "assemblies/bwa/CYTB/33/33_read_depth.tsv"
#input_file = "assemblies/chloroplast-bwa/33/33_read_depth.tsv"

library(ggplot2)

data = read.table(input_file, sep="\t")
names(data) = c("reference", "base", "depth")

print(ggplot(data, aes(base, depth)) + geom_line(aes(colour=depth)))
