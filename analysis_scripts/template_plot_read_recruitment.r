PROJPATH="" #change this to your pipeline results
POGS=paste(PROJPATH, "/completed_pogs.txt",sep="")
OUTDIR= paste(PROJPATH, "/analysis", sep="")
pog_filter = "No POG filter thresholds"
WIDTH = 15
HEIGHT = 10

library(ggplot2)
options(repr.plot.width=30, repr.plot.height=10) #if there's many samples you may need to expand width
theme_set(theme_classic())

res = read.delim(paste(OUTDIR,"select_pogs_read_recruitment.tsv",sep="/"))

#the color palette
pal =  defaultColors = c(
        '#97B065', '#B93838', '#83AAF0', '#E2AA48', '#A9A3D2',
        '#797900', '#F4BCE6', '#84F4F6', '#9AD63B', '#A285F5',
        '#D2B48C', '#4EA24E', '#465569', '#9C669C', '#6495ED')

plot.subtitle = pog_filter
a = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=mOTU)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to mOTU/genome.", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

a

ggsave(plot=a, filename=paste(OUTDIR, "/mOTUs_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

b = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Domain)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Domain", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

b

ggsave(plot=b, filename=paste(OUTDIR, "/Domain_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

c = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Phylum)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Phylum", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

c

ggsave(plot=c, filename=paste(OUTDIR, "/Phylum_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

d = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Class)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Class", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

d

ggsave(plot=d, filename=paste(OUTDIR, "/Class_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

e = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Order)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Order", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

e

ggsave(plot=e, filename=paste(OUTDIR, "/Order_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

f = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Genus)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Genus", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

f

ggsave(plot=f, filename=paste(OUTDIR, "/Genus_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)

g = ggplot(res, aes(x=Sample, y=prcnt_mapped, fill=Species)) + 
  geom_bar(stat = "identity", width=0.8) + scale_fill_manual(values = pal) +
  labs(title="Summarized % of reads mapped per sample to Species", subtitle=plot.subtitle) +
  theme(text=element_text(size=15), panel.grid.major = element_line(colour = "red", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "red", linetype = "dotted"))

g

ggsave(plot=g, filename=paste(OUTDIR, "/Species_read_recruitment.png", sep=""), width=WIDTH, height=HEIGHT)
