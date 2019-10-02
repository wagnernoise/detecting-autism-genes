# DISCOVERING GENES THAT CAUSE AUTISM USING NETWORK THEORY


## Intro
This analysis uses the data provided by the BIOGRID dataset (http://thebiogrid.org). The dataset used in this study is a simplified version provided by MIT 'BIOGRID.txt'

## The code
The whole analysis of dectecting graphical network correlation between genes that can cause autism is divided into 3 codes (for a better understanding).

The first part ('protein-interaction.R') includes the computation of the graphical network of the data and the  removal of multiples edges and self nodes.

```{r}
library("png")
pp <- readPNG("figs/fig1.png")
plot.new() 
rasterImage(pp,0,0,1,1)
```

The second code ('gene-autism.R') is used to identify the known genes that cause autism. Dataset for this code are “gene-id-table.txt” and “gene-score.csv”. The output graphic is saved as "ASD_interactome.pdf"

Finally, the third and last code make the analysis of autism interactome ('austism-interactome.R') and indentify new candidate genes.

