% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xEnrichBarplot.r
\name{xEnrichBarplot}
\alias{xEnrichBarplot}
\title{Function to visualise enrichment results using a barplot}
\usage{
xEnrichBarplot(
eTerm,
top_num = 10,
displayBy = c("fc", "adjp", "fdr", "zscore", "pvalue"),
FDR.cutoff = 0.05,
bar.label = TRUE,
bar.label.size = 3,
bar.color = "lightyellow-orange",
bar.width = 0.8,
wrap.width = NULL,
font.family = "sans",
signature = TRUE
)
}
\arguments{
\item{eTerm}{an object of class "eTerm"}

\item{top_num}{the number of the top terms (sorted according to FDR or
adjusted p-values). If it is 'auto', only the significant terms (see
below FDR.cutoff) will be displayed}

\item{displayBy}{which statistics will be used for displaying. It can
be "fc" for enrichment fold change (by default), "adjp" or "fdr" for
adjusted p value (or FDR), "pvalue" for p value, "zscore" for
enrichment z-score}

\item{FDR.cutoff}{FDR cutoff used to declare the significant terms. By
default, it is set to 0.05. This option only works when setting top_num
(see above) is 'auto'}

\item{bar.label}{logical to indicate whether to label each bar with
FDR. By default, it sets to true for bar labelling}

\item{bar.label.size}{an integer specifying the bar labelling text
size. By default, it sets to 3}

\item{bar.color}{either NULL or fill color names ('lightyellow-orange'
by default)}

\item{bar.width}{bar width. By default, 80% of the resolution of the
data}

\item{wrap.width}{a positive integer specifying wrap width of name}

\item{font.family}{the font family for texts}

\item{signature}{logical to indicate whether the signature is assigned
to the plot caption. By default, it sets TRUE showing which function is
used to draw this graph}
}
\value{
an object of class "ggplot"
}
\description{
\code{xEnrichBarplot} is supposed to visualise enrichment results using
a barplot. It returns an object of class "ggplot".
}
\note{
none
}
\examples{
\dontrun{
# Load the XGR package and specify the location of built-in data
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"

# 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
cis <- xRDataLoader(RData.customised='JKscience_TS2A',
RData.location=RData.location)
ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
data <- df_cis$variant

# 2) Enrichment analysis using Experimental Factor Ontology (EFO)
# Considering LD SNPs and respecting ontology tree
eTerm <- xEnricherSNPs(data, ontology="EF", include.LD="EUR",
LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)

# 3) Barplot of enrichment results
bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc")
#pdf(file="enrichment_barplot.pdf", height=6, width=12, compress=TRUE)
print(bp)
#dev.off()
}

# 4) use font family (Arial)
\dontrun{
BiocManager::install("extrafont")
library(extrafont)
font_import()
fonttable()
## creating PDF files with fonts
library(extrafont)
loadfonts()
bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc",
font.family="Arial Black")
pdf(file="enrichment_barplot_fonts.pdf", height=6, width=12,
family="Arial Black")
print(bp)
dev.off()
}
}
\seealso{
\code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}},
\code{\link{xEnrichViewer}}
}
