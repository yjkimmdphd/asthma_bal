# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003341","cilium movement",0.05743316239425406,6.31069114087638,0.5625728737840002,0.00707273,"cilium movement"),
c("GO:0003351","epithelial cilium movement involved in extracellular fluid movement",0.004839592939343236,6.661543506395395,0.5393822082222656,0.6263572,"cilium movement"),
c("GO:0006954","inflammatory response",0.11987887912499628,1.9773139743742745,0.9640828838890211,-0,"inflammatory response"),
c("GO:0007155","cell adhesion",0.7060783158562494,3.8477116556169433,0.9879697346395311,0.00864154,"cell adhesion"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.216251639141856,4.244887733604929,0.9889401879342715,0.00881469,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0007368","determination of left/right symmetry",0.025431977741730852,3.3736596326249577,0.8734628498264415,-0,"determination of left/right symmetry"),
c("GO:0001775","cell activation",0.09565097880185117,1.7306066900341266,0.865370040722323,0.4799785,"determination of left/right symmetry"),
c("GO:0007399","nervous system development",0.5207002861047604,1.386481060283993,0.8580729381832477,0.64248608,"determination of left/right symmetry"),
c("GO:0032731","positive regulation of interleukin-1 beta production",0.006116846333644131,3.3736596326249577,0.8669363771477541,-0,"positive regulation of interleukin-1 beta production"),
c("GO:0002224","toll-like receptor signaling pathway",0.03528745119827656,3.3736596326249577,0.828905744619872,0.43054164,"positive regulation of interleukin-1 beta production"),
c("GO:0051649","establishment of localization in cell",2.1876556652574304,1.8528344961749192,0.9549898417108053,0.01113915,"establishment of localization in cell"),
c("GO:0060271","cilium assembly",0.17316429416848675,9.387216143280265,0.666635236794599,0,"cilium assembly"),
c("GO:0060285","cilium-dependent cell motility",0.04143755933906394,4.841637507904751,0.8670189594150165,0.00690953,"cilium-dependent cell motility"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
