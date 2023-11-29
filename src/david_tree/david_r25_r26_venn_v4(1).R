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
revigo.data <- rbind(c("GO:0003341","cilium movement",0.9515229998310906,6.31069114087638,0.561071686353571,0.00519511,"cilium movement"),
c("GO:0003351","epithelial cilium movement involved in extracellular fluid movement",0.23647317155565564,6.661543506395395,0.5405257741399334,0.69846878,"cilium movement"),
c("GO:0007018","microtubule-based movement",2.1901919936940484,2.884173892174819,0.577414665838287,0.68842059,"cilium movement"),
c("GO:0006954","inflammatory response",3.1135634254827993,1.9773139743742745,0.9853570058750698,-0,"inflammatory response"),
c("GO:0007155","cell adhesion",5.360058555261528,3.8477116556169433,0.9911521884472024,0.00649384,"cell adhesion"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.9402623726141547,4.244887733604929,0.9926817439821877,0.00518796,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0007368","determination of left/right symmetry",0.715049828275435,3.3736596326249577,0.9492894604770048,-0,"determination of left/right symmetry"),
c("GO:0001775","cell activation",4.031304543663082,1.7306066900341266,0.9486999096750758,0.1511869,"determination of left/right symmetry"),
c("GO:0007399","nervous system development",12.358538370587242,1.386481060283993,0.9397279645155444,0.33701693,"determination of left/right symmetry"),
c("GO:0032731","positive regulation of interleukin-1 beta production",0.3434491301165475,3.3736596326249577,0.925469279436007,-0,"positive regulation of interleukin-1 beta production"),
c("GO:0002224","toll-like receptor signaling pathway",0.24210348516412367,3.3736596326249577,0.9004763770178037,0.1183438,"positive regulation of interleukin-1 beta production"),
c("GO:0032760","positive regulation of tumor necrosis factor production",0.5686616744552672,2.2680261002136812,0.9246495185594399,0.64678495,"positive regulation of interleukin-1 beta production"),
c("GO:0051649","establishment of localization in cell",10.027588536681494,1.8528344961749192,0.945951090323669,0.00859243,"establishment of localization in cell"),
c("GO:0060271","cilium assembly",1.8523731771859693,9.387216143280265,0.6924022287939741,0,"cilium assembly"),
c("GO:0060285","cilium-dependent cell motility",0.7544620235347109,4.841637507904751,0.8308551521500195,0.00505928,"cilium-dependent cell motility"));

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
