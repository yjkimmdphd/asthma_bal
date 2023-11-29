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
revigo.data <- rbind(c("GO:0003341","cilium movement",0.9515229998310906,6.183758700008217,0.5606415255568971,0.00519511,"cilium movement"),
c("GO:0003351","epithelial cilium movement involved in extracellular fluid movement",0.23647317155565564,6.384049948343599,0.540175322220431,0.69846878,"cilium movement"),
c("GO:0007018","microtubule-based movement",2.1901919936940484,2.777492404504671,0.5769173010610327,0.68842059,"cilium movement"),
c("GO:0006954","inflammatory response",3.1135634254827993,2.1886338810240797,0.9724137177714235,-0,"inflammatory response"),
c("GO:0007155","cell adhesion",5.360058555261528,2.850718673726238,0.9903854979023881,0.00649384,"cell adhesion"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.9402623726141547,4.141462802430361,0.9921105448155548,0.00518796,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0007368","determination of left/right symmetry",0.715049828275435,3.237321436272564,0.9780572647702253,-0,"determination of left/right symmetry"),
c("GO:0001775","cell activation",4.031304543663082,1.6810397923089273,0.9677721794126372,0.1511869,"determination of left/right symmetry"),
c("GO:0032731","positive regulation of interleukin-1 beta production",0.3434491301165475,3.237321436272564,0.9225078775240496,-0,"positive regulation of interleukin-1 beta production"),
c("GO:0002224","toll-like receptor signaling pathway",0.24210348516412367,3.237321436272564,0.8811132846534178,0.1183438,"positive regulation of interleukin-1 beta production"),
c("GO:0007166","cell surface receptor signaling pathway",11.727943246438826,1.4527415505596304,0.9332579967761274,0.23593378,"positive regulation of interleukin-1 beta production"),
c("GO:0032760","positive regulation of tumor necrosis factor production",0.5686616744552672,2.1886338810240797,0.9214854955437136,0.64678495,"positive regulation of interleukin-1 beta production"),
c("GO:0051649","establishment of localization in cell",10.027588536681494,1.7793221304913946,0.9450965213424619,0.00859243,"establishment of localization in cell"),
c("GO:0060271","cilium assembly",1.8523731771859693,9.924453038607469,0.6918724544441087,0,"cilium assembly"));

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
