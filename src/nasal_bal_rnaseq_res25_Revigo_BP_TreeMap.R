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
revigo.data <- rbind(c("GO:0003351","epithelial cilium movement involved in extracellular fluid movement",0.004839592939343236,4.581698708680254,0.6776831232616413,0.00599898,"epithelial cilium movement involved in extracellular fluid movement"),
c("GO:0003341","cilium movement",0.05743316239425406,4.023650020996727,0.6933241078864577,0.6263572,"epithelial cilium movement involved in extracellular fluid movement"),
c("GO:0006954","inflammatory response",0.11987887912499628,7.023191662661934,0.8088102355060482,-0,"inflammatory response"),
c("GO:0001782","B cell homeostasis",0.0036155584364715447,2.1650148261381545,0.8084821747692732,0.57225905,"inflammatory response"),
c("GO:0002224","toll-like receptor signaling pathway",0.03528745119827656,4.512861624522814,0.6175217941849953,0.22315369,"inflammatory response"),
c("GO:0006952","defense response",1.0410213876080876,2.5630930091220585,0.8121565751104448,0.50776052,"inflammatory response"),
c("GO:0007165","signal transduction",8.135788134528843,1.6292014468287996,0.7269472786373754,0.65912233,"inflammatory response"),
c("GO:0007166","cell surface receptor signaling pathway",1.1935334257213288,3.4056074496245734,0.7631827392424497,0.37240978,"inflammatory response"),
c("GO:0032731","positive regulation of interleukin-1 beta production",0.006116846333644131,4.183758700008217,0.7679111173638482,0.43054164,"inflammatory response"),
c("GO:0007155","cell adhesion",0.7060783158562494,3.8416375079047502,0.990887624109389,0.00881469,"cell adhesion"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.216251639141856,4.337242168318426,0.9916650366284384,0.00782319,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0007368","determination of left/right symmetry",0.025431977741730852,2.280971581526965,0.9179614525972312,-0,"determination of left/right symmetry"),
c("GO:0001775","cell activation",0.09565097880185117,1.919569422627293,0.9070315719105528,0.4799785,"determination of left/right symmetry"),
c("GO:0034440","lipid oxidation",0.1496016299882067,1.9849217049426684,0.9918821766682623,0.0075991,"lipid oxidation"),
c("GO:0051092","positive regulation of NF-kappaB transcription factor activity",0.020316311803098615,2.3846266914355088,0.8992874438329446,-0,"positive regulation of NF-kappaB transcription factor activity"),
c("GO:0060271","cilium assembly",0.17316429416848675,9.170053304058364,0.7847363449071298,0,"cilium assembly"),
c("GO:0060285","cilium-dependent cell motility",0.04143755933906394,3.493494967595128,0.8973642847742443,0.00690953,"cilium-dependent cell motility"));

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

