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
revigo.data <- rbind(c("GO:0000281","mitotic cytokinesis",0.03798165757688001,1.3384015624748078,0.9938329515387426,0.00686726,"mitotic cytokinesis"),
c("GO:0006468","protein phosphorylation",4.340892012483024,10.331614083309999,0.9756039451485845,0.01029427,"protein phosphorylation"),
c("GO:0018105","peptidyl-serine phosphorylation",0.0408920657019146,1.8927162065650522,0.979569634091496,0.43875663,"protein phosphorylation"),
c("GO:0046777","protein autophosphorylation",0.0549352006777958,3.3477536589966768,0.9792365477376578,0.45092599,"protein phosphorylation"),
c("GO:0006909","phagocytosis",0.026872213991033685,3.5072396109731625,0.9511026809326598,-0,"phagocytosis"),
c("GO:0003341","cilium movement",0.05743316239425406,2.111573539725549,0.9194928424096498,0.6263572,"phagocytosis"),
c("GO:0003351","epithelial cilium movement involved in extracellular fluid movement",0.004839592939343236,2.6384963555559575,0.9001747345771269,0.15192879,"phagocytosis"),
c("GO:0008104","protein localization",3.045078529796193,1.6466690285355703,0.9474155299917646,0.5257779,"phagocytosis"),
c("GO:0031623","receptor internalization",0.010141524997977682,1.3901285309977716,0.9533922483986125,0.65090729,"phagocytosis"),
c("GO:0035720","intraciliary anterograde transport",0.0006020387092928699,2.142088348049649,0.8555448406496075,0.67997856,"phagocytosis"),
c("GO:0061512","protein localization to cilium",0.008428541930100179,3.428291168191312,0.962029124464016,0.15424312,"phagocytosis"),
c("GO:0006915","apoptotic process",0.25677782496519486,6.174573882232177,0.9929341211145476,0.01073986,"apoptotic process"),
c("GO:0006954","inflammatory response",0.11987887912499628,15.882728704344236,0.8144046243887921,-0,"inflammatory response"),
c("GO:0006935","chemotaxis",0.5617320513966766,4.528708288941061,0.7517864038909725,0.56134016,"inflammatory response"),
c("GO:0006952","defense response",1.0410213876080876,1.4992242862509764,0.8028013722067043,0.69823798,"inflammatory response"),
c("GO:0006968","cellular defense response",0.0012074035993000652,1.5002059649207027,0.8569198381129437,0.50052312,"inflammatory response"),
c("GO:0006974","DNA damage response",2.609069457039582,2.9945500681822694,0.7695157625348175,0.55273719,"inflammatory response"),
c("GO:0007166","cell surface receptor signaling pathway",1.1935334257213288,1.3103872144663724,0.6852006139743227,0.59865377,"inflammatory response"),
c("GO:0007173","epidermal growth factor receptor signaling pathway",0.008425215749385854,2.0955105846823385,0.7458083117188719,0.60258737,"inflammatory response"),
c("GO:0007224","smoothened signaling pathway",0.018696461795222222,1.3551654475608241,0.7373737459175868,0.49264193,"inflammatory response"),
c("GO:0007249","I-kappaB kinase/NF-kappaB signaling",0.027850111121045303,4.920818753952375,0.7487803208681156,0.4023106,"inflammatory response"),
c("GO:0007259","receptor signaling pathway via JAK-STAT",0.010324464937265571,1.649620627451385,0.7458468666938561,0.37024472,"inflammatory response"),
c("GO:0007264","small GTPase mediated signal transduction",0.33345626897253483,3.9956786262173574,0.7058580981664736,0.51359779,"inflammatory response"),
c("GO:0007265","Ras protein signal transduction",0.07726717799377551,1.8230046130184485,0.7327377894128535,0.44164372,"inflammatory response"),
c("GO:0007267","cell-cell signaling",0.420302847423567,1.895762953134955,0.798954083859648,0.52294769,"inflammatory response"),
c("GO:0019221","cytokine-mediated signaling pathway",0.14873017064105348,4.749579997691106,0.6271588857519393,0.64664597,"inflammatory response"),
c("GO:0019722","calcium-mediated signaling",0.06199668233430831,1.8991761984324766,0.7363696083280621,0.43251915,"inflammatory response"),
c("GO:0030593","neutrophil chemotaxis",0.006319743357217971,4.632644078973981,0.7350468764446034,0.69343872,"inflammatory response"),
c("GO:0035556","intracellular signal transduction",3.7944870018179575,7.168130225719498,0.655021442852773,0.32756022,"inflammatory response"),
c("GO:0042542","response to hydrogen peroxide",0.01811105398950098,2.045654719694031,0.7770390147176681,0.59569768,"inflammatory response"),
c("GO:0042832","defense response to protozoan",0.0027075111014607522,2.7358329258644205,0.8188232484404377,0.52396412,"inflammatory response"),
c("GO:0045087","innate immune response",0.16539433601982298,1.7933082094833497,0.7435803443843246,0.6880562,"inflammatory response"),
c("GO:0048015","phosphatidylinositol-mediated signaling",0.05648852707138569,1.3901285309977716,0.737874951173876,0.42877579,"inflammatory response"),
c("GO:0048678","response to axon injury",0.0074872327879461345,1.9485448122319633,0.8572185182168465,0.35333865,"inflammatory response"),
c("GO:0060397","growth hormone receptor signaling pathway via JAK-STAT",0.0005987125285785447,2.5528407278181766,0.718914254721869,0.52962494,"inflammatory response"),
c("GO:0071222","cellular response to lipopolysaccharide",0.01857671928950651,7.068542129310995,0.7451813685040436,0.21380622,"inflammatory response"),
c("GO:0071260","cellular response to mechanical stimulus",0.007340880836515823,2.662190648808808,0.8275482535580228,0.42689592,"inflammatory response"),
c("GO:1990090","cellular response to nerve growth factor stimulus",0.0023116955964560476,4.285670240254767,0.8680660810891044,0.18820155,"inflammatory response"),
c("GO:0007049","cell cycle",1.795399173617054,1.373098964247092,0.9916989170777561,0.01365818,"cell cycle"),
c("GO:0007155","cell adhesion",0.7060783158562494,6.913640169325252,0.9923416657212272,0.01208264,"cell adhesion"),
c("GO:0007159","leukocyte cell-cell adhesion",0.006130151056501432,1.6466690285355703,0.9944984690913174,0.00608729,"leukocyte cell-cell adhesion"),
c("GO:0016477","cell migration",0.182534145240741,3.630784142589857,0.9343652520595747,0.01035172,"cell migration"),
c("GO:0030216","keratinocyte differentiation",0.01565965880304327,3.5316526695878427,0.8970531352170081,0.00646486,"keratinocyte differentiation"),
c("GO:0001775","cell activation",0.09565097880185117,1.366349098267021,0.9150812016733332,0.4799785,"keratinocyte differentiation"),
c("GO:0001843","neural tube closure",0.0130685640265839,1.4992242862509764,0.8787062746889309,0.6962195,"keratinocyte differentiation"),
c("GO:0002467","germinal center formation",0.0024414166443147323,1.6466690285355703,0.7835573273448826,0.52511154,"keratinocyte differentiation"),
c("GO:0007368","determination of left/right symmetry",0.025431977741730852,1.6466690285355703,0.8910169430595272,0.44699979,"keratinocyte differentiation"),
c("GO:0007507","heart development",0.08811385330319016,2.5945819179913268,0.8739913543893447,0.64300476,"keratinocyte differentiation"),
c("GO:0008544","epidermis development",0.04074571375048429,1.6977296299576303,0.908299440001792,0.65876759,"keratinocyte differentiation"),
c("GO:0033077","T cell differentiation in thymus",0.010218027154407162,1.8097859415182171,0.8470619610905734,0.51285802,"keratinocyte differentiation"),
c("GO:0060041","retina development in camera-type eye",0.03059088402964931,1.4007638317241864,0.8851353624772086,0.67146902,"keratinocyte differentiation"),
c("GO:0035633","maintenance of blood-brain barrier",0.0018526826578791635,1.8922079909343281,0.9306585373049595,-0,"maintenance of blood-brain barrier"),
c("GO:0097009","energy homeostasis",0.0060736059843579035,1.7426707047181322,0.9267771575935498,0.68469124,"maintenance of blood-brain barrier"),
c("GO:0043616","keratinocyte proliferation",0.0040180263029049,1.572787374784582,0.9946323813381993,0.00593126,"keratinocyte proliferation"),
c("GO:0045944","positive regulation of transcription by RNA polymerase II",0.36927923526581774,6.040481623027001,0.7296321215629707,-0,"positive regulation of transcription by RNA polymerase II"),
c("GO:0007204","positive regulation of cytosolic calcium ion concentration",0.015832620200188184,3.070581074285707,0.8673807321976285,0.55877914,"positive regulation of transcription by RNA polymerase II"),
c("GO:0008360","regulation of cell shape",0.6890449444181898,3.9706162223147903,0.818940424414273,0.22059988,"positive regulation of transcription by RNA polymerase II"),
c("GO:0030335","positive regulation of cell migration",0.0679505458129505,4.749579997691106,0.7613021402025808,0.62837217,"positive regulation of transcription by RNA polymerase II"),
c("GO:0032695","negative regulation of interleukin-12 production",0.0021420603800254598,2.225191826164486,0.8091344659315616,0.66863997,"positive regulation of transcription by RNA polymerase II"),
c("GO:0032731","positive regulation of interleukin-1 beta production",0.006116846333644131,4.489454989793388,0.7590632246157973,0.55427737,"positive regulation of transcription by RNA polymerase II"),
c("GO:0032930","positive regulation of superoxide anion generation",0.0012672748521579196,1.408247169795908,0.8061494119877615,0.53188653,"positive regulation of transcription by RNA polymerase II"),
c("GO:0032956","regulation of actin cytoskeleton organization",0.1640671899148072,1.7426707047181322,0.8436002501115061,0.69763288,"positive regulation of transcription by RNA polymerase II"),
c("GO:0042127","regulation of cell population proliferation",0.23910250064927047,4.328827157284917,0.8479444153292075,0.21740759,"positive regulation of transcription by RNA polymerase II"),
c("GO:0042572","retinol metabolic process",0.012639486714435943,1.728346109625536,0.8634495544268833,0.5497591,"positive regulation of transcription by RNA polymerase II"),
c("GO:0043065","positive regulation of apoptotic process",0.0747093450244594,4.647817481888637,0.7506696825819061,0.63302702,"positive regulation of transcription by RNA polymerase II"),
c("GO:0043087","regulation of GTPase activity",0.17705259942353296,3.036212172654445,0.8450094770213465,0.6427394,"positive regulation of transcription by RNA polymerase II"),
c("GO:0048661","positive regulation of smooth muscle cell proliferation",0.006865236994367313,1.640797474546678,0.781592209283119,0.66114728,"positive regulation of transcription by RNA polymerase II"),
c("GO:0050680","negative regulation of epithelial cell proliferation",0.017881547520212534,2.051564195236017,0.8352881164716099,0.17779103,"positive regulation of transcription by RNA polymerase II"),
c("GO:0050727","regulation of inflammatory response",0.04520612208839445,3.7520267336381936,0.832465522092911,0.17542755,"positive regulation of transcription by RNA polymerase II"),
c("GO:0051056","regulation of small GTPase mediated signal transduction",0.13655634922662307,3.316052869248488,0.8170775930991547,0.60022411,"positive regulation of transcription by RNA polymerase II"),
c("GO:0051092","positive regulation of NF-kappaB transcription factor activity",0.020316311803098615,3.428291168191312,0.8379723857338064,0.29438591,"positive regulation of transcription by RNA polymerase II"),
c("GO:1902017","regulation of cilium assembly",0.009985194504404396,1.7426707047181322,0.8674223622817792,0.17079775,"positive regulation of transcription by RNA polymerase II"),
c("GO:0060271","cilium assembly",0.17316429416848675,17.6458915608526,0.8838008867340083,0,"cilium assembly"),
c("GO:0030036","actin cytoskeleton organization",0.3211826621366746,9.856985199745905,0.9114299003502003,0.53252619,"cilium assembly"));

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

