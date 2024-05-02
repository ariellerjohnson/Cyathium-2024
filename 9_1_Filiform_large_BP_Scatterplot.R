# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0001558","regulation of cell growth",0.15345225803226778,-4.045844930869528,3.1902410819222715,4.794209116346496,-2.3010299956639813,0.8412555923460938,-0),
c("GO:0003002","regionalization",0.20289494577423747,6.516083089391241,2.058328990177226,4.915505361754376,-1.7144426909922261,0.6618682483107127,0.64094789),
c("GO:0006346","DNA methylation-dependent heterochromatin formation",0.018290836776328874,-0.8978343724799159,7.341814266641133,3.8705209500127644,-1.8153085691824011,0.7053351625177444,0.82087165),
c("GO:0006558","L-phenylalanine metabolic process",0.12549471032253884,3.6855899324508146,-5.626712207276226,4.7068628073573535,-1.6382721639824072,0.9430019803429512,0.47756098),
c("GO:0006886","intracellular protein transport",1.5838297073554848,-6.132514850488149,-0.13135940515688044,5.807938017848966,-1.4111682744057927,0.9943744110001392,0.01163917),
c("GO:0007389","pattern specification process",0.24841623191114867,6.465316458307156,1.5602786626401322,5.003413136276344,-1.4365189146055892,0.6918008904221694,0.55885723),
c("GO:0009072","aromatic amino acid metabolic process",0.7377977862098182,3.3005417965542034,-5.651863744862038,5.4761676559933745,-1.518557371497695,0.944820420382996,0.00870532),
c("GO:0009266","response to temperature stimulus",0.4441067421784182,-3.423571097911933,-4.133398420223469,5.255718634099335,-2.431798275933005,0.9621747570822731,0),
c("GO:0009409","response to cold",0.07047431558355632,-3.6162502206863127,-3.947292473260246,4.4562749128401675,-1.5114492834995557,0.9621747570822731,0.7065124),
c("GO:0009555","pollen development",0.012900450031977538,6.2201951007526395,0.8960620289877729,3.718916686014861,-1.5934598195660448,0.7127498410876127,0.44433299),
c("GO:0009653","anatomical structure morphogenesis",1.6026307453462993,5.908903559698763,1.5434822397237393,5.813062995420341,-1.764471553092451,0.6805641557273501,0.66725502),
c("GO:0010093","specification of floral organ identity",0.0003844994659893955,6.222157415323006,3.079491415528795,2.1958996524092336,-2.167491087293764,0.5882854937832545,0.4532593),
c("GO:0010152","pollen maturation",0.0005175954349857247,6.266691089386378,0.4941985808381901,2.3242824552976926,-1.6382721639824072,0.7478123142326125,0.82012329),
c("GO:0016049","cell growth",0.07489112996210118,-1.643316132194734,-1.6037835324347738,4.48267353350333,-1.7569619513137056,0.9954851199558575,-0),
c("GO:0031047","regulatory ncRNA-mediated gene silencing",0.26322685557224024,-2.3329761458690124,6.247973996674981,5.02856311976092,-1.7423214251308154,0.7800657793556054,0.16576864),
c("GO:0031048","regulatory ncRNA-mediated heterochromatin formation",0.03419087559105702,-1.1514613357356696,7.201417818697978,4.142170386276096,-1.8153085691824011,0.7025518873587674,0.63785474),
c("GO:0040007","growth",0.14923755234738403,-5.137552874254319,-1.7829742579591252,4.782114147479071,-1.728158393463501,1,-0),
c("GO:0040008","regulation of growth",0.2209540969749061,-4.627629700994294,4.986950618016621,4.952535760655082,-2.0268721464003012,0.9147275716533091,0.15250965),
c("GO:0040029","epigenetic regulation of gene expression",0.14583128143714463,-1.3000449008832542,6.913791570419948,4.772086889479079,-1.489454989793388,0.7562624604222571,0.66347235),
c("GO:0048448","stamen morphogenesis",0.0003968231668223889,5.844099398725171,3.4806968603354846,2.2095150145426308,-1.7144426909922261,0.6353290532104999,0.86560637),
c("GO:0048449","floral organ formation",0.0007640694516455937,5.958922939452267,3.1251349319692556,2.4927603890268375,-1.7144426909922261,0.5998026474970356,0.86448028),
c("GO:0048646","anatomical structure formation involved in morphogenesis",0.4177389518761454,5.8357472089501785,2.0493002407992416,5.229136392540384,-2.3010299956639813,0.6631122531692216,-0),
c("GO:0051128","regulation of cellular component organization",1.6291932501217334,-3.236992075806403,4.784225677744801,5.820202116511591,-1.5406075122407692,0.8963014588373961,0.20470873),
c("GO:0051510","regulation of unidimensional cell growth",0.00424921204721614,-4.55215412230862,2.397154801300757,3.2367890994092927,-1.4365189146055892,0.866263993685493,0.73752684),
c("GO:0055062","phosphate ion homeostasis",0.04733779963969443,-0.46346297576871043,-5.6875566906851365,4.283459536376994,-2.167491087293764,0.9726201306855664,0.51004304),
c("GO:0055081","monoatomic anion homeostasis",0.027077635470253197,-0.358792605466418,-5.353278149823162,4.0408791245157865,-2.167491087293764,0.9726201306855664,-0),
c("GO:0070828","heterochromatin organization",0.12110500808582658,0.5141965933880268,8.271983449129392,4.691399799098994,-1.489454989793388,0.8081848151009099,0.00906521),
c("GO:0090698","post-embryonic plant morphogenesis",0.015421879222407998,5.236636034091634,1.81670763939099,3.7964355588101744,-1.518557371497695,0.7145821332892383,0.5706376),
c("GO:0090701","specification of plant organ identity",0.00039435842665579027,6.493169740692264,2.779201786883495,2.2068258760318495,-2.167491087293764,0.6450367662688867,0.82062815),
c("GO:1902221","erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process",0.12549471032253884,4.028860970408261,-5.572620918113012,4.7068628073573535,-1.6382721639824072,0.9430019803429512,0.52275515),
c("GO:1905393","plant organ formation",0.0036354917457330667,5.348604991700774,2.6940866165781774,3.1690863574870227,-1.7144426909922261,0.7052159421904909,0.51811641));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = dispensability), alpha = I(0.6) );
p1 <- p1 + scale_colour_gradientn( colours = c("red", "gray", "blue"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = dispensability), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.2, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/12,max(one.data$plot_X)+one.x_range/8);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/16,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

filiform_BP <- p1

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

