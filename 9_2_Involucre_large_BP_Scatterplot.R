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
revigo.data <- rbind(c("GO:0000902","cell morphogenesis",0.6354075502089755,-3.9895282255618345,-4.7859661030155305,5.411282913017384,-1.914709421769935,0.6909654934597287,0.69769973),
c("GO:0001558","regulation of cell growth",0.15345225803226778,0.4409156994036262,5.217685763256858,4.794209116346496,-2.0877779434675845,0.8920330771530458,0.72249856),
c("GO:0002376","immune system process",0.9427113541805001,1.0174177134817612,-8.954783238374535,5.582608730687733,-1.4461169733561257,1,-0),
c("GO:0006596","polyamine biosynthetic process",0.17095684269545167,6.8242199615677395,2.4310911349550377,4.841121607196544,-1.530767257493388,0.9882057554623386,0.10821682),
c("GO:0006950","response to stress",6.919654498638822,4.985885502359611,-3.7443355515089802,6.448313422563495,-1.4038429190838275,0.8417918538430418,0.50338616),
c("GO:0006955","immune response",0.6433982378290884,5.846627247283121,-3.869078461381926,5.416710388184187,-1.4461169733561257,0.8723509626235395,0.29788224),
c("GO:0006970","response to osmotic stress",0.1720610462900879,4.219098111819338,-4.426781686533095,4.843917638006392,-1.853871964321762,0.8432287663361118,0.65932678),
c("GO:0009607","response to biotic stimulus",1.0265223788055218,5.223818705827184,-4.684789475864123,5.619598321840807,-2.0599818449923366,0.8670312966150227,0.24307207),
c("GO:0009627","systemic acquired resistance",0.011522660278848872,5.882119447315008,-2.331521736850085,3.6698745024898023,-2.0877779434675845,0.8658670185354438,0.30530823),
c("GO:0009651","response to salt stress",0.07343446852364134,4.22522428267175,-4.89256583748561,4.474143389760806,-2.1493537648169334,0.8513762693826683,-0),
c("GO:0009653","anatomical structure morphogenesis",1.6026307453462993,-4.337130231455072,-5.037305612661848,5.813062995420341,-2.55129368009492,0.7214650151429309,-0),
c("GO:0009826","unidimensional cell growth",0.016627137163874758,-4.720667495547815,-3.891725489996526,3.8291107101552946,-1.6802695056697754,0.6847087640159675,0.83635703),
c("GO:0010043","response to zinc ion",0.10276734124633234,4.3543028455457335,-2.0776554186463874,4.620094394032491,-1.9200955323332793,0.8899664265616133,0.2502007),
c("GO:0010053","root epidermal cell differentiation",0.006001642305667808,-3.824883464274691,-5.487830786015883,3.3866772839608377,-1.4769041617474323,0.6250845768633835,0.77151248),
c("GO:0010184","cytokinin transport",7.394220499796067E-05,5.241260441145661,3.8430996814521805,1.4913616938342726,-1.530767257493388,0.9591525656875045,0.10211233),
c("GO:0010191","mucilage metabolic process",0.001777077660117655,-2.737016918186832,1.6506281952423472,2.858537197569639,-1.530767257493388,0.9929835360598858,0.03809888),
c("GO:0010192","mucilage biosynthetic process",0.0017228533764524837,-5.6979589942728275,2.963175563382194,2.845098040014257,-1.530767257493388,0.9894730261738934,0.00443713),
c("GO:0016049","cell growth",0.07489112996210118,-6.392035300028961,-1.3209683087644843,4.48267353350333,-1.4945786724167192,0.8989513970803041,0.89877773),
c("GO:0022603","regulation of anatomical structure morphogenesis",0.9381565143526257,1.9536111021376046,4.603300583374904,5.580505296928462,-1.8544928285903375,0.8730282749840231,0.54797281),
c("GO:0022604","regulation of cell morphogenesis",0.6684621805832305,2.026077228928679,4.883550085443135,5.433307300026958,-1.8544928285903375,0.8676458517563548,0.85389249),
c("GO:0042221","response to chemical",4.8560360057130705,4.553828080918956,-3.4370524385368877,6.2945109760284925,-2.048176964684088,0.8467597543594402,0.37706273),
c("GO:0048653","anther development",0.006450225015988769,-3.57969945775161,-5.783067113704152,3.417969642214737,-1.530767257493388,0.6106187174509824,0.77458391),
c("GO:0048655","anther wall tapetum morphogenesis",0.00025140349699306625,-3.354201968736178,-5.631116326904885,2.012837224705172,-1.530767257493388,0.6022264868862206,0.84787201),
c("GO:0050896","response to stimulus",17.567785530535815,-7.018323890084432,1.4665883142466385,6.852945939203233,-2.2276782932770804,1,-0),
c("GO:0051510","regulation of unidimensional cell growth",0.00424921204721614,1.2188505050313243,4.814572404768544,3.2367890994092927,-2.3010299956639813,0.8512457456763143,0.73752684),
c("GO:0051513","regulation of monopolar cell growth",0.0030119124835835983,1.2973825215889612,5.104202942924347,3.0874264570362855,-3.0604807473813813,0.8535582400031556,0),
c("GO:0055085","transmembrane transport",12.81017399389553,-2.974369668965129,5.051720975067988,6.715783969524812,-1.4610494379856385,0.9754080840644579,0.18512772),
c("GO:0060560","developmental growth involved in morphogenesis",0.07016868780289807,-4.603014019299351,-4.067642889477935,4.454387467146955,-1.6802695056697754,0.6642247056419217,0.56279732),
c("GO:0080170","hydrogen peroxide transmembrane transport",0.00028098037899225057,-4.214886760130968,4.7752065335252105,2.060697840353612,-1.9200955323332793,0.9812347858027399,-0),
c("GO:0098542","defense response to other organism",0.8128269416212488,5.333581259223512,-3.0375565186846982,5.518228264418013,-1.4461169733561257,0.8254773317627071,0.66056519),
c("GO:1905392","plant organ morphogenesis",0.021231271795081108,-3.646968387867472,-5.193480162544161,3.9352552817840474,-1.7602001815529014,0.6312913887226849,0.50936355));

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
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/15,max(one.data$plot_X)+one.x_range/12);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/18,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

involucre_BP <- p1

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

