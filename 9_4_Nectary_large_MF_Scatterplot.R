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
revigo.data <- rbind(c("GO:0004029","aldehyde dehydrogenase (NAD+) activity",0.15095807017417034,0.8578801568724232,6.105935614273045,4.83300040911894,-2.5421181032660076,0.8317181449581453,-0),
c("GO:0004030","aldehyde dehydrogenase [NAD(P)+] activity",0.16055538094028188,0.881916182013904,6.411291278991541,4.859768557874373,-2.5421181032660076,0.8311552866858052,0.80608029),
c("GO:0004386","helicase activity",1.24556595594382,0.13464166394064303,0.4352187217382287,5.749505197054165,-1.5606673061697374,0.7891487328548674,0.03957792),
c("GO:0004527","exonuclease activity",0.7727986135167833,6.28920879893663,0.8167603338012376,5.542205274791192,-1.5409092103994135,0.8114910188965833,0.78646589),
c("GO:0004532","RNA exonuclease activity",0.24062019205886398,5.569533755820939,0.998993025405868,5.035473765808079,-1.6581699430794896,0.7839695685697119,0.83610407),
c("GO:0004540","RNA nuclease activity",0.9410486960186195,5.38880381609263,0.7658196891471721,5.627750814997703,-1.5972229303896526,0.7827708118997793,0.73625869),
c("GO:0004645","1,4-alpha-oligoglucan phosphorylase activity",0.04346952008966833,2.228374342725698,-5.393881628715051,4.292344693840592,-1.3159629625134803,0.825009919414179,0.59537292),
c("GO:0005346","purine ribonucleotide transmembrane transporter activity",0.06273731889898977,-6.068359811269219,-1.8427730116015903,4.451678999573242,-2.684029654543082,0.5634640775092842,0),
c("GO:0008184","glycogen phosphorylase activity",0.026534524174512637,2.695637271943583,-5.4315399131068895,4.077985291029502,-1.4277093938485823,0.8290182917150358,0.57206112),
c("GO:0008379","thioredoxin peroxidase activity",0.08722731430123927,-0.5704386224338298,6.384788032140077,4.594801235778386,-1.3417023496918101,0.8561195807547959,0.79255133),
c("GO:0010294","abscisic acid glucosyltransferase activity",0.0002683166826939687,1.488582397725432,-5.815568358782729,2.0863598306747484,-1.3417023496918101,0.849319830088768,0.46611246),
c("GO:0015095","magnesium ion transmembrane transporter activity",0.11206545614301731,-6.327325256407695,0.4735160114502026,4.703618051075181,-1.3417023496918101,0.7166705477085538,0.41029382),
c("GO:0015444","P-type magnesium transporter activity",0.004432768997563995,-5.263518805685632,0.5528653706127323,3.3010299956639813,-1.3417023496918101,0.6967830189508664,0.79425954),
c("GO:0015605","organophosphate ester transmembrane transporter activity",0.11816799755437105,-6.217864707242028,-1.5696462709340813,4.726645720240912,-1.3159629625134803,0.6046532990002592,0.78947045),
c("GO:0015662","P-type ion transporter activity",0.2160503668997799,-4.9508803928702685,-0.10268568532183235,4.988697160029578,-1.3417023496918101,0.6570298004770673,0.48305443),
c("GO:0016620","oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",0.5304576466994853,0.5847831923348995,6.284820976077479,5.378790224142972,-1.6283780728239787,0.8223950293383303,0.84272965),
c("GO:0016717","oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water",0.1033795351007671,1.8206249186137575,6.120525568944156,4.668581584710601,-1.4277093938485823,0.8880020354532742,0.29025415),
c("GO:0016758","hexosyltransferase activity",1.2623656847710059,2.529030430600198,-5.7993345556826235,5.755323637746441,-1.4155556928348239,0.8165364328084201,0.04473212),
c("GO:0016796","exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters",0.4035194633595386,6.24120607815936,1.1611204319171735,5.260004568309291,-1.6581699430794896,0.8003071019738288,0.03262267),
c("GO:0016896","RNA exonuclease activity, producing 5'-phosphomonoesters",0.23929191360486657,5.473125431865983,1.2136717152040317,5.033069741661653,-1.6581699430794896,0.7840484474410742,0.83570328),
c("GO:0019829","ATPase-coupled monoatomic cation transmembrane transporter activity",0.4582427616696358,-5.273588228542782,-0.16763277681561775,5.315235409617727,-1.3417023496918101,0.6393097040347915,0.71589944),
c("GO:0022804","active transmembrane transporter activity",3.418048523453295,-6.331668824549807,-0.688567182215685,6.187916199908961,-1.4586704223333065,0.7004983616732254,0.38815234),
c("GO:0022853","active monoatomic ion transmembrane transporter activity",1.2952488921071121,-6.131722039402982,-0.16552577035806235,5.766491667382492,-1.3128278954052004,0.6610377869596482,0.62204132),
c("GO:0033946","xyloglucan-specific endo-beta-1,4-glucanase activity",8.426474332537859E-05,-3.2984619708163225,-6.6043243758096235,1.591064607026499,-1.3417023496918101,0.9510358561732292,0.11987935),
c("GO:0047807","cytokinin 7-beta-glucosyltransferase activity",4.2132371662689295E-05,1.2551211328778835,-6.201858329915125,1.3010299956639813,-1.3417023496918101,0.8591851440074177,0.3785305),
c("GO:0051920","peroxiredoxin activity",0.17148318765363627,-0.3486637008187473,6.487929876701791,4.888364858208229,-1.3417023496918101,0.851206332093804,0.30167262),
c("GO:0052610","beta-cryptoxanthin hydroxylase activity",4.4349864908094E-06,-2.763841994439781,6.468801585009136,0.47712125471966244,-1.3417023496918101,0.931528219857981,0.16565557),
c("GO:0097363","protein O-acetylglucosaminyltransferase activity",0.017662333699648435,2.0865375324026987,-5.940735797560582,3.901240302073309,-1.3417023496918101,0.8272215884973738,0.55416948),
c("GO:0102659","UDP-glucose: 4-methylthiobutylhydroximate S-glucosyltransferase activity",8.8699729816188E-06,3.4468510641751915,-5.13804138298079,0.6989700043360189,-1.3417023496918101,0.8769393948077832,0.34991333),
c("GO:0140358","P-type transmembrane transporter activity",0.22194224645282018,-4.939036267646892,-0.5447868815777004,5.00038201108384,-1.3417023496918101,0.6678581010501607,0.6719511));

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
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/6,max(one.data$plot_X)+one.x_range/15);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/12,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

nectary_MF <- p1
# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

