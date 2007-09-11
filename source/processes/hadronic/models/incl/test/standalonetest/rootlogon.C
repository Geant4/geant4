{
// TODO:
// teksti snadimmaksi
// tekstin boldaus? pois
// tekstin fontti helveticaksi, palatinoksi tms
// stats-boksiin under/overflow
// fontti/boldaus stats-boksissa, otsikossa, akseleilla
// grid-mahdollisuus?
// 
// --------- tiedostoihin/skripteihin
// hyv/hyl-viivat skripteihin
// pvm normaaliin muotoon 

cout <<"Applying customized plotting style settings..." << endl;
gROOT->SetStyle("Plain");

gStyle->SetPadBorderMode(0);  // is this needed
//gStyle->SetPadFillColor(0);
gStyle->SetPadBorderSize(1);

gStyle->SetFrameBorderMode(0);
gStyle->SetFrameBorderMode(0);


// Set Canvas default size - small/large
//
//gStyle->SetCanvasDefW(350);
//gStyle->SetCanvasDefH(250);
gStyle->SetCanvasDefW(900); // Old: 500
gStyle->SetCanvasDefH(670); // Old: 500
//gStyle->SetCanvasDefW(1400);
//gStyle->SetCanvasDefH(1000);

// COLORS
//
// Change background beige to white
gStyle->SetFillColor(0);
//gStyle->SetFillStyle(3001);
gStyle->SetCanvasColor(10);
gStyle->SetTitleFillColor(0);
gStyle->SetStatColor(10);

// FONTS: Text size for independent texts
gStyle->SetTextSize(0.04);
//gStyle->SetTextFont(42);  // default 62
gStyle->SetTextFont(102);  // default 62

//gStyle->SetTitleFont(42,"pad");  // default 62
gStyle->SetTitleFont(102,"pad");  // default 62
gStyle->SetTitleFont(102,"X");  // default 62
gStyle->SetTitleFont(102,"Y");  // default 62
gStyle->SetLabelFont(102,"X");  // default 62
gStyle->SetLabelFont(102,"Y");  // default 62
gStyle->SetStatFont(102);  // default 62



// Histograms
//gStyle->SetHistFillStyle(3001); // 3001 solid, opaque
//gStyle->SetHistFillColor(1);  // black
// remember +100 for darker, +150 for lighter version
//gStyle->SetHistLineColor(38); 
// fill styles: transparent 4000 - 4100 opaque
// fill patterns: page 140
//
// usual light blue
//gStyle->SetHistFillColor(38+150); // 18 light grey, 38 blue
//
// nice b/w pattern
gStyle->SetHistFillColor(kWhite);  // black
//gStyle->SetHistFillStyle(3017); // 3001 solid, opaque
//
// grey-looking b/w pattern
//gStyle->SetHistFillColor(3);  // black
//gStyle->SetHistFillStyle(3002); 

// STATS Box position and options
//
//gStyle->SetStatX(0.882653);  // nicely almost in the corner
//gStyle->SetStatY(0.882263);
gStyle->SetStatX(0.89988); // fX2NDC
gStyle->SetStatY(0.90);// fY2NDC
gStyle->SetOptStat(0);

gStyle->SetMarkerStyle(2);

//If the statistics box is drawn, you can select the type of 
//information displayed with gStyle->SetOptStat(mode). 
//The mode has up to seven digits that can be set 
//to on (1) or off (0). 
//       mode = iourmen ( default = 0001111 ) = 1 
//n: the name of histogram is printed = 1 
//e: the number of entries printed = 1 
//m: the mean value printed = 1 
//r: the root mean square printed = 1 
//u: the number of underflows printed = 1 
//o: the number of overflows printed = 1 
//i: the integral of bins printed
//
//WARNING: never call SetOptStat(000111); but SetOptStat(1111) ,
//0001111 will be taken as an octal number. With the option "same",
//the statistic box is not redrawn. With the option "sames" , the
//statistic box is drawn. If it hides the previous statistics box, you
//can change its position with these lines (if h is the pointer to the
//histogram):



// Axis definitions
gStyle->SetTickLength(0.015,"X");
gStyle->SetTickLength(0.01,"Y");

// Main title of graph
//
// these worked almost ok
//gStyle->SetTitleX(0.0997537);
//gStyle->SetTitleY(0.94645);  // fY2NDC, top
//
// new try
//
gStyle->SetTitleX(0.10);  //fX1NDC
//gStyle->SetTitleY(0.9486);  // fY2NDC, top
gStyle->SetTitleY(0.9786);  // fY2NDC, top
//
// some good values
//.099182 fX1NDL
//.90015 fY1NDC
//.1809 fX2NDc
//.9464 fY2NDC
//
//SetTitlePS
//SetTitleColor
//SetTitleFont
//SetTitleOffset
//SetTitleSize
//SetTitleFillColor
//SetTitleTextColor
//SetTitleStyle
//SetTitleFontSize
gStyle->SetTitleBorderSize(0.0);
//SetTitleXOffset
//SetTitleXSize
//SetTitleYOffset
//SetTitleYSize
//SetTitleX
//SetTitleY
//SetTitleW
//SetTitleH
//SetTitle
gStyle->SetLegendBorderSize(0.0);

gStyle->SetLabelSize(0.035,"X"); 
gStyle->SetLabelSize(0.035,"y");  // fY2NDC
gStyle->SetTitleSize(0.035,"X"); 
gStyle->SetTitleSize(0.035,"y"); 

gStyle->SetOptTitle(kFALSE);
cout <<"Done setting up plotting style." << endl;
cout << endl;
}


