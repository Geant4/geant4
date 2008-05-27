{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);
gStyle->SetNdivisions(210, "x");
gStyle->SetNdivisions(10, "y");
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(0.05, "x");
gStyle->SetLabelSize(0.05, "y");
gStyle->SetTitleOffset(0.9, "x");
gStyle->SetTitleOffset(1.4, "y");
gStyle->SetTitleSize(0.05, "x");
gStyle->SetTitleSize(0.05, "y"); 
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.1);
gStyle->SetTickLength(0.05, "x");
gStyle->SetTickLength(0.03, "y"); 
gStyle->SetPadBorderMode(0);

Int_t nmod = 11;
TString nam[nmod]  = {"BIC","BERT","LHEP","PRECO","FTFP","QGSP","QGSC","RPG","LIEGE","QGSB","FTFC"};
TString fl[nmod]   = {"bic","bert","lhep","prec","ftfp","qgsp","qgsc","rpg","incl","qgsb","ftfc"}; 
Int_t coll[nmod]   = {4, 2, 3, 6, 7, 4, 8, 5, 3, 1, 9};
Int_t active[nmod] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

TString tit[6] = {"E (MeV)","d#sigma/dEd#Omega (mb/srad/MeV)","d#sigma/dE (mb/MeV)","log_{10}(E/MeV)","d#sigma/dEd#Omega (b/srad/MeV)","1/dEd#Omega (1/MeV/srad/#muC)"}; 
TString hist[12]= {"h1","h2","h3","h4","h5","h6","h7","h8","h9","h10","h11","h12"}; 
TH1F* hh[12];
Int_t nplots = 6;

TString txtx = "Data";
TString nfout = "afig.gif";
TString nfout1 = "afig_em.gif";
TString nfout2 = "afig_log.gif";
TString nfout4 = "afig_log0.gif";

TLegend* leg[11];
TLegend* leg1[11];
TGraphErrors* gr[11];
}
