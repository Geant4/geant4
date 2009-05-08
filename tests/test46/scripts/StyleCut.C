{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);
gStyle->SetNdivisions(210, "x");
gStyle->SetNdivisions(10, "y");
gStyle->SetTextFont(1);
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(0.03, "x");
gStyle->SetLabelSize(0.03, "y");
gStyle->SetTitleOffset(1.4, "x");
gStyle->SetTitleOffset(1.6, "y");
gStyle->SetTitleSize(0.05, "x");
gStyle->SetTitleSize(0.05, "y");
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadLeftMargin(0.18);
gStyle->SetPadRightMargin(0.1);
gStyle->SetTickLength(0.05, "x");
gStyle->SetTickLength(0.03, "y");
gStyle->SetPadBorderMode(0);

gStyle->SetMarkerSize(1.5);
TCanvas c1("c1"," ",0.5, 5, 800, 800);

Int_t col[6] = {4, 2, 5, 5, 6, 7};
Int_t mar[6] = {21, 20, 22, 25, 24, 26};

TString tit[3]  = {"#gamma 10 GeV in HCAL","#pi^{-} 30 GeV in HCAL","#pi^{-} 30 GeV in ECAL + HCAL"};
TString axtit[6] = {"cut (mm)","E_{vis}/E_{0}","#sigma_{E}/E_{vis} (%)","E_{vis} (GeV)","E_{HCAL} (GeV)","E_{rec}/E_{0}"};
}
