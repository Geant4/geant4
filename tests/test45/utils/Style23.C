{ 
gROOT->ProcessLine(".x   $TEST45/utils/StyleBase.C");
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(.055, "x");
gStyle->SetLabelSize(.055, "y");
gStyle->SetTitleOffset(.95, "x");
gStyle->SetTitleOffset(.95, "y");
gStyle->SetTitleSize(.07, "x");
gStyle->SetTitleSize(.07, "y");
gStyle->SetPadBottomMargin(.15);          
gStyle->SetPadTopMargin(.05);
gStyle->SetPadLeftMargin(.15);
gStyle->SetPadRightMargin(.03);

gStyle->SetMarkerStyle(21); 
gStyle->SetMarkerColor(1); 
gStyle->SetMarkerSize(.6);
TCanvas c1("c1"," ",0.5, 5, 800, 600);
}
