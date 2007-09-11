{
  // Compare the results of the INCL demo:

  TFile *protonFile = new TFile("Pb208Proton1GeV.root");
  TFile *alphaFile = new TFile("Pb208Alpha1GeV.root");

  TFile *protonFileF77 = new TFile("../data/pb208proton1gev_incl42_noabla.root");
  TFile *alphaFileF77 = new TFile("../data/pb208alpha1gev_incl42_noabla.root");

  TTree *proton = (TTree*) protonFile->Get("h101");
  TTree *alpha = (TTree*) protonFile->Get("h101");

  TTree *protonF77 = (TTree*) protonFileF77->Get("h101");
  TTree *alphaF77 = (TTree*) alphaFileF77->Get("h101");
  
  TText *xlabel = new TText();
  xlabel-> SetNDC();
  xlabel -> SetTextColor(1);
  xlabel -> SetTextSize(0.04);
  xlabel -> SetTextAlign(22);
  xlabel -> SetTextAngle(0);

  TCanvas *c1 = new TCanvas();

  protonF77->SetLineColor(kRed);
  alphaF77->SetLineColor(kRed);

  c1->SetLogy();
  protonF77->Draw("Avv >> fhisto");
  proton->Draw("Avv >> histo", "", "same");
  histo->SetTitle("C++ INCL 4.2");
  fhisto->SetTitle("FORTRAN INCL 4.2");
  c1->BuildLegend(0.2, 0.7, 0.7, 0.85);
  xlabel -> DrawText(0.5, 0.95, "p(1 GeV) + 208Pb");
  c1->SaveAs("Pb208Proton1GeV.eps");


  TCanvas *c2 = new TCanvas();
  c2->SetLogy();
  alphaF77->Draw("Avv >> fhisto");
  alpha->Draw("Avv >> histo", "", "same");
  histo->SetTitle("C++ INCL 4.2");
  fhisto->SetTitle("FORTRAN INCL 4.2");
  c2->BuildLegend(0.2, 0.7, 0.7, 0.85);
  xlabel -> DrawText(0.5, 0.95, "Alpha(1 GeV) + 208Pb");
  c2->SaveAs("Pb208Alpha1GeV.eps");
}

