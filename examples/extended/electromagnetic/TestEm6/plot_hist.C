{
  gROOT->Reset();

  // Draw histos filled by Geant4 simulation 
  //   

  double minE  = 40;
  double maxE  = 10000; 
  double minY  = 1e-2; //c1
  double maxY  = 1e+2; //c1
  double minY1 = 1e-8; //c2
  double maxY1 = 1;    //c2
 
  //file open
  //
  TFile f("testem6_0.root");

  //CrossSectionPerAtom
  //
  TH1D* h7 = (TH1D*)f.Get("h7");  //ee to MuMu
  TH1D* h8 = (TH1D*)f.Get("h8");  //ee to gammagamma
  TH1D* h9 = (TH1D*)f.Get("h9");  //ee to Hadrons
  TH1D* h10 = (TH1D*)f.Get("h10"); //ee to gammagamma(theory)
  TH1D* h11 = (TH1D*)f.Get("h11"); //ee to MuMu(theory)
  //CrossSectionPerVolume
  //
  TH1D* h12 = (TH1D*)f.Get("h12"); //eBrem
  TH1D* h13 = (TH1D*)f.Get("h13"); //eIoni
  TH1D* h14 = (TH1D*)f.Get("h14"); //ee to MuMu
  TH1D* h15 = (TH1D*)f.Get("h15"); //ee to gammagamma
  TH1D* h16 = (TH1D*)f.Get("h16"); //ee to Hadrons
  //R ratio
  //
  TH1D* h17 = (TH1D*)f.Get("h17"); //R

  //c1
  //
  TCanvas* c1 = new TCanvas("c1", "  "); 
  gPad->SetLogy();
  h7->SetMinimum(minY);
  h7->SetMaximum(maxY);
  h7->GetXaxis()->SetLimits(std::log10(minE),std::log10(maxE));
  h7->SetTitle("Comparison of e+ Annihilation Processes(totcrsPerAtom)");
  h7->GetXaxis()->SetTickLength(0);
  h7->GetXaxis()->SetLabelOffset(999);
  h7->GetYaxis()->SetTitle("microbarn");
  
  h7->SetLineColor(kRed);
  h9->SetLineColor(kBlue);
  h11->SetMarkerColor(kRed);
  h11->SetMarkerStyle(22);
  h10->SetMarkerStyle(22);
   
  c1->cd();
  h7->Draw("HIST");
  h8->Draw("HIST SAME");
  h9->Draw("HIST SAME");
  h10->Draw("HIST P SAME");
  h11->Draw("HIST P SAME");

  TGaxis *xaxis = new TGaxis(std::log10(minE),minY,std::log10(maxE),minY,minE,maxE,510,"G");
  xaxis->SetTitle("GeV");
  xaxis->Draw();

  gStyle->SetOptStat(0);

  TLegend* leg = new TLegend(0.78,0.59,0.98,0.76);
  leg->AddEntry(h8,"to 2 gammas","l");
  leg->AddEntry(h7,"to MuPair","l");
  leg->AddEntry(h9,"to Hadrons","l");
  leg->Draw();

  //c2
  //
  TCanvas* c2 = new TCanvas("c2", "  ");
  c2->cd();

  gPad->SetLogy();
  h15->SetMinimum(minY1);
  h15->SetMaximum(maxY1);
  h15->GetXaxis()->SetLimits(std::log10(minE),std::log10(maxE));
  h15->SetTitle("Comparison of EM Processes(totcrsPerVolume)");
  h15->GetXaxis()->SetTickLength(0);
  h15->GetXaxis()->SetLabelOffset(999);
  h15->GetYaxis()->SetTitle("1/mm");

  h12->SetLineColor(6);
  h13->SetLineColor(kGreen);
  h14->SetLineColor(kRed);
  h16->SetLineColor(kBlue);

  c2->cd();
  h15->Draw("HIST");
  h13->Draw("HIST SAME");
  h14->Draw("HIST SAME");
  h12->Draw("HIST SAME");
  h16->Draw("HIST SAME");

  TGaxis *xaxis1 = new TGaxis(std::log10(minE),minY1,std::log10(maxE),minY1,minE,maxE,510,"G");
  xaxis1->SetTitle("GeV");
  xaxis1->Draw();

  gStyle->SetOptStat(0);

  TLegend* leg1 = new TLegend(0.78,0.59,0.98,0.76);
  leg1->AddEntry(h15,"to 2 gammas","l");
  leg1->AddEntry(h14,"to MuPair","l");
  leg1->AddEntry(h16,"to Hadrons","l");
  leg1->AddEntry(h12,"Bremsstrahlung","l");
  leg1->AddEntry(h13,"Ionization","l");
  leg1->Draw(); 
}  
