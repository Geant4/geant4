{
  gROOT->Reset();
  
  // Draw histos 1-6 filled by Geant4 simulation with run02.mac
  //   
  TFile f0("testem6_0.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->Divide(2,3);
  for(int i=1; i < 7; ++i)
  {
    c1->cd(i);
    const char* hname = std::string("h"+std::to_string(i)).c_str(); // h1, h2, .. h6
    f0.Get(hname)->Draw("HIST");
  }
}
