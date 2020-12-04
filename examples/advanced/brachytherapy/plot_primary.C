{
// Plot the energy spectrum of the primary particles
gROOT -> Reset();

TFile f("primary.root");
  
// Draw histos filled by Geant4 simulation 
//   
  
TCanvas* c1 = new TCanvas("c1", " ");

h10->Draw();
}
