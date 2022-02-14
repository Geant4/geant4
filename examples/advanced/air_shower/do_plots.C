{
// Plot the energy spectrum of the primary particles
gROOT -> Reset();

TFile f("ultra.root");
  
TCanvas *c1 = new TCanvas("Photons energy");
PhotonEnergy->Draw();

TCanvas *c2 = new TCanvas("Number of photons");
NumberDetectedPhotons->Draw();
}
