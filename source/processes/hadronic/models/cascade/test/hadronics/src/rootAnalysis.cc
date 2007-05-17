//#ifdef  G4ANALYSIS_USE_ROOT

#include "rootAnalysis.hh"

rootAnalysis* rootAnalysis::instance = 0;

rootAnalysis::rootAnalysis(): hfile(0),ntuple(0),h1(0) {  
  G4cout << "::: rootAnalysis::rootAnalysis()" << G4endl;
}

rootAnalysis::~rootAnalysis() { 
  G4cout << "::: rootAnalysis::~rootAnalysis()" << G4endl;
  //  delete hfile;
  //ntuple = 0;

  delete ntuple;
  ntuple = 0;

  delete h1;
  h1 = 0;
}

rootAnalysis* rootAnalysis::getInstance() {
  if (instance == 0) instance = new rootAnalysis;
  return instance;
}

void rootAnalysis::book() {
  G4cout << "::: rootAnalysis::book()" << G4endl;
  hfile  = new TFile("hadronics.root","RECREATE","Demo");
  h1 = new TH1F("h1" , "Energy deposited / step", 100, 0 , 1);
  ntuple = new TNtuple("ntuple" , "Demo ntuple","iSlice:jSlice:kSlice:fEnergy");
}

void rootAnalysis::FillEnergyDeposit(G4int i, G4int j, G4int k, G4double energy) {
  if (ntuple)
    {
      //     ntuple->Fill(i, j, k, energy);
    }
}

void rootAnalysis::BraggPeak(G4int slice, G4double energy) {
  //  h1->Fill(slice,energy);
}

void rootAnalysis::SecondaryProtonEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryGammaEnergyDeposit(G4int slice, G4double energy) {
  G4cout << "::: rootAnalysis::SecondaryGammaEnergyDeposit" << G4endl;
  //  ntuple->Fill(slice, 1, 1, energy);
}

void rootAnalysis::SecondaryElectronEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryTritonEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::SecondaryPionEnergyDeposit(G4int slice, G4double energy) {
}

void rootAnalysis::electronEnergyDistribution(G4double energy) {
  G4cout << "::: rootAnalysis::electronEnergyDistribution(G4double energy)" << G4endl;  
  h1->Fill(energy);
}

void rootAnalysis::gammaEnergyDistribution(G4double energy) {
  G4cout << "::: rootAnalysis::gammaEnergyDistribution(G4double energy = " << energy/MeV << " MeV)" << G4endl;  
  h1->Fill(energy);
}

void rootAnalysis::deuteronEnergyDistribution(G4double energy) {
}

void rootAnalysis::tritonEnergyDistribution(G4double energy) {
}

void rootAnalysis::alphaEnergyDistribution(G4double energy) {
}

void rootAnalysis::genericIonInformation(G4int a, G4double z, G4int electronOccupancy, G4double energy)  {
  // if (ionTuple)
}

void rootAnalysis::finish() {  
  G4cout << "::: rootAnalysis::finish()" << G4endl;
  hfile->Print();
  hfile->Write();
}
//#endif











