//#ifdef G4ANALYSIS_USE_ROOT
#ifndef ROOTANALYSIS_HH
#define ROOTANALYSIS_HH 1

#include "globals.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"

class TFile;
class TNtuple;
class TH1F;

class rootAnalysis
{
  private:
  rootAnalysis(); // Notice private constructor

public:
  ~rootAnalysis();
  
  static rootAnalysis* getInstance();
  
  void book();
  // Book the ntuples and histograms in a .hbk file
  
  void FillEnergyDeposit(G4int voxelXId, G4int voxelYId, G4int voxelZId, 
                         G4double energyDeposit);
  // Fill the ntuple with the energy deposit in the phantom 

  void BraggPeak(G4int, G4double);
  // Fill 1D histogram with the Bragg peak in the phantom

  void SecondaryProtonEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary protons

   void SecondaryNeutronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary neutrons

  void SecondaryAlphaEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary alpha particles

  void SecondaryGammaEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary gamma

  void SecondaryElectronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary electrons

  void SecondaryTritonEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary tritons

  void SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary deuterons

  void SecondaryPionEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary pions

  void electronEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary electrons originated in the phantom

  void gammaEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary gamma originated in the phantom

  void deuteronEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary deuterons originated in the phantom

  void tritonEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary tritons originated in the phantom

  void alphaEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary alpha originated in the phantom

  void genericIonInformation(G4int, G4double, G4int, G4double);
 
  void finish();
  // Close the .hbk file with the histograms and the ntuples

private:
  static rootAnalysis* instance;

  TFile *hfile;
  TNtuple *ntuple;
  TH1F *h1;
};
#endif
//#endif



