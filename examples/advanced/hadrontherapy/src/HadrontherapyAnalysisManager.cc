//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: HadrontherapyAnalisysManager.cc;
// See more at: http://geant4infn.wikispaces.com/HadrontherapyExample

// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
//
// (a) Laboratori Nazionali del Sud
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
//
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyAnalysisFileMessenger.hh"
#include <time.h>
#ifdef ANALYSIS_USE
HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::instance = 0;

#ifdef G4ROOTANALYSIS_USE
#undef G4ANALYSIS_USE
#endif

/////////////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_USE
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() :
  AnalysisFileName("DoseDistribution.root"), aFact(0), theTree(0), histFact(0), tupFact(0), h1(0), h2(0), h3(0),
  h4(0), h5(0), h6(0), h7(0), h8(0), h9(0), h10(0), h11(0), h12(0), h13(0), h14(0), ntuple(0),
  ionTuple(0)
{
	fMess = new HadrontherapyAnalysisFileMessenger(this);
}
#endif
#ifdef G4ROOTANALYSIS_USE
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() :
  AnalysisFileName("DoseDistribution.root"),theTFile(0), th1(0), th2(0), th3(0),
  th4(0), th5(0), th6(0), th7(0), th8(0), th9(0), th10(0), th11(0), th12(0), th13(0), th14(0), theROOTNtuple(0),
  theROOTIonTuple(0)
{
  fMess = new HadrontherapyAnalysisFileMessenger(this);
}
#endif
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager()
{
delete(fMess); //kill the messenger
#ifdef G4ANALYSIS_USE
  delete ionTuple;
  ionTuple = 0;

  delete ntuple;
  ntuple = 0;

  delete h14;
  h14 = 0;

  delete h13;
  h13 = 0;

  delete h12;
  h12 = 0;

  delete h11;
  h11 = 0;

  delete h10;
  h10 = 0;

  delete h9;
  h9 = 0;

  delete h8;
  h8 = 0;

  delete h7;
  h7 = 0;

  delete h6;
  h6 = 0;

  delete h5;
  h5 = 0;

  delete h4;
  h4 = 0;

  delete h3;
  h3 = 0;

  delete h2;
  h2 = 0;

  delete h1;
  h1 = 0;

  delete tupFact;
  tupFact = 0;

  delete histFact;
  histFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
#endif
#ifdef G4ROOTANALYSIS_USE
  delete theROOTIonTuple;
  theROOTIonTuple = 0;

  delete theROOTNtuple;
  theROOTNtuple = 0;

  delete th14;
  th14 = 0;

  delete th13;
  th13 = 0;

  delete th12;
  th12 = 0;

  delete th11;
  th11 = 0;

  delete th10;
  th10 = 0;

  delete th9;
  th9 = 0;

  delete th8;
  th8 = 0;

  delete th7;
  th7 = 0;

  delete th6;
  th6 = 0;

  delete th5;
  th5 = 0;

  delete th4;
  th4 = 0;

  delete th3;
  th3 = 0;

  delete th2;
  th2 = 0;

  delete th1;
  th1 = 0;
#endif
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysisManager;
  return instance;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SetAnalysisFileName(G4String aFileName)
{
  this->AnalysisFileName = aFileName;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::book()
{
#ifdef G4ANALYSIS_USE
  // Build up  the  analysis factory
  aFact = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* treeFact = aFact -> createTreeFactory();

  // Create the .hbk or the .root file
  G4String fileName = "DoseDistribution.hbk";

  std::string opts = "export=root";

  theTree = treeFact -> create(fileName,"hbook",false,true);
  theTree = treeFact -> create(AnalysisFileName,"ROOT",false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFact;

  // Create the histogram and the ntuple factory
  histFact = aFact -> createHistogramFactory(*theTree);
  tupFact = aFact -> createTupleFactory(*theTree);

  // Create the histograms with the enrgy deposit along the X axis
  h1 = histFact -> createHistogram1D("10","slice, energy", 400, 0., 400. );

  h2 = histFact -> createHistogram1D("20","Secondary protons - slice, energy", 400, 0., 400. );

  h3 = histFact -> createHistogram1D("30","Secondary neutrons - slice, energy", 400, 0., 400. );

  h4 = histFact -> createHistogram1D("40","Secondary alpha - slice, energy", 400, 0., 400. );

  h5 = histFact -> createHistogram1D("50","Secondary gamma - slice, energy", 400, 0., 400. );

  h6 = histFact -> createHistogram1D("60","Secondary electron - slice, energy", 400, 0., 400. );

  h7 = histFact -> createHistogram1D("70","Secondary triton - slice, energy", 400, 0., 400. );

  h8 = histFact -> createHistogram1D("80","Secondary deuteron - slice, energy", 400, 0., 400. );

  h9 = histFact -> createHistogram1D("90","Secondary pion - slice, energy", 400, 0., 400. );

  h10 = histFact -> createHistogram1D("100","Energy distribution of secondary electrons", 70, 0., 70. );

  h11 = histFact -> createHistogram1D("110","Energy distribution of secondary photons", 70, 0., 70. );

  h12 = histFact -> createHistogram1D("120","Energy distribution of secondary deuterons", 70, 0., 70. );

  h13 = histFact -> createHistogram1D("130","Energy distribution of secondary tritons", 70, 0., 70. );

  h14 = histFact -> createHistogram1D("140","Energy distribution of secondary alpha particles", 70, 0., 70. );

  // Create the ntuple
  G4String columnNames = "int i; int j; int k; double energy;";
  G4String options = "";
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);

  // Create the ntuple
  G4String columnNames2 = "int a; double z;  int occupancy; double energy;";
  G4String options2 = "";
  if (tupFact) ionTuple = tupFact -> create("2","2", columnNames2, options2);
#endif
#ifdef G4ROOTANALYSIS_USE
  // Use ROOT
  theTFile = new TFile(AnalysisFileName, "RECREATE");

  // Create the histograms with the enrgy deposit along the X axis
  th1 = createHistogram1D("braggPeak","slice, energy", 400, 0., 400.);
  th2 = createHistogram1D("h20","Secondary protons - slice, energy", 400, 0., 400.);
  th3 = createHistogram1D("h30","Secondary neutrons - slice, energy", 400, 0., 400.);
  th4 = createHistogram1D("h40","Secondary alpha - slice, energy", 400, 0., 400.);
  th5 = createHistogram1D("h50","Secondary gamma - slice, energy", 400, 0., 400.);
  th6 = createHistogram1D("h60","Secondary electron - slice, energy", 400, 0., 400.);
  th7 = createHistogram1D("h70","Secondary triton - slice, energy", 400, 0., 400.);
  th8 = createHistogram1D("h80","Secondary deuteron - slice, energy", 400, 0., 400.);
  th9 = createHistogram1D("h90","Secondary pion - slice, energy", 400, 0., 400.);
  th10 = createHistogram1D("h100","Energy distribution of secondary electrons", 70, 0., 70.);
  th11 = createHistogram1D("h110","Energy distribution of secondary photons", 70, 0., 70.);
  th12 = createHistogram1D("h120","Energy distribution of secondary deuterons", 70, 0., 70.);
  th13 = createHistogram1D("h130","Energy distribution of secondary tritons", 70, 0., 70.);
  th14 = createHistogram1D("h140","Energy distribution of secondary alpha particles", 70, 0., 70.);

  theROOTNtuple = new TNtuple("theROOTNtuple", "Energy deposit by slice", "i/I:j/I:k/I:energy/F");
  theROOTIonTuple = new TNtuple("theROOTIonTuple", "Generic ion information", "a/I:z/F:occupancy/I:energy/F");
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::FillEnergyDeposit(G4int i,
						     G4int j,
						     G4int k,
						     G4double energy)
{
#ifdef G4ANALYSIS_USE
 if (ntuple)     {
       G4int iSlice = ntuple -> findColumn("i");
       G4int jSlice = ntuple -> findColumn("j");
       G4int kSlice = ntuple -> findColumn("k");
       G4int iEnergy = ntuple -> findColumn("energy");

       ntuple -> fill(iSlice,i);
       ntuple -> fill(jSlice,j);
       ntuple -> fill(kSlice,k);
       ntuple -> fill(iEnergy, energy);     }

  ntuple -> addRow();
#endif
#ifdef G4ROOTANALYSIS_USE
  if (theROOTNtuple) {
    theROOTNtuple->Fill(i, j, k, energy);
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::BraggPeak(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h1 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th1->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryProtonEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h2 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th2->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h3 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th3->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h4 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th4->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryGammaEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h5 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th5->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryElectronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h6 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th6->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryTritonEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h7 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th7->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h8 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th8->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryPionEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h9 -> fill(slice,energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th9->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::electronEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h10 -> fill(energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th10->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::gammaEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h11 -> fill(energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th11->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::deuteronEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h12 -> fill(energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th12->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::tritonEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h13 -> fill(energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th13->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::alphaEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h14 -> fill(energy);
#endif
#ifdef G4ROOTANALYSIS_USE
  th14->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::genericIonInformation(G4int a,
							 G4double z,
							 G4int electronOccupancy,
							 G4double energy)
{
#ifdef G4ANALYSIS_USE
  if (ionTuple)    {
       G4int aIndex = ionTuple -> findColumn("a");
       G4int zIndex = ionTuple -> findColumn("z");
       G4int electronIndex = ionTuple -> findColumn("occupancy");
       G4int energyIndex = ionTuple -> findColumn("energy");

       ionTuple -> fill(aIndex,a);
      ionTuple -> fill(zIndex,z);
      ionTuple -> fill(electronIndex, electronOccupancy);
       ionTuple -> fill(energyIndex, energy);
     }
   ionTuple -> addRow();
#endif
#ifdef G4ROOTANALYSIS_USE
   if (theROOTIonTuple) {
     theROOTIonTuple->Fill(a, z, electronOccupancy, energy);
   }
#endif
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::flush()
{
  HadrontherapyMatrix* matrix = HadrontherapyMatrix::getInstance();

  matrix->TotalEnergyDeposit();
#ifdef G4ANALYSIS_USE
  theTree -> commit();
  theTree ->close();
#endif
#ifdef G4ROOTANALYSIS_USE
  theROOTNtuple->Write();
  theROOTIonTuple->Write();
  theTFile->Write();
  theTFile->Clear();
  theTFile->Close();
#endif
  matrix->flush();
  //this->book();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::finish()
{
#ifdef G4ANALYSIS_USE
  // Write all histograms to file
  theTree -> commit();
  // Close (will again commit)
  theTree ->close();
#endif
#ifdef G4ROOTANALYSIS_USE
  theROOTNtuple->Write();
  theROOTIonTuple->Write();
  theTFile->Write();
  theTFile->Close();
#endif
}

#endif









