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

#ifdef G4ANALYSIS_USE_ROOT
#undef G4ANALYSIS_USE
#endif

/////////////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_USE
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() :
  analysisFileName("DoseDistribution.root"), aFact(0), theTree(0), histFact(0), tupFact(0), h1(0), h2(0), h3(0),
  h4(0), h5(0), h6(0), h7(0), h8(0), h9(0), h10(0), h11(0), h12(0), h13(0), h14(0), h15(0), h16(0), ntuple(0),
  ionTuple(0),
  fragmentTuple(0),
  eventCounter(0)
{
	fMess = new HadrontherapyAnalysisFileMessenger(this);
}
#endif
#ifdef G4ANALYSIS_USE_ROOT
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() :
  analysisFileName("DoseDistribution.root"),theTFile(0), histo1(0), histo2(0), histo3(0),
  histo4(0), histo5(0), histo6(0), histo7(0), histo8(0), histo9(0), histo10(0), histo11(0), histo12(0), histo13(0), histo14(0), histo15(0), histo16(0),
  theROOTNtuple(0),
  theROOTIonTuple(0),
  fragmentNtuple(0),
  metaData(0),
  eventCounter(0)
{
  fMess = new HadrontherapyAnalysisFileMessenger(this);
}
#endif
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager()
{
delete(fMess); //kill the messenger
#ifdef G4ANALYSIS_USE
  delete fragmentTuple;
  fragmentTuple = 0;

  delete ionTuple;
  ionTuple = 0;

  delete ntuple;
  ntuple = 0;

  delete h16;
  h16 = 0;

  delete h15;
  h15 = 0;

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
#ifdef G4ANALYSIS_USE_ROOT
  delete metaData;
  metaData = 0;

  delete fragmentNtuple;
  fragmentNtuple = 0;

  delete theROOTIonTuple;
  theROOTIonTuple = 0;

  delete theROOTNtuple;
  theROOTNtuple = 0;
  
  delete histo16;
  histo14 = 0;

  delete histo15;
  histo14 = 0;

  delete histo14;
  histo14 = 0;

  delete histo13;
  histo13 = 0;

  delete histo12;
  histo12 = 0;

  delete histo11;
  histo11 = 0;

  delete histo10;
  histo10 = 0;

  delete histo9;
  histo9 = 0;

  delete histo8;
  histo8 = 0;

  delete histo7;
  histo7 = 0;

  delete histo6;
  histo6 = 0;

  delete histo5;
  histo5 = 0;

  delete histo4;
  histo4 = 0;

  delete histo3;
  histo3 = 0;

  delete histo2;
  histo2 = 0;

  delete histo1;
  histo1 = 0;
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
  this->analysisFileName = aFileName;
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
  theTree = treeFact -> create(analysisFileName,"ROOT",false,true,opts);

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

  h15 = histFact -> createHistogram1D("150","Energy distribution of helium fragments after the phantom", 70, 0., 500.);

  h16 = histFact -> createHistogram1D("160","Energy distribution of hydrogen fragments after the phantom", 70, 0., 500.);

  // Create the ntuple
  G4String columnNames = "int i; int j; int k; double energy;";
  G4String options = "";
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);

  // Create the ntuple
  G4String columnNames2 = "int a; double z;  int occupancy; double energy;";
  G4String options2 = "";
  if (tupFact) ionTuple = tupFact -> create("2","2", columnNames2, options2);

  // Create the fragment ntuple
  G4String columnNames3 = "int a; double z; double energy; double posX; double posY; double posZ;";
  G4String options3 = "";
  if (tupFact) fragmentTuple = tupFact -> create("3","3", columnNames3, options3);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  // Use ROOT
  theTFile = new TFile(analysisFileName, "RECREATE");

  // Create the histograms with the enrgy deposit along the X axis
  histo1 = createHistogram1D("braggPeak","slice, energy", 400, 0., 400.);
  histo2 = createHistogram1D("h20","Secondary protons - slice, energy", 400, 0., 400.);
  histo3 = createHistogram1D("h30","Secondary neutrons - slice, energy", 400, 0., 400.);
  histo4 = createHistogram1D("h40","Secondary alpha - slice, energy", 400, 0., 400.);
  histo5 = createHistogram1D("h50","Secondary gamma - slice, energy", 400, 0., 400.);
  histo6 = createHistogram1D("h60","Secondary electron - slice, energy", 400, 0., 400.);
  histo7 = createHistogram1D("h70","Secondary triton - slice, energy", 400, 0., 400.);
  histo8 = createHistogram1D("h80","Secondary deuteron - slice, energy", 400, 0., 400.);
  histo9 = createHistogram1D("h90","Secondary pion - slice, energy", 400, 0., 400.);
  histo10 = createHistogram1D("h100","Energy distribution of secondary electrons", 70, 0., 70.);
  histo11 = createHistogram1D("h110","Energy distribution of secondary photons", 70, 0., 70.);
  histo12 = createHistogram1D("h120","Energy distribution of secondary deuterons", 70, 0., 70.);
  histo13 = createHistogram1D("h130","Energy distribution of secondary tritons", 70, 0., 70.);
  histo14 = createHistogram1D("h140","Energy distribution of secondary alpha particles", 70, 0., 70.);
  histo15 = createHistogram1D("heliumEnergyAfterPhantom","Energy distribution of secondary helium fragments after the phantom",
			   70, 0., 500.);
  histo16 = createHistogram1D("hydrogenEnergyAfterPhantom","Energy distribution of secondary helium fragments after the phantom",
			   70, 0., 500.);

  theROOTNtuple = new TNtuple("theROOTNtuple", "Energy deposit by slice", "i:j:k:energy");
  theROOTIonTuple = new TNtuple("theROOTIonTuple", "Generic ion information", "a:z:occupancy:energy");
  fragmentNtuple = new TNtuple("fragmentNtuple", "Fragments", "A:Z:energy:posX:posY:posZ");
  metaData = new TNtuple("metaData", "Metadata", "events");
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
#ifdef G4ANALYSIS_USE_ROOT
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
#ifdef G4ANALYSIS_USE_ROOT
  histo1->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryProtonEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h2 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo2->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h3 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo3->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h4 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo4->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryGammaEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h5 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo5->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryElectronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h6 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo6->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryTritonEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h7 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo7->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h8 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo8->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryPionEnergyDeposit(G4int slice, G4double energy)
{
#ifdef G4ANALYSIS_USE
  h9 -> fill(slice,energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo9->Fill(slice, energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::electronEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h10 -> fill(energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo10->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::gammaEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h11 -> fill(energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo11->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::deuteronEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h12 -> fill(energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo12->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::tritonEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h13 -> fill(energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo13->Fill(energy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::alphaEnergyDistribution(G4double energy)
{
#ifdef G4ANALYSIS_USE
  h14 -> fill(energy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo14->Fill(energy);
#endif
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::heliumEnergy(G4double secondaryParticleKineticEnergy)
{
#ifdef G4ANALYSIS_USE
  h15->fill(secondaryParticleKineticEnergy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo15->Fill(secondaryParticleKineticEnergy);
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::hydrogenEnergy(G4double secondaryParticleKineticEnergy)
{
#ifdef G4ANALYSIS_USE
  h16->fill(secondaryParticleKineticEnergy);
#endif
#ifdef G4ANALYSIS_USE_ROOT
  histo16->Fill(secondaryParticleKineticEnergy);
#endif
}



void HadrontherapyAnalysisManager::fillFragmentTuple(G4int A, G4double Z, G4double energy, G4double posX, G4double posY, G4double posZ)
{
#ifdef G4ANALYSIS_USE
  if (fragmentTuple)    {
       G4int aIndex = fragmentTuple -> findColumn("a");
       G4int zIndex = fragmentTuple -> findColumn("z");
       G4int energyIndex = fragmentTuple -> findColumn("energy");
       G4int posXIndex = fragmentTuple -> findColumn("posX");
       G4int posYIndex = fragmentTuple -> findColumn("posY");
       G4int posZIndex = fragmentTuple -> findColumn("posZ");

       fragmentTuple -> fill(aIndex,A);
       fragmentTuple -> fill(zIndex,Z);
       fragmentTuple -> fill(energyIndex, energy);
       fragmentTuple -> fill(posXIndex, posX);
       fragmentTuple -> fill(posYIndex, posY);
       fragmentTuple -> fill(posZIndex, posZ);
       fragmentTuple -> addRow();
  }
#endif
#ifdef G4ANALYSIS_USE_ROOT
  //G4cout <<" A = " << A << "  Z = " << Z << " energy = " << energy << G4endl;
  fragmentNtuple->Fill(A, Z, energy, posX, posY, posZ);
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
#ifdef G4ANALYSIS_USE_ROOT
   if (theROOTIonTuple) {
     theROOTIonTuple->Fill(a, z, electronOccupancy, energy);
   }
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::startNewEvent()
{
  eventCounter++;
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
#ifdef G4ANALYSIS_USE_ROOT
  metaData->Fill((Float_t) eventCounter);
  metaData->Write();
  theROOTNtuple->Write();
  theROOTIonTuple->Write();
  fragmentNtuple->Write();
  theTFile->Write();
  //  theTFile->Clear();
  theTFile->Close();
#endif
  eventCounter = 0;
  matrix->flush();
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
#ifdef G4ANALYSIS_USE_ROOT
  metaData->Fill((Float_t) eventCounter);
  metaData->Write();
  theROOTNtuple->Write();
  theROOTIonTuple->Write();
  fragmentNtuple->Write();
  theTFile->Write();
  theTFile->Close();
#endif
  eventCounter = 0;
}

#endif









