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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), 
//  G.A.P. Cirrone(d), F.Romano(d)
//
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTAnalysisManager.hh"
#include "IORTMatrix.hh"
#include "IORTAnalysisFileMessenger.hh"
#include <time.h>

IORTAnalysisManager* IORTAnalysisManager::instance = 0;

IORTAnalysisManager::IORTAnalysisManager() : 
  analysisFileName("DoseDistribution"),
  eventCounter(0)
{
  fMess = new IORTAnalysisFileMessenger(this);
}
/////////////////////////////////////////////////////////////////////////////

IORTAnalysisManager::~IORTAnalysisManager()
{
  if (fMess)
    delete fMess; 
  delete G4AnalysisManager::Instance(); 
}

/////////////////////////////////////////////////////////////////////////////

IORTAnalysisManager* IORTAnalysisManager::GetInstance()
{
  if (instance == 0) instance = new IORTAnalysisManager;
  return instance;
}

/////////////////////////////////////////////////////////////////////////////

void IORTAnalysisManager::SetAnalysisFileName(G4String aFileName)
{
  analysisFileName = aFileName;
}

/////////////////////////////////////////////////////////////////////////////

void IORTAnalysisManager::book()
{
  // Create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(1);
  man->SetFirstHistoId(1);
  man->SetFirstNtupleId(1);

  man->OpenFile(analysisFileName);
  
  // Create the histograms with the energy deposit along the X axis
  //ID=1 <different waterthicknesses are accoutned for in ROOT-analysis stage>
  man->CreateH1("braggPeak","slice, energy", 400, 0., 80); //
  //ID=2
  man->CreateH1("h20","Secondary protons - slice, energy", 400, 0., 400.);
  //ID=3
  man->CreateH1("h30","Secondary neutrons - slice, energy", 400, 0., 400.);
  //ID=4
  man->CreateH1("h40","Secondary alpha - slice, energy", 400, 0., 400.);
  //ID=5
  man->CreateH1("h50","Secondary gamma - slice, energy", 400, 0., 400.);
  //ID=6
  man->CreateH1("h60","Secondary electron - slice, energy", 400, 0., 400.);
  //ID=7
  man->CreateH1("h70","Secondary triton - slice, energy", 400, 0., 400.);
  //ID=8
  man->CreateH1("h80","Secondary deuteron - slice, energy", 400, 0., 400.);
  //ID=9
  man->CreateH1("h90","Secondary pion - slice, energy", 400, 0., 400.);
  //ID=10
  man->CreateH1("h100","Energy distribution of secondary electrons", 70, 0., 70.);
  //ID=11
  man->CreateH1("h110","Energy distribution of secondary photons", 70, 0., 70.);
  //ID=12
  man->CreateH1("h120","Energy distribution of secondary deuterons", 70, 0., 70.);
  //ID = 13
  man->CreateH1("h130","Energy distribution of secondary tritons", 70, 0., 70.);
  //ID = 14
  man->CreateH1("h140","Energy distribution of secondary alpha particles", 70, 0., 70.);
  //ID = 15
  man->CreateH1("heliumEnergyAfterPhantom",
		"Energy distribution of secondary helium fragments after the phantom",
		70, 0., 500.);
  //ID= 16
  man->CreateH1("hydrogenEnergyAfterPhantom",
		"Energy distribution of secondary helium fragments after the phantom",
		70, 0., 500.);

  //Now the ntuples
  //ID = 1
  man->CreateNtuple("kinFragNtuple", 
		    "Kinetic energy by voxel & fragment");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleIColumn("A");
  man->CreateNtupleDColumn("Z");
  man->CreateNtupleDColumn("kineticEnergy");
  man->FinishNtuple();

  //ID = 2
  man->CreateNtuple("kineticEnergyPrimaryNtuple", 
		    "Kinetic energy by voxel of primary");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleDColumn("kineticEnergy");
  man->FinishNtuple();

  //ID = 3
  man->CreateNtuple("doseFragNtuple",
		    "Energy deposit by voxel & fragment");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleIColumn("A");
  man->CreateNtupleDColumn("Z");
  man->CreateNtupleDColumn("energy");
  man->FinishNtuple();

  // ID =4
  man->CreateNtuple("fluenceFragNtuple", 
		    "Fluence by voxel & fragment");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleIColumn("A");
  man->CreateNtupleDColumn("Z");
  man->CreateNtupleDColumn("fluence");
  man->FinishNtuple();

  // ID=5
  man->CreateNtuple("letFragNtuple", 
		    "Let by voxel & fragment");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleIColumn("A");
  man->CreateNtupleDColumn("Z");
  man->CreateNtupleDColumn("letT");
  man->CreateNtupleDColumn("letD");
  man->FinishNtuple();

  //ID=6
  man->CreateNtuple("theROOTNtuple", 
		    "Energy deposit by slice");
  man->CreateNtupleIColumn("i");
  man->CreateNtupleIColumn("j");
  man->CreateNtupleIColumn("k");
  man->CreateNtupleDColumn("energy");
  man->FinishNtuple();

  //ID=7
  man->CreateNtuple("theROOTIonTuple",
		    "Generic ion information");
  man->CreateNtupleIColumn("a");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleIColumn("occupancy");
  man->CreateNtupleDColumn("energy");
  man->FinishNtuple();

  //ID=8
  man->CreateNtuple("fragmentNtuple",
		    "Fragments");
  man->CreateNtupleIColumn("A");
  man->CreateNtupleDColumn("Z");
  man->CreateNtupleDColumn("energy");
  man->CreateNtupleDColumn("posX");
  man->CreateNtupleDColumn("posY");
  man->CreateNtupleDColumn("posZ");
  man->FinishNtuple();

}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::FillEnergyDeposit(G4int i,
					    G4int j,
					    G4int k,
					    G4double energy)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(6,0,i);
  man->FillNtupleIColumn(6,1,j);
  man->FillNtupleIColumn(6,2,k);
  man->FillNtupleDColumn(6,3,energy);
  man->AddNtupleRow(6);  
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::BraggPeak(G4int slice, G4double energy)
{
  //FIXME
  G4AnalysisManager::Instance()->FillH1(1,slice,energy);
  //histo1->SetBinContent(slice, energy); //This uses setbincontent instead of fill to get labels correct
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryProtonEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(2,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(3,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(4,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryGammaEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(5,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryElectronEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(6,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryTritonEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(7,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(8,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryPionEnergyDeposit(G4int slice, G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(9,slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::electronEnergyDistribution(G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(10,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::gammaEnergyDistribution(G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(11,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::deuteronEnergyDistribution(G4double energy)
{
 G4AnalysisManager::Instance()->FillH1(12,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::tritonEnergyDistribution(G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(13,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::alphaEnergyDistribution(G4double energy)
{
  G4AnalysisManager::Instance()->FillH1(14,energy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::heliumEnergy(G4double secondaryParticleKineticEnergy)
{
  G4AnalysisManager::Instance()->FillH1(15,secondaryParticleKineticEnergy);
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::hydrogenEnergy(G4double secondaryParticleKineticEnergy)
{
  G4AnalysisManager::Instance()->FillH1(16,secondaryParticleKineticEnergy);
}

/////////////////////////////////////////////////////////////////////////////
// FillKineticFragmentTuple create an ntuple where the voxel indexs, the atomic number and mass and the kinetic
// energy of all the particles interacting with the phantom, are stored
void IORTAnalysisManager::FillKineticFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, G4double kinEnergy)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(1,0,i);
  man->FillNtupleIColumn(1,1,j);
  man->FillNtupleIColumn(1,2,k);
  man->FillNtupleIColumn(1,3,A);
  man->FillNtupleDColumn(1,4,Z);
  man->FillNtupleDColumn(1,5,kinEnergy);
  man->AddNtupleRow(1);  
}

/////////////////////////////////////////////////////////////////////////////
// FillKineticEnergyPrimaryNTuple creates a ntuple where the voxel indexs and the kinetic
// energies of ONLY primary particles interacting with the phantom, are stored
void IORTAnalysisManager::FillKineticEnergyPrimaryNTuple(G4int i, G4int j, G4int k, G4double kinEnergy)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(2,0,i);
  man->FillNtupleIColumn(2,1,j);
  man->FillNtupleIColumn(2,2,k);
  man->FillNtupleDColumn(2,3,kinEnergy);
  man->AddNtupleRow(2);  
}

/////////////////////////////////////////////////////////////////////////////
// This function is called only if ROOT is activated.
// It is called by the IORTMatric.cc class file and it is used to create two ntuples containing 
// the total energy deposited and the fluence values, in each voxel and per any particle (primary 
// and secondary particles beam) 
void IORTAnalysisManager::FillVoxelFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, 
						 G4double energy, G4double fluence)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(3,0,i);
  man->FillNtupleIColumn(3,1,j);
  man->FillNtupleIColumn(3,2,k);
  man->FillNtupleIColumn(3,3,A);
  man->FillNtupleDColumn(3,4,Z);
  man->FillNtupleDColumn(3,5,energy);
  man->AddNtupleRow(3);
  
  
  // Fill the ntuple containing the voxel, mass and atomic number and the fluence
  if (i==1 && Z==1) {
    man->FillNtupleIColumn(4,0,i);
    man->FillNtupleIColumn(4,1,j);
    man->FillNtupleIColumn(4,2,k);
    man->FillNtupleIColumn(4,3,A);
    man->FillNtupleDColumn(4,4,Z);
    man->FillNtupleDColumn(4,5,fluence);
    man->AddNtupleRow(4);    
  }
}

void IORTAnalysisManager::FillLetFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, 
					       G4double letT, G4double letD)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(5,0,i);
  man->FillNtupleIColumn(5,1,j);
  man->FillNtupleIColumn(5,2,k);
  man->FillNtupleIColumn(5,3,A);
  man->FillNtupleDColumn(5,4,Z);
  man->FillNtupleDColumn(5,5,letT);
  man->FillNtupleDColumn(5,6,letD);
  man->AddNtupleRow(5);  
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::FillFragmentTuple(G4int A, G4double Z, G4double energy, 
					    G4double posX, G4double posY, G4double posZ)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(8,0,A);
  man->FillNtupleDColumn(8,1,Z);
  man->FillNtupleDColumn(8,2,energy);
  man->FillNtupleDColumn(8,3,posX);
  man->FillNtupleDColumn(8,4,posY);
  man->FillNtupleDColumn(8,5,posZ);
  man->AddNtupleRow(8); 
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::genericIonInformation(G4int a,
						G4double z,
						G4int electronOccupancy,
						G4double energy)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(7,0,a);
  man->FillNtupleDColumn(7,1,z);
  man->FillNtupleIColumn(7,2,electronOccupancy);
  man->FillNtupleDColumn(7,3,energy);
  man->AddNtupleRow(7);  
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::startNewEvent()
{
  eventCounter++;
}
/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::setGeometryMetaData(G4double endDetectorPosition, G4double waterThickness, 
					      G4double phantomCenter)
{
  detectorDistance = endDetectorPosition;
  phantomDepth = waterThickness;
  phantomCenterDistance = phantomCenter;
}

/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::setBeamMetaData(G4double meanKineticEnergy,G4double sigmaEnergy)
{
  beamEnergy = meanKineticEnergy;
  energyError = sigmaEnergy;
}
/////////////////////////////////////////////////////////////////////////////
// Flush data & close the file
void IORTAnalysisManager::flush()
{
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  eventCounter = 0;
}

