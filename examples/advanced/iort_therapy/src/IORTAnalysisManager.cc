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
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
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

IORTAnalysisManager::IORTAnalysisManager() 
#ifdef G4ANALYSIS_USE_ROOT 
:
analysisFileName("DoseDistribution.root"),theTFile(0), histo1(0), histo2(0), histo3(0),
histo4(0), histo5(0), histo6(0), histo7(0), histo8(0), histo9(0), histo10(0), histo11(0), histo12(0), histo13(0), histo14(0), histo15(0), histo16(0),
kinFragNtuple(0),
kineticEnergyPrimaryNtuple(0),
doseFragNtuple(0),
fluenceFragNtuple(0),
letFragNtuple(0),
theROOTNtuple(0),
theROOTIonTuple(0),
fragmentNtuple(0),
metaData(0),
eventCounter(0)
#endif
{
#ifdef G4ANALYSIS_USE_ROOT
    fMess = new IORTAnalysisFileMessenger(this);
#else
    fMess = new IORTAnalysisFileMessenger();
#endif
}
/////////////////////////////////////////////////////////////////////////////

IORTAnalysisManager::~IORTAnalysisManager()
{
    delete fMess; 
#ifdef G4ANALYSIS_USE_ROOT
    Clear();	
#endif
}

IORTAnalysisManager* IORTAnalysisManager::GetInstance()
{
	if (instance == 0) instance = new IORTAnalysisManager;
	return instance;
}
#ifdef G4ANALYSIS_USE_ROOT
void IORTAnalysisManager::Clear()
{
    if (theTFile)
    {
	delete metaData;
	metaData = 0;

	delete fragmentNtuple;
	fragmentNtuple = 0;

	delete theROOTIonTuple;
	theROOTIonTuple = 0;

	delete theROOTNtuple;
	theROOTNtuple = 0;

	delete histo16;
	histo16 = 0;

	delete histo15;
	histo15 = 0;

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
    }
}
/////////////////////////////////////////////////////////////////////////////

void IORTAnalysisManager::SetAnalysisFileName(G4String aFileName)
{
	this->analysisFileName = aFileName;
}

	/////////////////////////////////////////////////////////////////////////////
G4bool IORTAnalysisManager::IsTheTFile()
{
    return (theTFile) ? true:false; 
}
void IORTAnalysisManager::book()
{
	delete theTFile; // this is similar to theTFile->Close() => delete all associated variables created via new, moreover it delete itself.  

	theTFile = new TFile(analysisFileName, "RECREATE");

	// Create the histograms with the energy deposit along the X axis
	histo1 = createHistogram1D("braggPeak","slice, energy", 400, 0., 80); //<different waterthicknesses are accoutned for in ROOT-analysis stage
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

	kinFragNtuple  = new TNtuple("kinFragNtuple", 
		"Kinetic energy by voxel & fragment", 
		"i:j:k:A:Z:kineticEnergy");
	kineticEnergyPrimaryNtuple= new TNtuple("kineticEnergyPrimaryNtuple", 
		"Kinetic energy by voxel of primary", 
		"i:j:k:kineticEnergy");
	doseFragNtuple = new TNtuple("doseFragNtuple",
		"Energy deposit by voxel & fragment",
		"i:j:k:A:Z:energy");

	fluenceFragNtuple = new TNtuple("fluenceFragNtuple", 
		"Fluence by voxel & fragment",
		"i:j:k:A:Z:fluence");

	letFragNtuple = new TNtuple("letFragNtuple", 
		"Let by voxel & fragment",
		"i:j:k:A:Z:letT:letD");

	theROOTNtuple =   new TNtuple("theROOTNtuple", 
		"Energy deposit by slice",
		"i:j:k:energy");

	theROOTIonTuple = new TNtuple("theROOTIonTuple",
		"Generic ion information",
		"a:z:occupancy:energy");

	fragmentNtuple =  new TNtuple("fragmentNtuple",
		"Fragments",
		"A:Z:energy:posX:posY:posZ");

	metaData =        new TNtuple("metaData",
		"Metadata",
		"events:detectorDistance:waterThickness:beamEnergy:energyError:phantomCenterDistance");
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::FillEnergyDeposit(G4int i,
						     G4int j,
						     G4int k,
						     G4double energy)
{
	if (theROOTNtuple) 
	{
		theROOTNtuple->Fill(i, j, k, energy);
	}
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::BraggPeak(G4int slice, G4double energy)
{
	histo1->SetBinContent(slice, energy); //This uses setbincontent instead of fill to get labels correct
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryProtonEnergyDeposit(G4int slice, G4double energy)
{
	histo2->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy)
{
	histo3->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy)
{
	histo4->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryGammaEnergyDeposit(G4int slice, G4double energy)
{
	histo5->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryElectronEnergyDeposit(G4int slice, G4double energy)
{
	histo6->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryTritonEnergyDeposit(G4int slice, G4double energy)
{
	histo7->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy)
{
	histo8->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::SecondaryPionEnergyDeposit(G4int slice, G4double energy)
{
	histo9->Fill(slice, energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::electronEnergyDistribution(G4double energy)
{
	histo10->Fill(energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::gammaEnergyDistribution(G4double energy)
{
	histo11->Fill(energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::deuteronEnergyDistribution(G4double energy)
{
	histo12->Fill(energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::tritonEnergyDistribution(G4double energy)
{
	histo13->Fill(energy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::alphaEnergyDistribution(G4double energy)
{
	histo14->Fill(energy);
}
	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::heliumEnergy(G4double secondaryParticleKineticEnergy)
{
	histo15->Fill(secondaryParticleKineticEnergy);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::hydrogenEnergy(G4double secondaryParticleKineticEnergy)
{
	histo16->Fill(secondaryParticleKineticEnergy);
}

	/////////////////////////////////////////////////////////////////////////////
	// FillKineticFragmentTuple create an ntuple where the voxel indexs, the atomic number and mass and the kinetic
	// energy of all the particles interacting with the phantom, are stored
void IORTAnalysisManager::FillKineticFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, G4double kinEnergy)
{
    kinFragNtuple -> Fill(i, j, k, A, Z, kinEnergy);
}

	/////////////////////////////////////////////////////////////////////////////
	// FillKineticEnergyPrimaryNTuple creates a ntuple where the voxel indexs and the kinetic
	// energies of ONLY primary particles interacting with the phantom, are stored
void IORTAnalysisManager::FillKineticEnergyPrimaryNTuple(G4int i, G4int j, G4int k, G4double kinEnergy)
{
    kineticEnergyPrimaryNtuple -> Fill(i, j, k, kinEnergy);
}

	/////////////////////////////////////////////////////////////////////////////
	// This function is called only if ROOT is activated.
	// It is called by the IORTMatric.cc class file and it is used to create two ntuples containing 
	// the total energy deposited and the fluence values, in each voxel and per any particle (primary 
	// and secondary particles beam) 
void IORTAnalysisManager::FillVoxelFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, G4double energy, G4double fluence)
{
		// Fill the ntuple containing the voxel, mass and atomic number and the energy deposited
    doseFragNtuple ->    Fill( i, j, k, A, Z, energy );
	
		// Fill the ntuple containing the voxel, mass and atomic number and the fluence
	if (i==1 && Z==1) {
		fluenceFragNtuple -> Fill( i, j, k, A, Z, fluence );

	}
}

void IORTAnalysisManager::FillLetFragmentTuple(G4int i, G4int j, G4int k, G4int A, G4double Z, G4double letT, G4double letD)
{
	letFragNtuple -> Fill( i, j, k, A, Z, letT, letD);

}
	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::FillFragmentTuple(G4int A, G4double Z, G4double energy, G4double posX, G4double posY, G4double posZ)
{
	fragmentNtuple->Fill(A, Z, energy, posX, posY, posZ);
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::genericIonInformation(G4int a,
														 G4double z,
														 G4int electronOccupancy,
														 G4double energy)
{
	if (theROOTIonTuple) {
		theROOTIonTuple->Fill(a, z, electronOccupancy, energy);
	}
}

	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::startNewEvent()
{
	eventCounter++;
}
	/////////////////////////////////////////////////////////////////////////////
void IORTAnalysisManager::setGeometryMetaData(G4double endDetectorPosition, G4double waterThickness, G4double phantomCenter)
{
	this->detectorDistance = endDetectorPosition;
	this->phantomDepth = waterThickness;
	this->phantomCenterDistance = phantomCenter;
}
void IORTAnalysisManager::setBeamMetaData(G4double meanKineticEnergy,G4double sigmaEnergy)
{
	this->beamEnergy = meanKineticEnergy;
	this->energyError = sigmaEnergy;
}
/////////////////////////////////////////////////////////////////////////////
// Flush data & close the file
void IORTAnalysisManager::flush()
{
    if (theTFile)
    {
	theTFile -> Write(); 
	theTFile -> Close();
    }
    theTFile = 0;
    eventCounter = 0;
}

#endif
