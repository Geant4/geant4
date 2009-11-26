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
// $Id: HadrontherapyDetectorSD.cc; 
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy

#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyDetectorHit.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "HadrontherapyMatrix.hh"

HadrontherapyDetectorSD::HadrontherapyDetectorSD(G4String name):G4VSensitiveDetector(name)
{ 
    G4String HCname;
    collectionName.insert(HCname="HadrontherapyDetectorHitsCollection");
    HitsCollection = NULL; 
    G4String sensitiveDetectorName = name;

}

HadrontherapyDetectorSD::~HadrontherapyDetectorSD()
{ 
}

void HadrontherapyDetectorSD::Initialize(G4HCofThisEvent*)
{
    HitsCollection = new HadrontherapyDetectorHitsCollection(sensitiveDetectorName,
	    collectionName[0]);
}

G4bool HadrontherapyDetectorSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
    //The code doesn't seem to get here if we use the IAEA geometry. FIXME
    if(!ROhist)
	return false;

    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "DetectorPhys")
	return false;
    // Get kinetic energy
    G4Track * theTrack = aStep  ->  GetTrack();
    G4double KineticEnergy =  theTrack -> GetKineticEnergy();     

    G4ParticleDefinition *def = theTrack -> GetDefinition();
    //Get particle name  
    G4String particleName =  def -> GetParticleName();  
    // G4cout << def -> GetParticleType() << '\n';  
    // Get unique track_id (in an event)
    G4int TrackID = theTrack -> GetTrackID();

    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
    //if(energyDeposit == 0.) return false;

    // Read voxel indexes: i is the x index, k is the z index

    G4int k  = ROhist -> GetReplicaNumber(0);
    G4int i  = ROhist -> GetReplicaNumber(2);
    G4int j  = ROhist -> GetReplicaNumber(1);
    /*
       if (i==25) 
       {

    //if (ofs && particleName!="e-" && TrackID!=1) ofs << particleName << " " << TrackID << '\n'; }
    if (ofs) ofs << particleName << " " << TrackID << '\n'; }
    */
    HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
    // Increment Fluences
    // Hit voxel (marked with track_id)
    G4int* hitTrack = matrix -> GetHitTrack(i,j,k); // hitTrack must be cleared at every eventAction
if ( *hitTrack != TrackID )
{
    //G4cout << "TrackID " << TrackID << " Voxel " << i << '\t' << j << '\t' << k << G4endl;
    *hitTrack = TrackID;

    // Directly fill fluence data 
    matrix-> Fill(def, i, j, k, 0, true);
    // Fill kinetic energy histo
    //matrix -> fillEnergySpectrum(kineticEnergy, particleName,i,j,k);
}	 
if(energyDeposit != 0)                       
{  
    // Create a hit with the information of position is in the detector     
    HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();       
    detectorHit -> SetEdepAndPosition(i, j, k, energyDeposit); 
    HitsCollection -> insert(detectorHit);
}

// Energy deposit of secondary particles along X (integrated on Y and Z)

#ifdef ANALYSIS_USE
HadrontherapyAnalysisManager* analysis = 
HadrontherapyAnalysisManager::getInstance();
if(energyDeposit != 0)                       
{  

    // Fill DOSE matrix for single nuclide 
    matrix->Fill(def, i, j, k, energyDeposit/MeV);

    if(TrackID != 1)
    {
	if (particleName == "proton")
	    analysis -> SecondaryProtonEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "neutron")
	    analysis -> SecondaryNeutronEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "alpha")
	    analysis -> SecondaryAlphaEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "gamma")
	    analysis -> SecondaryGammaEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "e-")
	    analysis -> SecondaryElectronEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "triton")
	    analysis -> SecondaryTritonEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "deuteron")
	    analysis -> SecondaryDeuteronEnergyDeposit(i, energyDeposit/MeV);

	else if (particleName == "pi+" || particleName == "pi-" ||  particleName == "pi0")
	    analysis -> SecondaryPionEnergyDeposit(i, energyDeposit/MeV);   	
    }
}
#endif

return true;
}

void HadrontherapyDetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    static G4int HCID = -1;
    if(HCID < 0)
    { 
	HCID = GetCollectionID(0); 
    }
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

