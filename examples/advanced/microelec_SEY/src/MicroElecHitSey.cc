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
/// \file MicroElecHitSey.cc
/// \brief Implementation of the MicroElecHitSey class
//
//-


#include "MicroElecHitSey.hh"
#include "G4UnitsTable.hh"
#ifdef	_G4CNF_USE_VIS_
	#include "G4VVisManager.hh"
#endif	// _G4CNF_USE_VIS_
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"

G4ThreadLocal G4Allocator<MicroElecHitSey>* MicroElecHitSeyAllocator=0;

MicroElecHitSey::MicroElecHitSey(){
    NbPrim = 999;
    NbSec = 999;
    NbSup50 = 999;
    trackID = 999;
    ParentID = 999;
    ParticleType = "initialization";
    ParticleName = "initialization";
    VolumeName = "initialization";
    Z = 999;
    A = 999;
    VertexKineticEnergy = 999.;
    PreStepKineticEnergy = 999.;
    PostStepKineticEnergy = 999.;
    Edep = 999.;
    Ni_Edep = 999.;
    PrePos = G4ThreeVector(0.0,0.0,0.0);
    PostPos = G4ThreeVector(0.0, 0.0, 0.0);
    PreStepMomentum = G4ThreeVector(0.0, 0.0, 0.0);
    PostStepMomentum = G4ThreeVector(0.0, 0.0, 0.0);
    StepLength = 999.;
}

MicroElecHitSey::~MicroElecHitSey(){
}

MicroElecHitSey::MicroElecHitSey(const MicroElecHitSey& right):G4VHit(){    
    NbPrim          = right.NbPrim;
    NbSec           = right.NbSec;
    NbSup50         = right.NbSup50;
    trackID    		= right.trackID;
    ParentID		= right.ParentID;
    ParticleType 	= right.ParticleType;
    ParticleName 	= right.ParticleName;
    VolumeName 		= right.VolumeName;
    Z				= right.Z;
    A				= right.A;
    VertexKineticEnergy	= right.VertexKineticEnergy;
    PreStepKineticEnergy = right.PreStepKineticEnergy;
    PostStepKineticEnergy = right.PostStepKineticEnergy;
    Edep       		= right.Edep;
	Ni_Edep       	= right.Ni_Edep;
    PrePos        	= right.PrePos;
    PostPos        	= right.PostPos;
    PreStepMomentum	= right.PreStepMomentum;
    PostStepMomentum = right.PostStepMomentum;
    StepLength 		= right.StepLength;

//-------------------
//New due to the use of touchables to identify the different sensitive detectors
	VolumeCopyNumber = right.VolumeCopyNumber;
}

const MicroElecHitSey& MicroElecHitSey::operator=(const MicroElecHitSey& right){
    NbPrim          = right.NbPrim;
    NbSec           = right.NbSec;
    NbSup50         = right.NbSup50;
    trackID    		= right.trackID;
    ParentID		= right.ParentID;
    ParticleType 	= right.ParticleType;
    ParticleName 	= right.ParticleName;
    VolumeName		= right.VolumeName;
    Z				= right.Z;
    A				= right.A;
    VertexKineticEnergy	= right.VertexKineticEnergy;
    PreStepKineticEnergy = right.PreStepKineticEnergy;
    PostStepKineticEnergy = right.PostStepKineticEnergy;
    Edep       		= right.Edep;
	Ni_Edep       	= right.Ni_Edep;
    PrePos        	= right.PrePos;
    PostPos        	= right.PostPos;
    PreStepMomentum	= right.PreStepMomentum;
    PostStepMomentum = right.PostStepMomentum;
    StepLength 		= right.StepLength;

	VolumeCopyNumber = right.VolumeCopyNumber;

    return *this;
}

int MicroElecHitSey::operator==(const MicroElecHitSey& ) const{
    return 0; 
}



void MicroElecHitSey::Print(){
    
    G4cout << "  [trackID, Parent ID]: [" << trackID << ", " << ParentID << "],  ParticleName :" << ParticleName<<", [PrimID, SecID]=["<<NbPrim<<", "<<NbSec << "]"<<G4endl;
    G4cout << "   Kinetic energy : " << G4BestUnit(PreStepKineticEnergy, "Energy") << ", position : " << G4BestUnit(PrePos, "Length") << G4endl;

}
