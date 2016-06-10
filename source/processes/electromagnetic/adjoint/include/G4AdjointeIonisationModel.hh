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
// $Id: G4AdjointeIonisationModel.hh 68044 2013-03-13 14:29:07Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointeIonisationModel
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	 September 2009 creation by L. Desorgher. Separate the concrete ionisation stuff from G4VEMAdjointModel  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint EM model for discrete reverse e- ionisation
//

#ifndef G4AdjointeIonisationModel_h
#define G4AdjointeIonisationModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4AdjointeIonisationModel: public G4VEmAdjointModel
{

public: //methods

//Constructor, destructor
  G4AdjointeIonisationModel();

  virtual ~G4AdjointeIonisationModel();

//Concrete implementation or virtual methods
  
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  
  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd, // kinetic energy of the secondary particle 
				      G4double Z, 
                                      G4double A = 0.);

				
private:  
  G4double DiffCrossSectionMoller(G4double kinEnergyProj,G4double kinEnergyProd);
private: //attributes   
   G4bool WithRapidSampling;
   

};
#endif

