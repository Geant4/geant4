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
// $Id: G4AdjointhMultipleScattering.hh 90259 2015-05-22 08:57:37Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointhMultipleScattering
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// GEANT4 Class header file
//
// File name:     G4AdjointhMultipleScattering
//
// Author:        Desorgher Laurent
//
// Creation date: 03.06.2009 cloned from G4hMultipleScattering by U.Laszlo with slight modification for adjoint_ion. 
//
//
//------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	03.06.2009  Creation by L. Desorgher. Cloned from G4hMultipleScattering by U.Laszlo with slight modifications.	  		
//		09.11.2009  Remove 	AlongStepGetPhysicalInteractionLength, to call the one of the base class.
//-------------------------------------------------------------
//	Documentation:
//		The class simulates the multiple scattering for adjoint proton of charged particle. In this approximate implementation the reverse multiple scattering 
//		is the same than the foward one. This should be changed in the future to have the MultipleScaterring computed for the energy  at the end of the step 
//		and not before the step. 
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4AdjointhMultipleScattering_h
#define G4AdjointhMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VMscModel;

class G4AdjointhMultipleScattering : public G4VMultipleScattering

{
public:    // with description

  G4AdjointhMultipleScattering(const G4String& processName="msc");

  virtual ~G4AdjointhMultipleScattering();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  // PrG4int few lines of informations about the process: validity range,
  void PrintInfo();

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

private:        // data members

  G4bool   isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

