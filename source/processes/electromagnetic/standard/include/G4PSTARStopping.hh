//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4PSTARStopping.hh,v 1.2 2006-05-18 18:26:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4PSTARStopping_h
#define G4PSTARStopping_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4PSTARStopping
//
// Description: Data on stopping power
//
// Author:      A.Ivantchenko 18.04.2006
//
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            CSMAN-5288
//
// Modifications:
// 
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the NIST Data Base  
// http://physics.nist.gov/PhysRefData/STAR/Text/PSTAR.html
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"

class G4Material;

class G4PSTARStopping 
{ 
public: 

  G4PSTARStopping();

  ~G4PSTARStopping();

  G4int GetIndex(const G4Material*);

  G4double GetElectronicDEDX(G4int idx, G4double energy);
  G4double GetElectronicDEDX(const G4Material*, G4double energy);

private:

  void Initialise();

  // hide assignment operator
  G4PSTARStopping & operator=(const  G4PSTARStopping &right);
  G4PSTARStopping(const  G4PSTARStopping&);

  const G4Material* currentMaterial;
  G4int index, matIndex;
  G4String name[74];
  G4double currentE, res;
  G4double e[74][60], kinE[60]; 
  G4int Znum[74];
};

inline G4double G4PSTARStopping::GetElectronicDEDX(const G4Material* mat, 
						   G4double energy)
{
  return GetElectronicDEDX(GetIndex(mat), energy);
}

#endif
