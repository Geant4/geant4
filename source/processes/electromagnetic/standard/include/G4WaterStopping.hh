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
// $Id: G4WaterStopping.hh,v 1.1 2006-05-13 17:55:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4WaterStopping_h
#define G4WaterStopping_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4WaterStopping
//
// Description: Data on stopping power
//
// Author:      V.Ivanchenko 12.05.2006
//
// Modifications:
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the ICRU73 report
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>

class G4EmCorrections;

class G4WaterStopping
{
public:

  G4WaterStopping(G4EmCorrections* corr = 0);

  ~G4WaterStopping();

  G4double GetElectronicDEDX(G4int Z, G4double energy);

private:

  void Initialise(G4EmCorrections*);

  // hide assignment operator
  G4WaterStopping & operator=(const  G4WaterStopping &right);
  G4WaterStopping(const  G4WaterStopping&);

  G4int    Z[8];
  G4int    A[8];
  std::vector<G4LPhysicsFreeVector*>  dedx;
};

#endif
