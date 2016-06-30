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
// $Id: G4WaterStopping.hh 96934 2016-05-18 09:10:41Z gcosmo $

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

  explicit G4WaterStopping(G4EmCorrections* corr = 0, G4bool splineFlag = true);

  ~G4WaterStopping();

  G4double GetElectronicDEDX(G4int Z, G4double energy);

private:

  void Initialise(G4EmCorrections*);

  void AddData(const G4double* energy, const G4double* stoppower, 
	       G4double factor);

  // hide assignment operator
  G4WaterStopping & operator=(const  G4WaterStopping &right) = delete;
  G4WaterStopping(const  G4WaterStopping&) = delete;

  G4bool   spline;
  static const G4int Z[17];
  static const G4double A[17];
  G4double emin;
  std::vector<G4LPhysicsFreeVector*>  dedx;
};

#endif
