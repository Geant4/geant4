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
// $Id: G4IronStopping.hh,v 1.2 2008-08-27 08:52:19 antoni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4IronStopping_h
#define G4IronStopping_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4IronStopping
//
// Description: Data on stopping powers for light ions in compounds
//
// Author:      A.Ivantchenko 8.08.2008
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

class G4IronStopping
{
public:

  G4IronStopping(G4EmCorrections* corr = 0, G4bool splineFlag = true);

  ~G4IronStopping();

  G4double GetDEDX(G4int idxMaterial, G4double kinEnergy);

  G4double GetDEDX(const G4String& NameMaterial, G4double kinEnergy);

  G4int GetMaterialIndex(const G4String& NameMaterial);

  G4double GetDensity(G4int idx);

  G4String GetMaterialName(G4int idx);

  G4PhysicsVector* GetPhysicsVector(G4int idx);

  G4PhysicsVector* GetPhysicsVector(const G4String& NameMaterial);

private:

  void AddData(G4double* energy, G4double* stoppower, G4double factor);

  void Initialise(G4EmCorrections*);

  // hide assignment operator
  G4IronStopping & operator=(const  G4IronStopping &right);
  G4IronStopping(const  G4IronStopping&);

  G4bool spline;
  G4String MatName[16];
  G4double Density[16];

  std::vector<G4LPhysicsFreeVector*>  dedx;
};

#endif
 
