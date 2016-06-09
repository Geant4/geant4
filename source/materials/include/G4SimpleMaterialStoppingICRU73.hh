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
// $Id: G4SimpleMaterialStoppingICRU73.hh,v 1.8 2009/11/09 16:51:07 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $

#ifndef G4SimpleMaterialStoppingICRU73_h
#define G4SimpleMaterialStoppingICRU73_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4SimpleMaterialStoppingICRU73
//
// Description: Data on stopping powers for light ions in compounds
//
// Author:      Ivantchenko 10.07.2008
//
// in the framework of the ESA Technology Research Programme
// (ESA contract 21435/08/NL/AT)
//
// Modifications:
// 03.11.2009 A. Lechner:  Added new methods BuildPhysicsVector according
//            to interface changes in base class G4VIonDEDXTable.
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the ICRU73 report
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"
#include "G4VIonDEDXTable.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>
#include <utility>
#include <map>


class G4SimpleMaterialStoppingICRU73 : public G4VIonDEDXTable
{
public:

  G4SimpleMaterialStoppingICRU73(G4bool splineFlag = true);

  ~G4SimpleMaterialStoppingICRU73();

  G4bool BuildPhysicsVector(G4int ionZ, 
                            const G4String& matName);

  G4bool BuildPhysicsVector(G4int ionZ, 
                            G4int matZ);

  G4bool IsApplicable(G4int ionZ,  
                      G4int matZ);

  G4bool IsApplicable(G4int ionZ, 
                      const G4String& matName);

  G4PhysicsVector* GetPhysicsVector(G4int ionZ, 
                                    G4int matZ);

  G4PhysicsVector* GetPhysicsVector(G4int ionZ, 
                               	    const G4String& matName);

  G4double GetDEDX(G4double kinEnergyPerNucleon,
                   G4int ionZ,
                   const G4String& matName);

  G4double GetDEDX(G4double kinEnergyPerNucleon,
                   G4int ionZ,
                   G4int matZ);

private:
  // Function for creating a physics vector
  G4PhysicsVector* CreatePhysicsVector(G4double* energy, 
                                       G4double* stoppower, 
                                       G4double factor);

  // Assignment operator and copy constructor
  G4SimpleMaterialStoppingICRU73 & 
                       operator=(const G4SimpleMaterialStoppingICRU73 &right);
  G4SimpleMaterialStoppingICRU73(const G4SimpleMaterialStoppingICRU73&);

  // Flag indicating the use of spline interpolation for dE/dx vectors
  G4bool spline;

  // Minimum and maximum atomic number
  G4int minIonAtomicNmb;
  G4int maxIonAtomicNmb;

  // Vectors containing the atomic numbers and names of the materials
  std::vector<G4int> atomicNumbersMat;
  std::vector<G4String> namesMat;

  // Keys (ion atomic number and material names) corresponding to created 
  // dE/dx vectors
  typedef std::pair<G4int, G4String> DEDXKey;
  std::vector<DEDXKey> dedxKeys;

  // Vector of dE/dx vectors
  std::vector<G4PhysicsVector*>  dedx;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
