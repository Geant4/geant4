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
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4LowEWentzelVIModel
//
// Author:        V.Ivanchenko 
//
// Creation date: 11.02.2014 from G4WentzelVIModel
//
// Modifications:
//
// Class Description:
//
// Implementation of the model of multiple scattering for low-energy e-

// -------------------------------------------------------------------
//

#ifndef G4LowEWentzelVIModel_h
#define G4LowEWentzelVIModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LowEWentzelVIModel : public G4WentzelVIModel
{

public:

  G4LowEWentzelVIModel();

  virtual ~G4LowEWentzelVIModel();

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4double& currentMinimalStep);

private:

  //  hide assignment operator
  G4LowEWentzelVIModel & operator=(const  G4LowEWentzelVIModel &right);
  G4LowEWentzelVIModel(const  G4LowEWentzelVIModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

