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
/// \file biasing/ReverseMC01/src/RMC01DoubleWithWeightHit.cc
/// \brief Implementation of the RMC01DoubleWithWeightHit class
//
//
//////////////////////////////////////////////////////////////
//      Class Name:        RMC01DoubleWithWeightHit
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RMC01DoubleWithWeightHit.hh"

G4Allocator<RMC01DoubleWithWeightHit> RMC01DoubleWithWeightHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DoubleWithWeightHit::RMC01DoubleWithWeightHit(G4double aValue,
                                                   G4double aWeight)
: G4VHit(), fValue(aValue), fWeight(aWeight)
{;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DoubleWithWeightHit::~RMC01DoubleWithWeightHit()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DoubleWithWeightHit::RMC01DoubleWithWeightHit(
                                                 const RMC01DoubleWithWeightHit &right)
  : G4VHit()
{
  fValue = right.fValue;
  fWeight = right.fWeight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const RMC01DoubleWithWeightHit& RMC01DoubleWithWeightHit::operator=(
                                                 const RMC01DoubleWithWeightHit &right)
{
  fValue = right.fValue;
  fWeight = right.fWeight;
 return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RMC01DoubleWithWeightHit::operator==
                                         (const RMC01DoubleWithWeightHit &right) const
{
 return(fValue == right.fValue && fWeight == right.fWeight);
}

