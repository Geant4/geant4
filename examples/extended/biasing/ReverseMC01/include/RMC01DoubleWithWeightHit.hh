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
/// \file biasing/ReverseMC01/include/RMC01DoubleWithWeightHit.hh
/// \brief Definition of the RMC01DoubleWithWeightHit class
//
//
//////////////////////////////////////////////////////////////
//  Class Name:           RMC01DoubleWithWeightHit
//        Author:               L. Desorgher
//        Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//        Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//                 17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RMC01DoubleWithWeightHit_h
#define RMC01DoubleWithWeightHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RMC01DoubleWithWeightHit : public G4VHit
{
public:
  RMC01DoubleWithWeightHit(G4double value, G4double weight);
  virtual ~RMC01DoubleWithWeightHit();
 
  RMC01DoubleWithWeightHit(const RMC01DoubleWithWeightHit &right);
 
  const RMC01DoubleWithWeightHit& operator = (const RMC01DoubleWithWeightHit &right);
 
  G4bool operator == (const RMC01DoubleWithWeightHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods to get the information - energy deposit and associated
  // position in the phantom - of the hits stored in the hits collection  
 
  inline G4double GetValue() {return fValue;}
  
  inline G4double GetWeight() {return fWeight;}

private:
 
  G4double fValue; // It can be anything
  G4double fWeight; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<RMC01DoubleWithWeightHit>
                                            RMC01DoubleWithWeightHitsCollection;
extern G4Allocator<RMC01DoubleWithWeightHit> RMC01DoubleWithWeightHitAllocator;

inline void* RMC01DoubleWithWeightHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RMC01DoubleWithWeightHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void RMC01DoubleWithWeightHit::operator delete(void *aHit)
{
  RMC01DoubleWithWeightHitAllocator.FreeSingle((RMC01DoubleWithWeightHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

