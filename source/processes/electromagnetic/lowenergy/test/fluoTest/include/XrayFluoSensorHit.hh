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
//
// $Id: XrayFluoSensorHit.hh
// GEANT4 tag $Name: 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoSensorHit_h
#define XrayFluoSensorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <CLHEP/Random/Randomize.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoRunAction;
class XrayFluoSensorHit : public G4VHit
{
public:

  XrayFluoSensorHit();
  ~XrayFluoSensorHit();
  XrayFluoSensorHit(const XrayFluoSensorHit&);
      const XrayFluoSensorHit& operator=(const XrayFluoSensorHit&);
      int operator==(const XrayFluoSensorHit&) const;
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  void AddEnergy(G4double de)    {EdepTot += de;};

  void Draw();
  void Print();
     G4double GetEdepTot()      { return EdepTot;};
  G4double GetEdepDetect()   { return EdepDetect;};
 
private:

   G4double EdepTot;
 
  G4double EdepDetect; 
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<XrayFluoSensorHit> XrayFluoSensorHitsCollection;

extern G4Allocator<XrayFluoSensorHit> XrayFluoSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* XrayFluoSensorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) XrayFluoSensorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoSensorHit::operator delete(void* aHit)
{
  XrayFluoSensorHitAllocator.FreeSingle((XrayFluoSensorHit*) aHit);
}

#endif



