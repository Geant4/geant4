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
#include "tls.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoRunAction;

class XrayFluoSensorHit : public G4VHit
{
public:
  
  XrayFluoSensorHit();
  ~XrayFluoSensorHit();
  XrayFluoSensorHit(const XrayFluoSensorHit&);
  const XrayFluoSensorHit& operator=(const XrayFluoSensorHit&);
  G4bool operator==(const XrayFluoSensorHit&) const;
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline void AddEnergy(G4double de)    {EdepTot += de;};
  void Draw();
  void Print();
  inline G4double GetEdepTot()      { return EdepTot;};
  
private:
  
  G4double EdepTot;
  
  G4double EdepDetect; 
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<XrayFluoSensorHit> XrayFluoSensorHitsCollection;

extern G4ThreadLocal G4Allocator<XrayFluoSensorHit> *XrayFluoSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* XrayFluoSensorHit::operator new(size_t)
{
  if (!XrayFluoSensorHitAllocator)
    XrayFluoSensorHitAllocator = new G4Allocator<XrayFluoSensorHit>;
  return (void*) XrayFluoSensorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoSensorHit::operator delete(void* aHit)
{
  XrayFluoSensorHitAllocator->FreeSingle((XrayFluoSensorHit*) aHit);
}

#endif



