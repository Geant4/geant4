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
// $Id: PhotInRun.hh,v 1.4 2006/06/29 16:24:55 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#ifndef PhotInRun_h
#define PhotInRun_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "G4Allocator.hh"

#include "PhotInConstants.hh"
#include "PhotInCalorHit.hh"
#include "PhotInStackingAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

class PhotInRun : public G4Run
{
  public:
    PhotInRun();
    virtual ~PhotInRun();
    void *operator new(size_t);
    void operator delete(void* aRun);

    G4double GetTotalE(G4int i) const       { return totE[i]; }
    G4double GetTotalL(G4int i) const       { return totL[i]; }
    G4int    GetNStep(G4int i) const        { return nStep[i]; }
    G4int    GetNGamma(G4int i) const       { return nGamma[i]; }
    G4int    GetNElectron(G4int i) const    { return nElectron[i]; }
    G4int    GetNPositron(G4int i) const    { return nPositron[i]; }
    G4double GetEMinGamma(G4int i) const    { return eMinGamma[i]; }
    G4double GetEMinElectron(G4int i) const { return eMinElectron[i]; }
    G4double GetEMinPositron(G4int i) const { return eMinPositron[i]; }

    virtual void RecordEvent(const G4Event*);

  private:
    G4double totE[PhotInDiNSections];
    G4double totL[PhotInDiNSections];
    G4int nStep[PhotInDiNSections];
    G4int nGamma[PhotInDiNSections];
    G4int nElectron[PhotInDiNSections];
    G4int nPositron[PhotInDiNSections];
    G4int nProton[PhotInDiNSections];
    G4int nNeutron[PhotInDiNSections];
    G4int nPion[PhotInDiNSections];
    G4int nKaon[PhotInDiNSections];
    G4double eMinGamma[PhotInDiNSections];
    G4double eMinElectron[PhotInDiNSections];
    G4double eMinPositron[PhotInDiNSections];
    G4double eMinProton[PhotInDiNSections];
    G4double eMinNeutron[PhotInDiNSections];
    G4double eMinPion[PhotInDiNSections];
    G4double eMinKaon[PhotInDiNSections];
    G4double eMaxGamma[PhotInDiNSections];
    G4double eMaxElectron[PhotInDiNSections];
    G4double eMaxPositron[PhotInDiNSections];
    G4double eMaxProton[PhotInDiNSections];
    G4double eMaxNeutron[PhotInDiNSections];
    G4double eMaxPion[PhotInDiNSections];
    G4double eMaxKaon[PhotInDiNSections];
};

extern G4Allocator<PhotInRun> anPhotInRunAllocator;

inline void* PhotInRun::operator new(size_t)
{
  void* aRun;
  aRun = (void*)anPhotInRunAllocator.MallocSingle();
  return aRun;
}

inline void PhotInRun::operator delete(void* aRun)
{
  anPhotInRunAllocator.FreeSingle((PhotInRun*)aRun);
}

#endif

