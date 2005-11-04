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
// $Id: PhotInRun.hh,v 1.3 2005-11-04 13:51:36 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

