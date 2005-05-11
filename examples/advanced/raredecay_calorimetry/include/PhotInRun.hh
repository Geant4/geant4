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
// $Id: PhotInRun.hh,v 1.1 2005-05-11 10:37:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef PhotInRun_h
#define PhotInRun_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "G4Allocator.hh"

class G4Event;
class G4HCtable;
class G4DCtable;

class PhotInRun : public G4Run
{
  public:
    PhotInRun();
    virtual ~PhotInRun();
    inline void *operator new(size_t);
    inline void operator delete(void* aRun);

  public:
    virtual void RecordEvent(const G4Event*);

  private:
    G4double totE[6];
    G4double totL[6];
    G4int nStep[6];
    G4int nGamma[6];
    G4int nElectron[6];
    G4int nPositron[6];
    G4double eMinGamma[6];
    G4double eMinElectron[6];
    G4double eMinPositron[6];

  public:
    inline G4double GetTotalE(G4int i) const
    { return totE[i]; }
    inline G4double GetTotalL(G4int i) const
    { return totL[i]; }
    inline G4int GetNStep(G4int i) const
    { return nStep[i]; }
    inline G4int GetNGamma(G4int i) const
    { return nGamma[i]; }
    inline G4int GetNElectron(G4int i) const
    { return nElectron[i]; }
    inline G4int GetNPositron(G4int i) const
    { return nPositron[i]; }
    inline G4double GetEMinGamma(G4int i) const
    { return eMinGamma[i]; }
    inline G4double GetEMinElectron(G4int i) const
    { return eMinElectron[i]; }
    inline G4double GetEMinPositron(G4int i) const
    { return eMinPositron[i]; }
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

