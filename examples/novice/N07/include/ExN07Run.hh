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
// $Id: ExN07Run.hh,v 1.1 2003-03-10 01:43:36 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN07Run_h
#define ExN07Run_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "G4Allocator.hh"

class G4Event;
class G4HCtable;
class G4DCtable;

class ExN07Run : public G4Run
{
  public:
    ExN07Run();
    virtual ~ExN07Run();
    inline void *operator new(size_t);
    inline void operator delete(void* aRun);

  public:
    virtual void RecordEvent(G4Event*);

  private:
    G4double totE[6];
    G4double totL[6];
    G4int nStep[6];

  public:
    inline G4double GetTotalE(G4int i) const
    { return totE[i]; }
    inline G4double GetTotalL(G4int i) const
    { return totL[i]; }
    inline G4int GetNStep(G4int i) const
    { return nStep[i]; }
};

extern G4Allocator<ExN07Run> anExN07RunAllocator;

inline void* ExN07Run::operator new(size_t)
{
  void* aRun;
  aRun = (void*)anExN07RunAllocator.MallocSingle();
  return aRun;
}

inline void ExN07Run::operator delete(void* aRun)
{
  anExN07RunAllocator.FreeSingle((ExN07Run*)aRun);
}

#endif

