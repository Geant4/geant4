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
// $Id: AllocatorTest.cc,v 1.1 2003-04-07 13:00:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4Timer.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ios.hh"

#define USE_G4ALLOCATOR 1

class G4Type
{
  public:

    G4Type(){}
    ~G4Type(){}

#ifdef USE_G4ALLOCATOR
    inline void *operator new(size_t);
    inline void operator delete(void *aObj);
#endif

  private:
  
    G4int         i1, i2, i3, i4, i5, i6, i7;
    G4double      d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13;
    G4double*     p1, p2, p3, p4, p5, p6, p7;
    G4ThreeVector v1, v2, v3;
};

#ifdef USE_G4ALLOCATOR

  G4Allocator<G4Type> aAllocator;

  inline void* G4Type::operator new(size_t)
  {
    void *aValue;
    aValue = (void *) aAllocator.MallocSingle();
    return aValue;
  }

  inline void G4Type::operator delete(void *aValue)
  {
    aAllocator.FreeSingle((G4Type *) aValue);
  }

#endif

int main()
{
  G4Timer timer;
  G4Type* pObj = 0;
  G4Type* pDum = 0;

  timer.Start();
  for (size_t i=0; i<5000000; i++)
  {
    pObj = new G4Type();
    delete pDum;
    delete pObj;
  }
  timer.Stop();  

  G4Type ref;
  G4cout << "Object size: " << sizeof(ref) << G4endl
         << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time  : " << timer.GetUserElapsed() << G4endl;

  return 0;
}
