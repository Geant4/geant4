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
// $Id: DeallocatorTest.cc,v 1.2 2004-05-26 14:44:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4Timer.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Type
{
  public:

    G4Type(){}
    ~G4Type(){}

    inline void *operator new(size_t);
    inline void operator delete(void *aObj);

  private:
  
    G4int         i1, i2, i3, i4, i5, i6, i7;
    G4double      d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13;
    G4double*     p1, p2, p3, p4, p5, p6, p7;
    G4ThreeVector v1, v2, v3;
};

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

int main()
{
  G4Timer timer;
  G4Type* pObj = 0;
  const size_t maxiter = 10000000;
  const size_t modulo  = 1000000;

  timer.Start();
  for (size_t i=0; i<maxiter; i++)
  {
    pObj = new G4Type();
    if (i%modulo == 0)
    {
      G4cout << "Clearing storage ..." << G4endl
             << "   allocator size before: "
             << aAllocator.GetAllocatedSize()
             << " bytes" << G4endl;
      aAllocator.ResetStorage();
      G4cout << "   allocator size after : "
             << aAllocator.GetAllocatedSize()
             << " bytes" << G4endl;
    }
  }
  for (size_t j=0; j<maxiter; j++)
  {
    delete pObj;
  }
  timer.Stop();  

  G4Type ref;
  G4cout << "Object size: " << sizeof(ref) << G4endl
         << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time  : " << timer.GetUserElapsed() << G4endl;

  return 0;
}
