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
// $Id: Mars01Hit.hh,v 1.1 2001-12-13 14:58:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Mars01Hit_h
#define Mars01Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class Mars01Hit : public G4VHit
{
public:
  
  Mars01Hit();
  Mars01Hit(G4int id);
  ~Mars01Hit();
  Mars01Hit(const Mars01Hit &right);
  const Mars01Hit& operator=(const Mars01Hit &right);
  int operator==(const Mars01Hit &right) const;
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
  void operator delete(void *aHit,void*){}
#endif
  
  void Draw();
  void Print();
  
private:
  G4double edep;
  G4double npar;
  G4int id;
  
public:
  // tally : deposit energy
  inline void SetEdep(G4double de)
  { edep = de; };
  inline void AddEdep(G4double de)
  { edep += de; };
  inline G4double GetEdep()
  { return edep; };
  // tally : number of particle
  inline void SetNpar(G4double np)
  { npar = np; };
  inline void AddNpar(G4double np)
  { npar += np; };
  inline G4double GetNpar()
  { return npar; };
  // get ID
  inline G4double GetID()
  { return id; };

};

typedef G4THitsCollection<Mars01Hit> Mars01HitsCollection;

extern G4Allocator<Mars01Hit> Mars01HitAllocator;

inline void* Mars01Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Mars01HitAllocator.MallocSingle();
  return aHit;
}

inline void Mars01Hit::operator delete(void *aHit)
{
  Mars01HitAllocator.FreeSingle((Mars01Hit*) aHit);
}

#endif



