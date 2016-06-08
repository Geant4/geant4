// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHit.hh,v 1.4 1999/12/15 14:49:37 gunter Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#ifndef G4VHit_h
#define G4VHit_h 1

// class description:
//
//  This is the base class of hit object. The user should derive this
// base class to make his/her own hit class. Two virtual method Draw()
// and Print() can be implemented if the user wants these functionarities.
//  If a concrete hit class is used as a transient class, G4Allocator
// must be used.

class G4VHit 
{

  public:
      G4VHit();
      virtual ~G4VHit();

      int operator==(const G4VHit &right) const;

      virtual void Draw();
      virtual void Print();

};

#endif

