// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigi.hh,v 1.4.4.1 1999/12/07 20:47:46 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef G4VDigi_h
#define G4VDigi_h 1

// class description:
//
//  This is the base class of digi object. The user should derive this
// base class to make his/her own digi class. Two virtual method Draw()
// and Print() can be implemented if the user wants these functionarities.
//  If a concrete digi class is used as a transient class, G4Allocator
// must be used.

class G4VDigi 
{

  public:
      G4VDigi();
      virtual ~G4VDigi();

      int operator==(const G4VDigi &right) const;

      virtual void Draw();
      virtual void Print();

};

#endif

