// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHit.hh,v 2.1 1998/07/12 02:53:42 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4VHit_h
#define G4VHit_h 1

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

