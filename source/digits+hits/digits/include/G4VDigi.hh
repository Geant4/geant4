// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigi.hh,v 1.1 1999/01/07 16:06:28 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4VDigi_h
#define G4VDigi_h 1

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

