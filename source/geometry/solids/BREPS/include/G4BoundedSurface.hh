// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BoundedSurface.hh,v 1.1 1999-01-07 16:07:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Surface.hh"

class G4BoundedSurface: public G4Surface
{
public:
  G4BoundedSurface() {};

/* L. Broglia
  G4BoundedSurface(STEPentity& Ent, InstMgr&) {};  
*/

  ~G4BoundedSurface() {};  

  virtual char *Name() const
  {
    return "G4BoundedSurface"; 
  }  

};
