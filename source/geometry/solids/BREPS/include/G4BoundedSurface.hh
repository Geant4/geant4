// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BoundedSurface.hh,v 2.2 1998/10/20 16:31:12 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
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
