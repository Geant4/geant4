// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LineCreator.hh,v 1.2 2000-01-21 13:45:27 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4LineCreator
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4LINECREATOR_HH
#define G4LINECREATOR_HH

#include "G4GeometryCreator.hh"

class G4LineCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4LineCreator();
    ~G4LineCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    G4String Name() { return "Line"; }
    static G4LineCreator GetInstance() { return csc; }

  // Members

  private:

    static G4LineCreator csc;
};

#endif
