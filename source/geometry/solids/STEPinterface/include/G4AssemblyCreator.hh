// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AssemblyCreator.hh,v 1.2 2000-01-21 13:45:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4AssemblyCreator
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
#ifndef G4ASSEMBLYCREATOR_HH
#define G4ASSEMBLYCREATOR_HH

#include "G4GeometryCreator.hh"

class G4StepFileReader;

class G4AssemblyCreator: public G4GeometryCreator
{
  public:

  // Constructors
  
    G4AssemblyCreator();
    G4AssemblyCreator(G4String, G4String ="NIST");

  // Destructor
  
    ~G4AssemblyCreator();    

  // Member functions
  
    void ReadStepFile();
    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    G4String Name() { return "Assembly"; }
    static G4AssemblyCreator GetInstance() { return ci; }
  
  // Members
  
    int index; // add by L. Broglia, memory for get STEP entity
  
  private:
  
    G4StepFileReader* StepReader;
    G4String STEPfileName;
    static G4AssemblyCreator ci;
};

#endif
