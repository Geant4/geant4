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
// $Id: G4AssemblyCreator.hh,v 1.5 2002-11-21 16:49:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4AssemblyCreator
//
// Class description:
//
// Main class for detector geometry assembly read from a STEP file.
// Currently the default STEP reader is assumed to be the NIST reader
// from the SCL NIST public domain library.

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
  public:  // with description

  // Constructors
  
    G4AssemblyCreator();
    G4AssemblyCreator(const G4String&, const G4String& reader="NIST");

  // Destructor
  
    ~G4AssemblyCreator();    

  // Member functions
  
    void ReadStepFile();
    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    const char* Name() const { return "Assembly"; }
    static G4AssemblyCreator GetInstance();
    
  private:

    G4AssemblyCreator(const G4AssemblyCreator&);
    G4AssemblyCreator& operator=(const G4AssemblyCreator&);
      // Private copy constructor and assignment operator.

  private:

    G4int index;  // memory for get STEP entity - L.Broglia
    G4StepFileReader* StepReader;
    G4String STEPfileName;
    static G4AssemblyCreator ci;
};

#endif
