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
// $Id: G4GeometryTable.hh,v 1.7 2002-11-21 16:49:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometryTable
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
#ifndef G4GEOMETRYTABLE_HH
#define G4GEOMETRYTABLE_HH

#include "g4std/vector"
#include "globals.hh"
#include "G4GeometryCreator.hh"

typedef G4std::vector<G4GeometryCreator*> G4CreatorVector;

class G4GeometryTable
{
  public:

  // Destructor
  
    ~G4GeometryTable();    

  // Member functions

    static void RegisterObject(G4GeometryCreator*);
    static G4bool ExistsInTable(G4String&);
    static G4GeometryCreator* GetObject(G4String);
    static void* CreateObject(STEPentity&);
    static void* CreateSTEPObject(void*, G4String&);  
    static void PrintObjectNames();
    static const G4GeometryTable& GetInstance();

  public:

  // Constructor, singleton

    G4GeometryTable();

  // Members

  private:

    static G4CreatorVector RegisteredObjects;
    static G4GeometryTable gt;
};

#endif
