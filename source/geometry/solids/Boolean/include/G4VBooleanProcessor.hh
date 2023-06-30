//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4VBooleanProcessor
//
// Class description:
//
// Virtual base class for Boolean solid processing.

// 20.03.2023 Stewart Boogert - created
// --------------------------------------------------------------------
#ifndef G4VBOOLEANPROCESSOR_HH
#define G4VBOOLEANPROCESSOR_HH

class G4Polyhedron;
class G4PolyhedronArbitrary;

class G4VBooleanProcessor
{
  public:

    G4VBooleanProcessor() = default;
    virtual ~G4VBooleanProcessor() = default;

    virtual G4PolyhedronArbitrary* Intersection(G4Polyhedron* /*p1*/,
                                                G4Polyhedron* /*p2*/)
      { return nullptr; }
    virtual G4PolyhedronArbitrary* Union(G4Polyhedron* /*p1*/,
                                         G4Polyhedron* /*p2*/)
      { return nullptr; }
    virtual G4PolyhedronArbitrary* Subtraction(G4Polyhedron* /*p1*/,
                                               G4Polyhedron* /*p2*/)
      { return nullptr; }
};

#endif
