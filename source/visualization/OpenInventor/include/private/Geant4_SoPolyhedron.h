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
#ifndef Geant4_SoPolyhedron_h
#define Geant4_SoPolyhedron_h

// Inheritance :
#include <Inventor/nodes/SoShape.h>

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFNode.h>

class G4Polyhedron;

/**
 *  SoPolyhedron is an Inventor encapsulation of the G4Polyedron
 * class written by E.Chernyaev. In particular SoPolyhedron permits
 * to represent boolean operations over solids.
 *  To avoid clashes with other libraries (Geant4) where the G4Polyhedron
 * classes may be found, the G4Polyhedron (through usage of CPP macros)
 * had been renamed G4Polyhedron in the HEPVis lib.
 *  The solids are modeled with G4Polyedron objects. The G4Polyhedron
 * permits to produce a new G4Polyhedron according to the boolean 
 * operation done on them. The resulting G4Polyhedron  is then given
 * to an SoPolyhedron for rendering.
 *  Note that a boolean operation could be rendered in wire frame
 * by drawing the contour of the resulting solid (not by drawing
 * the wire frame of a triangulation).
 *  See the applications/Polyhedron example.
 */

class Geant4_SoPolyhedron : public SoShape {
  SO_NODE_HEADER(Geant4_SoPolyhedron);
private:
  Geant4_SoPolyhedron(const Geant4_SoPolyhedron&);
  Geant4_SoPolyhedron& operator=(const Geant4_SoPolyhedron&);
public:
  SoSFBool solid;
  SoSFBool reducedWireFrame;
  SoSFNode alternateRep;
public:
  Geant4_SoPolyhedron();
  Geant4_SoPolyhedron(const G4Polyhedron&);
  Geant4_SoPolyhedron(G4Polyhedron*);
  virtual void generateAlternateRep();
  virtual void clearAlternateRep();
public:
  static void initClass();
protected:
  virtual void computeBBox(SoAction*,SbBox3f&,SbVec3f&);
  virtual void generatePrimitives(SoAction*);
  virtual void doAction(SoAction*);
protected:
  virtual ~Geant4_SoPolyhedron();
private: 
  G4Polyhedron* fPolyhedron;
};

#endif
