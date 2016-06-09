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
#ifndef Geant4_SoPolyhedron_h
#define Geant4_SoPolyhedron_h

// Inheritance :
#include <Inventor/nodes/SoShape.h>

#include <Inventor/fields/SoSFBool.h>

class HepPolyhedron;

/**
 *  Geant4_SoPolyhedron is an Inventor encapsulation of the HepPolyedron
 * class written by E.Chernyaev. In particular Geant4_SoPolyhedron permits
 * to represent boolean operations over solids.
 *  The solids are modeled with HepPolyedron objects. The HepPolyhedron
 * permits to produce a new HepPolyhedron according to the boolean 
 * operation done on them. The resulting HepPolyhedron  is then given
 * to an Geant4_SoPolyhedron for rendering.
 *  Note that a boolean operation could be rendered in wire frame
 * by drawing the contour of the resulting solid (not by drawing
 * the wire frame of a triangulation).
 */

class Geant4_SoPolyhedron : public SoShape {
  SO_NODE_HEADER(Geant4_SoPolyhedron);
public:
  SoSFBool solid;
  SoSFBool reducedWireFrame;
public:
  Geant4_SoPolyhedron();
  Geant4_SoPolyhedron(const HepPolyhedron&);
  Geant4_SoPolyhedron(HepPolyhedron*);
public:
  static void initClass();
protected:
  virtual void computeBBox(SoAction*,SbBox3f&,SbVec3f&);
  virtual void generatePrimitives(SoAction*);
protected:
  virtual ~Geant4_SoPolyhedron();
private: 
  HepPolyhedron* fPolyhedron;
};

#endif
