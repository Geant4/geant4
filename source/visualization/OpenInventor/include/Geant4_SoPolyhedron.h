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
