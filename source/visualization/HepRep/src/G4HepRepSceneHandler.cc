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
// $Id: G4HepRepSceneHandler.cc,v 1.1 2001-08-24 23:27:34 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4HepRepGraphicsScene (was in modeling).

#include "G4HepRepSceneHandler.hh"
#include "G4HepRep.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ModelingParameters.hh"
#include "G4Polymarker.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"

//HepRep
#include "JHepRepAttribute.hh"
#include "JHepRepFactory.hh"
#include "JHepRep.hh"
#include "JHepRepType.hh"
#include "JHepRepInstance.hh"
#include "JHepRepPrimitive.hh"
#include "JHepRepPoint.hh"

G4int G4HepRepSceneHandler::fSceneIdCount = 0;
// Counter for HepRep scene handlers.

G4int G4HepRepSceneHandler::fSceneCount = 0;
// No. of extanct scene handlers.

G4HepRepSceneHandler::G4HepRepSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0)
{
  fSceneCount++;

  factory = ((G4HepRep*)(&system))->GetHepRepFactory();
  heprep = ((G4HepRep*)(&system))->GetHepRep();
  type = factory->CreateHepRepType(heprep, (const char *)name, "Geant 4 (HepRep)");
}

G4HepRepSceneHandler::~G4HepRepSceneHandler() {}

void G4HepRepSceneHandler::PrintThings() {
  G4cout <<
    "  with transformation "
	 << (void*)fpObjectTransformation;
  if (fpCurrentPV) {
    G4cout <<
      "\n  current physical volume: "
	   << fpCurrentPV->GetName() <<
      "\n  current logical volume: "
	   << fpCurrentLV->GetName() <<
      "\n  current depth of geometry tree: "
	   << fCurrentDepth;
  }
  G4cout << G4endl;
}

void G4HepRepSceneHandler::AddThis(const G4Box& box) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(box);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Cons& cons) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Cons& cons) called for "
	 << cons.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(cons);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Tubs& tubs) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Tubs& tubs) called for "
	 << tubs.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(tubs);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Trd& trd) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Trd& trd) called for "
	 << trd.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(trd);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Trap& trap) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Trap& trap) called for "
	 << trap.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(trap);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Sphere& sphere) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Sphere& sphere) called for "
	 << sphere.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(sphere);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Para& para) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Para& para) called for "
	 << para.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(para);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Torus& torus) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Torus& torus) called for "
	 << torus.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(torus);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Polycone& polycone) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Polycone& polycone) called for "
	 << polycone.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(polycone);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4Polyhedra& polyhedra) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Polyhedra& polyhedra) called for "
	 << polyhedra.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(polyhedra);  // Invoke default action.
}

void G4HepRepSceneHandler::AddThis(const G4VSolid& solid) {
  G4cout <<
    "G4HepRepSceneHandler::AddThis(const G4Solid& solid) called for "
	 << solid.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(solid);  // Invoke default action.
}


void G4HepRepSceneHandler::AddPrimitive(const G4Polyline& polyline) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	 << G4endl;
  PrintThings();

  instance = factory->CreateHepRepInstance(type);
  instance->AddValue("DrawAs", "Line");
  primitive = factory->CreateHepRepPrimitive(instance);
  SetColour(instance, GetColour(polyline));
  for (size_t i=0; i < polyline.size(); i++) {
      G4Point3D vertex = (*fpObjectTransformation) * polyline[i];
      factory->CreateHepRepPoint(primitive, vertex.x(), vertex.y(), vertex.z());
  }
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polymarker& line) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Polymarker& line) called"
	 << G4endl;
  PrintThings();

  instance = factory->CreateHepRepInstance(type);
  instance->AddValue("DrawAs", "Point");
    
  // FIXME: should be taken from G4
  instance->AddValue("MarkName", "Square");
  instance->AddValue("MarkSize", 5);
    
  primitive = factory->CreateHepRepPrimitive(instance);
  SetColour(instance, GetColour(line));
  for (size_t i=0; i < line.size(); i++) {
     G4Point3D vertex = (*fpObjectTransformation) * line[i];
     factory->CreateHepRepPoint(primitive, vertex.x(), vertex.y(), vertex.z());
  }
}

void G4HepRepSceneHandler::AddPrimitive(const G4Text& text) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	 << G4endl;
  PrintThings();
    G4cout << "G4HepRepSceneHandler::AddPrimitive G4Text : not yet implemented. " << G4endl;
}

void G4HepRepSceneHandler::AddPrimitive(const G4Circle& circle) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Circle& circle) called:"
    "\n  radius: " << circle.GetWorldRadius()
	 << G4endl;
  PrintThings();

  G4Point3D center = (*fpObjectTransformation) * circle.GetPosition();
  G4double  radius = circle.GetWorldSize();

  instance = factory->CreateHepRepInstance(type);
  SetColour (instance, GetColour(circle));
  instance->AddValue("DrawAs", "Point");
  instance->AddValue("MarkName", "Dot");
  instance->AddValue("MarkSize", radius);
    
  primitive = factory->CreateHepRepPrimitive(instance);
  factory->CreateHepRepPoint(primitive, center.x(), center.y(), center.z());
}

void G4HepRepSceneHandler::AddPrimitive(const G4Square& square) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Square& square) called:"
    "\n  side: " << square.GetWorldRadius()
	 << G4endl;
  PrintThings();

  G4Point3D center = (*fpObjectTransformation) * square.GetPosition();
  G4double size = square.GetWorldSize();
  instance = factory->CreateHepRepInstance(type);
  SetColour (instance, GetColour(square));
  instance->AddValue("DrawAs", "Point");
  instance->AddValue("MarkName", "Square");
  instance->AddValue("MarkSize", size);
    
  primitive = factory->CreateHepRepPrimitive(instance);
  factory->CreateHepRepPoint(primitive, center.x(), center.y(), center.z());
}

void G4HepRepSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	 << G4endl;
  PrintThings();

  G4Normal3D surfaceNormal;
  G4Point3D vertex;

  instance = factory->CreateHepRepInstance(type);
  instance->AddValue("DrawAs", "Polygon");
  SetColour(instance, GetColour(polyhedron));
    
  G4bool notLastFace;
  do {
      primitive = factory->CreateHepRepPrimitive(instance);
      notLastFace = polyhedron.GetNextNormal (surfaceNormal);

      G4int edgeFlag = 1;
      G4bool notLastEdge;
      do {
          notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
          vertex = (*fpObjectTransformation) * vertex;
//            G4cout << "Vertex: (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ")"
//                   << notLastEdge << ", " << edgeFlag << G4endl;
          factory->CreateHepRepPoint(primitive, vertex.x(), vertex.y(), vertex.z());
      } while (notLastEdge);

  } while (notLastFace);
}

void G4HepRepSceneHandler::AddPrimitive(const G4NURBS& nurbs) {
  G4cout <<
    "G4HepRepSceneHandler::AddPrimitive(const G4NURBS& nurbs) called."
	 << G4endl;
  PrintThings();
}

void G4HepRepSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace(&fCurrentDepth,
				       &fpCurrentPV,
				       &fpCurrentLV);
}

void G4HepRepSceneHandler::SetColour (JHepRepAttribute *attribute, const G4Colour& color) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::SetColour : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<  
                                  " blue : " << color.GetBlue ()   << G4endl;
#endif
    attribute->AddColor("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue(), color.GetAlpha());
}

JHepRep *G4HepRepSceneHandler::GetHepRep() {
    return heprep;
}

JHepRepFactory *G4HepRepSceneHandler::GetHepRepFactory() {
    return factory;
}
