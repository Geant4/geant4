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
// $Id: G4HepRepFileSceneHandler.cc,v 1.3 2001-11-19 15:07:28 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4HepRepFileGraphicsScene (was in modeling).

#include "G4HepRepFileSceneHandler.hh"
#include "G4HepRepFile.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ModelingParameters.hh"
#include "G4Polymarker.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"

//HepRep
#include "HepRepXMLWriter.hh"

G4int G4HepRepFileSceneHandler::fSceneIdCount = 0;
// Counter for HepRep scene handlers.

G4int G4HepRepFileSceneHandler::fSceneCount = 0;
// No. of extanct scene handlers.

G4HepRepFileSceneHandler::G4HepRepFileSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0)
{
  fSceneCount++;

  hepRepXMLWriter = ((G4HepRepFile*)(&system))->GetHepRepXMLWriter();
  G4cout << "  From scene handler const: opening file, G4HepRepFile.xml"
	 << G4endl;
}

G4HepRepFileSceneHandler::~G4HepRepFileSceneHandler() {}

#ifdef G4HEPREPFILEDEBUG
void G4HepRepFileSceneHandler::PrintThings() {
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
#endif

void G4HepRepFileSceneHandler::AddThis(const G4Box& box) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(box);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Cons& cons) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Cons& cons) called for "
	 << cons.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(cons);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Tubs& tubs) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Tubs& tubs) called for "
	 << tubs.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(tubs);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Trd& trd) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Trd& trd) called for "
	 << trd.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(trd);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Trap& trap) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Trap& trap) called for "
	 << trap.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(trap);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Sphere& sphere) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Sphere& sphere) called for "
	 << sphere.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(sphere);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Para& para) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Para& para) called for "
	 << para.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(para);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Torus& torus) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Torus& torus) called for "
	 << torus.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(torus);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Polycone& polycone) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Polycone& polycone) called for "
	 << polycone.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(polycone);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4Polyhedra& polyhedra) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Polyhedra& polyhedra) called for "
	 << polyhedra.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(polyhedra);  // Invoke default action.
}

void G4HepRepFileSceneHandler::AddThis(const G4VSolid& solid) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddThis(const G4Solid& solid) called for "
	 << solid.GetName()
	 << G4endl;
  PrintThings();
#endif
  hepRepXMLWriter->addType(fpCurrentPV->GetName());
  G4VSceneHandler::AddThis(solid);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	 << G4endl;
  PrintThings();
#endif

  hepRepXMLWriter->addInstance();
  hepRepXMLWriter->addAttValue("DrawAs","Line");

  G4Colour color = GetColour(polyline);
  hepRepXMLWriter->addAttValue("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue());

  hepRepXMLWriter->addPrimitive();

  for (size_t i=0; i < polyline.size(); i++) {
      G4Point3D vertex = (*fpObjectTransformation) * polyline[i];
      hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
  }
}


void G4HepRepFileSceneHandler::AddPrimitive (const G4Polymarker& line) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polymarker& line) called"
	 << G4endl;
  PrintThings();
#endif

  hepRepXMLWriter->addInstance();
  hepRepXMLWriter->addAttValue("DrawAs","Point");

  // FIXME: should be taken from G4
  hepRepXMLWriter->addAttValue("MarkName", "Square");
  hepRepXMLWriter->addAttValue("MarkSize", 5);

  G4Colour color = GetColour(line);
  hepRepXMLWriter->addAttValue("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue());

  hepRepXMLWriter->addPrimitive();
    
  for (size_t i=0; i < line.size(); i++) {
     G4Point3D vertex = (*fpObjectTransformation) * line[i];
     hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
  }
}

void G4HepRepFileSceneHandler::AddPrimitive(const G4Text& text) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	 << G4endl;
  PrintThings();
#endif
  G4cout << "G4HepRepFileSceneHandler::AddPrimitive G4Text : not yet implemented. " << G4endl;
}

void G4HepRepFileSceneHandler::AddPrimitive(const G4Circle& circle) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Circle& circle) called:"
    "\n  radius: " << circle.GetWorldRadius()
	 << G4endl;
  PrintThings();
#endif

  hepRepXMLWriter->addInstance();
  hepRepXMLWriter->addAttValue("DrawAs","Point");
  hepRepXMLWriter->addAttValue("MarkName", "Dot");

  G4Colour color = GetColour(circle);
  hepRepXMLWriter->addAttValue("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue());

  hepRepXMLWriter->addPrimitive();

  G4double  radius = circle.GetWorldSize();
  hepRepXMLWriter->addAttValue("MarkSize", radius);

  G4Point3D center = (*fpObjectTransformation) * circle.GetPosition();
  hepRepXMLWriter->addPoint(center.x(), center.y(), center.z());
}

void G4HepRepFileSceneHandler::AddPrimitive(const G4Square& square) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Square& square) called:"
    "\n  side: " << square.GetWorldRadius()
	 << G4endl;
  PrintThings();
#endif

  hepRepXMLWriter->addInstance();
  hepRepXMLWriter->addAttValue("DrawAs","Point");
  hepRepXMLWriter->addAttValue("MarkName", "Square");

  G4Colour color = GetColour(square);
  hepRepXMLWriter->addAttValue("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue());

  hepRepXMLWriter->addPrimitive();

  G4double size = square.GetWorldSize();
  hepRepXMLWriter->addAttValue("MarkSize", size);

  G4Point3D center = (*fpObjectTransformation) * square.GetPosition();
  hepRepXMLWriter->addPoint(center.x(), center.y(), center.z());
}

void G4HepRepFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	 << G4endl;
  PrintThings();
#endif
  if(polyhedron.GetNoFacets()==0)return;

  G4Normal3D surfaceNormal;
  G4Point3D vertex;

  hepRepXMLWriter->addInstance();
  hepRepXMLWriter->addAttValue("DrawAs","Polygon");

  G4Colour color = GetColour(polyhedron);
  hepRepXMLWriter->addAttValue("LineColor", color.GetRed(), color.GetGreen(), color.GetBlue());

  // Additional attributes.
  hepRepXMLWriter->addAttValue("LVol", fpCurrentLV->GetName());
  hepRepXMLWriter->addAttValue("Solid", fpCurrentLV->GetSolid()->GetName());
  hepRepXMLWriter->addAttValue("EType", fpCurrentLV->GetSolid()->GetEntityType());
  hepRepXMLWriter->addAttValue("Material", fpCurrentLV->GetMaterial()->GetName());
  hepRepXMLWriter->addAttValue("Density", fpCurrentLV->GetMaterial()->GetDensity());
  hepRepXMLWriter->addAttValue("State", fpCurrentLV->GetMaterial()->GetState());
  hepRepXMLWriter->addAttValue("Radlen", fpCurrentLV->GetMaterial()->GetRadlen());

  G4bool notLastFace;
  do {
      hepRepXMLWriter->addPrimitive();
      notLastFace = polyhedron.GetNextNormal (surfaceNormal);

      G4int edgeFlag = 1;
      G4bool notLastEdge;
      do {
          notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
          vertex = (*fpObjectTransformation) * vertex;
	  hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
      } while (notLastEdge);
  } while (notLastFace);
}

void G4HepRepFileSceneHandler::AddPrimitive(const G4NURBS& nurbs) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4NURBS& nurbs) called."
	 << G4endl;
  PrintThings();
#endif
    G4cout << "G4HepRepFileSceneHandler::AddPrimitive G4NURBS : not implemented. " << G4endl;
}

void G4HepRepFileSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace(&fCurrentDepth,
				       &fpCurrentPV,
				       &fpCurrentLV);
}

HepRepXMLWriter *G4HepRepFileSceneHandler::GetHepRepXMLWriter() {
    return hepRepXMLWriter;
}
