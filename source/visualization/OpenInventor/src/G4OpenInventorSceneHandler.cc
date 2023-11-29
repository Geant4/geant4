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
//
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor stored scene - creates OpenInventor display lists.

// this :
#include "G4OpenInventorSceneHandler.hh"

#include <Inventor/SoPath.h>
#include <Inventor/SoNodeKitPath.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoCoordinate4.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoResetTransform.h>
#include <Inventor/nodes/SoMatrixTransform.h>

#define USE_SOPOLYHEDRON

#ifndef USE_SOPOLYHEDRON
#include "HEPVis/nodes/SoBox.h"
#include "HEPVis/nodes/SoTubs.h"
#include "HEPVis/nodes/SoCons.h"
#include "HEPVis/nodes/SoTrd.h"
#include "HEPVis/nodes/SoTrap.h"
#endif
#include "HEPVis/nodes/SoMarkerSet.h"

using SoMarkerSet = HEPVis_SoMarkerSet;

#include "HEPVis/nodekits/SoDetectorTreeKit.h"
#include "HEPVis/misc/SoStyleCache.h"

#include "SoG4Polyhedron.h"
#include "SoG4LineSet.h"
#include "SoG4MarkerSet.h"

#include "G4Scene.hh"
#include "G4OpenInventor.hh"
#include "G4OpenInventorTransform3D.hh"
#include "G4ThreeVector.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4ModelingParameters.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#define G4warn G4cout

G4int G4OpenInventorSceneHandler::fSceneIdCount = 0;

G4OpenInventorSceneHandler::G4OpenInventorSceneHandler (G4OpenInventor& system,
                                          const G4String& name)
:G4VSceneHandler (system, fSceneIdCount++, name)
,fRoot(0)
,fDetectorRoot(0)
,fTransientRoot(0)
,fCurrentSeparator(0)
,fModelingSolid(false)
,fReducedWireFrame(true)
,fStyleCache(0)
,fPreviewAndFull(true)
{
  fStyleCache = new SoStyleCache;
  fStyleCache->ref();

  fRoot = new SoSeparator;
  fRoot->ref();
  fRoot->setName("Root");
  
  fDetectorRoot = new SoSeparator;
  fDetectorRoot->setName("StaticRoot");
  fRoot->addChild(fDetectorRoot);
  
  fTransientRoot = new SoSeparator;
  fTransientRoot->setName("TransientRoot");
  fRoot->addChild(fTransientRoot);
  
  fCurrentSeparator = fTransientRoot;
}

G4OpenInventorSceneHandler::~G4OpenInventorSceneHandler ()
{
  fRoot->unref();
  fStyleCache->unref();
}

void G4OpenInventorSceneHandler::ClearStore ()
{
  fDetectorRoot->removeAllChildren();
  fSeparatorMap.clear();

  fTransientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::ClearTransientStore ()
{
  fTransientRoot->removeAllChildren();
}

//
// Generates prerequisites for solids
//  
void G4OpenInventorSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation,
 const G4VisAttributes& visAttribs)
{
  G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);
  // Stores arguments away for future use, e.g., AddPrimitives.

  GeneratePrerequisites();
}

//
// Generates prerequisites for primitives
//  
void G4OpenInventorSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {

  G4VSceneHandler::BeginPrimitives (objectTransformation);

  // If thread of control has already passed through PreAddSolid,
  // avoid opening a graphical data base component again.
  if (!fProcessingSolid) {
    GeneratePrerequisites();
  }
}

//
// Method for handling G4Polyline objects (from tracking).
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyline& line)
{
  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4OpenInventorSceneHandler::AddPrimitive (const G4Polyline&)",
	 "OpenInventor-0001", JustWarning,
	 "2D polylines not implemented.  Ignored.");
    }
    return;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
  fpViewer -> GetApplicableVisAttributes (line.GetVisAttributes ());
  
  AddProperties(pVA);  // Colour, etc.
  AddTransform();      // Transformation

  G4int nPoints = (G4int)line.size();
  SbVec3f* pCoords = new SbVec3f[nPoints];

  for (G4int iPoint = 0; iPoint < nPoints ; ++iPoint) {
    pCoords[iPoint].setValue((G4float)line[iPoint].x(),
                             (G4float)line[iPoint].y(),
                             (G4float)line[iPoint].z());
  }

  //
  // Point Set
  // 
  SoCoordinate3 *polyCoords = new SoCoordinate3;
  polyCoords->point.setValues(0,nPoints,pCoords);
  fCurrentSeparator->addChild(polyCoords);
  
  //
  // Wireframe
  // 
  SoDrawStyle* drawStyle = fStyleCache->getLineStyle();
  fCurrentSeparator->addChild(drawStyle);

  SoG4LineSet *pLine = new SoG4LineSet;

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) LoadAtts(line, pLine);

#ifdef INVENTOR2_0
  pLine->numVertices.setValues(0,1,(const G4long *)&nPoints);
#else 
  pLine->numVertices.setValues(0,1,&nPoints);
#endif

  fCurrentSeparator->addChild(pLine);

  delete [] pCoords;
}

void G4OpenInventorSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4OpenInventorSceneHandler::AddPrimitive (const G4Polymarker&)",
	 "OpenInventor-0002", JustWarning,
	 "2D polymarkers not implemented.  Ignored.");
    }
    return;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
  fpViewer -> GetApplicableVisAttributes (polymarker.GetVisAttributes ());

  AddProperties(pVA);  // Colour, etc.
  AddTransform();      // Transformation

  G4int pointn = (G4int)polymarker.size();
  if(pointn<=0) return;

  SbVec3f* points = new SbVec3f[pointn];
  for (G4int iPoint = 0; iPoint < pointn ; ++iPoint) {
    points[iPoint].setValue((G4float)polymarker[iPoint].x(),
                            (G4float)polymarker[iPoint].y(),
                            (G4float)polymarker[iPoint].z());
  }

  SoCoordinate3* coordinate3 = new SoCoordinate3;
  coordinate3->point.setValues(0,pointn,points);
  fCurrentSeparator->addChild(coordinate3);

  MarkerSizeType sizeType;
  G4double screenSize = GetMarkerSize (polymarker, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.  OK.
    break;
  case world:
    // Draw in world coordinates.   Not implemented.  Use screenSize = 10.
    screenSize = 10.;
    break;
  }
  
  SoG4MarkerSet* markerSet = new SoG4MarkerSet;
  markerSet->numPoints = pointn;

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking())
    LoadAtts(polymarker, markerSet);

  G4VMarker::FillStyle style = polymarker.GetFillStyle();
  switch (polymarker.GetMarkerType()) {
  default:
    // Are available 5_5, 7_7 and 9_9
  case G4Polymarker::dots:
    if (screenSize <= 5.) {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
    } else if (screenSize <= 7.) {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
    } else {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
    }
    break;
  case G4Polymarker::circles:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_9_9;
      }
    }
    break;
  case G4Polymarker::squares:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_9_9;
      }
    }
  }
  fCurrentSeparator->addChild(markerSet);

  delete [] points;
}

// Method for handling G4Text objects
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Text& text)
{
  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4OpenInventorSceneHandler::AddPrimitive (const G4Text&)",
	 "OpenInventor-0003", JustWarning,
	 "2D text not implemented.  Ignored.");
    }
    return;
  }

  AddProperties(text.GetVisAttributes());  // Colour, etc.
  AddTransform(text.GetPosition());        // Transformation

  //
  // Color.  Note: text colour is worked out differently.  This
  // over-rides the colour added in AddProperties...
  //
  const G4Colour& c = GetTextColour (text);
  SoMaterial* material = 
    fStyleCache->getMaterial((G4float)c.GetRed(),
                             (G4float)c.GetGreen(),
                             (G4float)c.GetBlue(),
                             (G4float)(1-c.GetAlpha()));
  fCurrentSeparator->addChild(material);

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (text, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.  OK.
    break;
  case world:
    // Draw in world coordinates.   Not implemented.  Use size = 20.
    size = 20.;
    break;
  }

  //
  // Font
  // 
  SoFont *g4Font = new SoFont();
  g4Font->size = size;
  fCurrentSeparator->addChild(g4Font);

  //
  // Text
  // 
  SoText2 *g4String = new SoText2();
  g4String->string.setValue(text.GetText());
  g4String->spacing = 2.0;
  switch (text.GetLayout()) {
  default:
  case G4Text::left:
    g4String->justification = SoText2::LEFT; break;
  case G4Text::centre:
    g4String->justification = SoText2::CENTER; break;
  case G4Text::right:
    g4String->justification = SoText2::RIGHT; break;
  }
  fCurrentSeparator->addChild(g4String);
}

//
// Method for handling G4Circle objects
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Circle& circle) {
  AddCircleSquare(G4OICircle, circle);
}

//
// Method for handling G4Square objects - defaults to wireframe
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Square& square) {
  AddCircleSquare(G4OISquare, square);
}

void G4OpenInventorSceneHandler::AddCircleSquare
(G4OIMarker markerType, const G4VMarker& marker)
{
  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4OpenInventorSceneHandler::AddCircleSquare",
	 "OpenInventor-0004", JustWarning,
	 "2D circles and squares not implemented.  Ignored.");
    }
    return;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
  fpViewer -> GetApplicableVisAttributes (marker.GetVisAttributes ());

  AddProperties(pVA);  // Colour, etc.
  AddTransform();      // Transformation

  MarkerSizeType sizeType;
  G4double screenSize = GetMarkerSize (marker, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.  OK.
    break;
  case world:
    // Draw in world coordinates.   Not implemented.  Use size = 10.
    screenSize = 10.;
    break;
  }

  G4Point3D centre = marker.GetPosition();

  // Borrowed from AddPrimitive(G4Polymarker) - inefficient? JA
  SbVec3f* points = new SbVec3f[1];
  points[0].setValue((G4float)centre.x(),
		     (G4float)centre.y(),
		     (G4float)centre.z());
  SoCoordinate3* coordinate3 = new SoCoordinate3;
  coordinate3->point.setValues(0,1,points);
  fCurrentSeparator->addChild(coordinate3);

  SoG4MarkerSet* markerSet = new SoG4MarkerSet;
  markerSet->numPoints = 1;

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) LoadAtts(marker, markerSet);

  G4VMarker::FillStyle style = marker.GetFillStyle();
  switch (markerType) {
  case G4OICircle:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_9_9;
      }
    }
    break;
  case G4OISquare:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_9_9;
      }
    }
  break;
  }
  fCurrentSeparator->addChild(markerSet);

  delete [] points;
}

//
// Method for handling G4Polyhedron objects for drawing solids.
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron)
{
  if (polyhedron.GetNoFacets() == 0) return;

  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4OpenInventorSceneHandler::AddPrimitive (const G4Polyhedron&)",
	 "OpenInventor-0005", JustWarning,
	 "2D polyhedra not implemented.  Ignored.");
    }
    return;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
  fpViewer -> GetApplicableVisAttributes (polyhedron.GetVisAttributes ());

  AddProperties(pVA);  // Colour, etc.
  AddTransform();      // Transformation

  SoG4Polyhedron* soPolyhedron = new SoG4Polyhedron(polyhedron);

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking())
    LoadAtts(polyhedron, soPolyhedron);

  SbString name = "Non-geometry";
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel) {
    name = pPVModel->GetCurrentLV()->GetName().c_str();
  }
  SbName sbName(name);
  soPolyhedron->setName(sbName);
  soPolyhedron->solid.setValue(fModelingSolid);
  soPolyhedron->reducedWireFrame.setValue(fReducedWireFrame?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
}

void G4OpenInventorSceneHandler::AddCompound(const G4Mesh& mesh) {
  StandardSpecialMeshRendering(mesh);
}

void G4OpenInventorSceneHandler::GeneratePrerequisites()
{
  // Utility for PreAddSolid and BeginPrimitives.

  // This routines prepares for adding items to the scene database.  We
  // are expecting two kinds of solids: leaf parts and non-leaf parts.
  // For non-leaf parts, we create a detector tree kit.  This has two
  // separators.  The solid itself goes in the preview separator, the
  // full separator is forseen for daughters.  This separator is not
  // only created--it is also put in a dictionary for future use by
  // the leaf part.

  // For leaf parts, we first locate the mother volume and find its
  // separator through the dictionary.

  // The private member fCurrentSeparator is set to the proper
  // location on in the scene database so that when the solid is
  // actually added (in addthis), it is put in the right place.

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  
  if (pPVModel) {

    // This call comes from a G4PhysicalVolumeModel.  drawnPVPath is
    // the path of the current drawn (non-culled) volume in terms of
    // drawn (non-culled) ancestors.  Each node is identified by a
    // PVNodeID object, which is a physical volume and copy number.  It
    // is a vector of PVNodeIDs corresponding to the geometry hierarchy
    // actually selected, i.e., not culled.
    using PVNodeID = G4PhysicalVolumeModel::G4PhysicalVolumeNodeID;
    using PVPath = std::vector<PVNodeID>;
    const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
    //G4int currentDepth = pPVModel->GetCurrentDepth();
    G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
    G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
    //G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();
    // Note: pCurrentMaterial may be zero (parallel world).

    // The simplest algorithm, used by the Open Inventor Driver
    // developers, is to rely on the fact the G4PhysicalVolumeModel
    // traverses the geometry hierarchy in an orderly manner.  The last
    // mother, if any, will be the node to which the volume should be
    // added.  So it is enough to keep a map of scene graph nodes keyed
    // on the volume path ID.  Actually, it is enough to use the logical
    // volume as the key.  (An alternative would be to keep the PVNodeID
    // in the tree and match the PVPath from the root down.)

    // Find mother.  ri points to mother, if any...
    PVPath::const_reverse_iterator ri;
    G4LogicalVolume* MotherVolume = 0;
    ri = ++drawnPVPath.rbegin();
    if (ri != drawnPVPath.rend()) {
      // This volume has a mother.
      MotherVolume = ri->GetPhysicalVolume()->GetLogicalVolume();
    }

    if (pCurrentLV->GetNoDaughters()!=0 ||
	pCurrentPV->IsReplicated()) {  //????Don't understand this???? JA
      // This block of code is executed for non-leaf parts:

      // Make the detector tree kit:
      SoDetectorTreeKit* detectorTreeKit = new SoDetectorTreeKit();  

      SoSeparator* previewSeparator   =  
	(SoSeparator*) detectorTreeKit->getPart("previewSeparator",TRUE);
      previewSeparator->renderCaching = SoSeparator::OFF;

      SoSeparator* fullSeparator =  
	(SoSeparator*) detectorTreeKit->getPart("fullSeparator",   TRUE);
      fullSeparator->renderCaching = SoSeparator::OFF;

      if(fPreviewAndFull) detectorTreeKit->setPreviewAndFull();
      else detectorTreeKit->setPreview(TRUE);

      // Colour, etc., for SoDetectorTreeKit.  Treated differently to
      // othere SoNodes(?).  Use fpVisAttribs stored away in
      // PreAddSolid...
      const G4VisAttributes* pApplicableVisAttribs =
	fpViewer->GetApplicableVisAttributes (fpVisAttribs);

      // First find the color attributes...
      const G4Colour& g4Col =  pApplicableVisAttribs->GetColour ();
      const G4double red = g4Col.GetRed ();
      const G4double green = g4Col.GetGreen ();
      const G4double blue = g4Col.GetBlue ();
      G4double transparency = 1 - g4Col.GetAlpha();

      // Drawing style...
      G4ViewParameters::DrawingStyle drawing_style =
	GetDrawingStyle(pApplicableVisAttribs);
      switch (drawing_style) {
      case (G4ViewParameters::wireframe):    
	fModelingSolid = false;
	break;
      case (G4ViewParameters::hlr):
      case (G4ViewParameters::hsr):
      case (G4ViewParameters::hlhsr):
	fModelingSolid = true;
	break;
      case (G4ViewParameters::cloud):
        fModelingSolid = false;
        break;
      }

      SoMaterial* material = 
	fStyleCache->getMaterial((G4float)red,
				 (G4float)green,
				 (G4float)blue,
				 (G4float)transparency);
      detectorTreeKit->setPart("appearance.material",material);

      SoLightModel* lightModel = 
	fModelingSolid ? fStyleCache->getLightModelPhong() : 
	fStyleCache->getLightModelBaseColor();
      detectorTreeKit->setPart("appearance.lightModel",lightModel);

      // Add the full separator to the dictionary; it is indexed by the 
      // address of the logical volume!
      fSeparatorMap[pCurrentLV] = fullSeparator;

      // Find out where to add this volume.
      // If no mother can be found, it goes under root.
      if (MotherVolume) {
	if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
	  //printf("debug : PreAddSolid : mother %s found in map\n",
	  //     MotherVolume->GetName().c_str());
	  fSeparatorMap[MotherVolume]->addChild(detectorTreeKit);
	} else {
	  // Mother not previously encountered.  Shouldn't happen, since
	  // G4PhysicalVolumeModel sends volumes as it encounters them,
	  // i.e., mothers before daughters, in its descent of the
	  // geometry tree.  Error!
	  G4warn <<
	    "ERROR: G4OpenInventorSceneHandler::GeneratePrerequisites: Mother "
		 << ri->GetPhysicalVolume()->GetName()
		 << ':' << ri->GetCopyNo()
		 << " not previously encountered."
	    "\nShouldn't happen!  Please report to visualization coordinator."
		 << G4endl;
	  // Continue anyway.  Add to root of scene graph tree...
	  //printf("debug : PreAddSolid : mother %s not found in map !!!\n",
	  //     MotherVolume->GetName().c_str());
	  fDetectorRoot->addChild(detectorTreeKit);
	}
      } else {
	//printf("debug : PreAddSolid : has no mother\n");
	fDetectorRoot->addChild(detectorTreeKit);
      }

      fCurrentSeparator = previewSeparator;

    } else {
      // This block of code is executed for leaf parts.

      if (MotherVolume) {
	if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
	  fCurrentSeparator = fSeparatorMap[MotherVolume];
	} else {
	  // Mother not previously encountered.  Shouldn't happen, since
	  // G4PhysicalVolumeModel sends volumes as it encounters them,
	  // i.e., mothers before daughters, in its descent of the
	  // geometry tree.  Error!
	  G4warn << "ERROR: G4OpenInventorSceneHandler::PreAddSolid: Mother "
		 << ri->GetPhysicalVolume()->GetName()
		 << ':' << ri->GetCopyNo()
		 << " not previously encountered."
	    "\nShouldn't happen!  Please report to visualization coordinator."
		 << G4endl;
	  // Continue anyway.  Add to root of scene graph tree...
	  fCurrentSeparator = fDetectorRoot;
	}
      } else {
	fCurrentSeparator = fDetectorRoot;
      }
    }

  } else {
    // Not from G4PhysicalVolumeModel, so add to root as leaf part...

    if (fReadyForTransients) {
      fCurrentSeparator = fTransientRoot;
    } else {
      fCurrentSeparator = fDetectorRoot;
    }
  }
}

void G4OpenInventorSceneHandler::AddProperties(const G4VisAttributes* visAtts)
{
  // Use the applicable vis attributes...
  const G4VisAttributes* pApplicableVisAttribs =
    fpViewer->GetApplicableVisAttributes (visAtts);

  // First find the color attributes...
  const G4Colour& g4Col =  pApplicableVisAttribs->GetColour ();
  const G4double red = g4Col.GetRed ();
  const G4double green = g4Col.GetGreen ();
  const G4double blue = g4Col.GetBlue ();
  G4double transparency = 1 - g4Col.GetAlpha();

  // Drawing style...
  G4ViewParameters::DrawingStyle drawing_style =
    GetDrawingStyle(pApplicableVisAttribs);
  switch (drawing_style) {
  case (G4ViewParameters::wireframe):    
    fModelingSolid = false;
    break;
  case (G4ViewParameters::hlr):
  case (G4ViewParameters::hsr):
  case (G4ViewParameters::hlhsr):
    fModelingSolid = true;
    break;
  case (G4ViewParameters::cloud):
    fModelingSolid = false;
    break;
  }

  // Edge visibility...
  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pApplicableVisAttribs);
  fReducedWireFrame = !isAuxEdgeVisible;

  SoMaterial* material = 
    fStyleCache->getMaterial((G4float)red,
			     (G4float)green,
			     (G4float)blue,
			     (G4float)transparency);
  fCurrentSeparator->addChild(material);

  SoLightModel* lightModel = 
    fModelingSolid ? fStyleCache->getLightModelPhong() : 
    fStyleCache->getLightModelBaseColor();
  fCurrentSeparator->addChild(lightModel);
}

void G4OpenInventorSceneHandler::AddTransform(const G4Point3D& translation)
{
  // AddTransform takes fObjectTransformation and "adds" a translation.
  // Set up the geometrical transformation for the coming
  fCurrentSeparator->addChild(fStyleCache->getResetTransform());

  SoMatrixTransform* matrixTransform = new SoMatrixTransform;
  G4OpenInventorTransform3D oiTran
  (fObjectTransformation * G4Translate3D(translation));
  SbMatrix* sbMatrix = oiTran.GetSbMatrix();

  const G4Vector3D scale = fpViewer->GetViewParameters().GetScaleFactor();
  SbMatrix sbScale;
  sbScale.setScale
    (SbVec3f((G4float)scale.x(),(G4float)scale.y(),(G4float)scale.z()));
  sbMatrix->multRight(sbScale);

  matrixTransform->matrix.setValue(*sbMatrix);
  delete sbMatrix;
  fCurrentSeparator->addChild(matrixTransform);
}
