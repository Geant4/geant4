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
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#include "G4VtkSceneHandler.hh"

#include "G4Box.hh"
#include "G4Circle.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4Material.hh"
#include "G4Mesh.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4PseudoScene.hh"
#include "G4Square.hh"
#include "G4SystemOfUnits.hh"
#include "G4Text.hh"
#include "G4UnitsTable.hh"
#include "G4VNestedParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VtkStore.hh"
#include "G4VtkVisContext.hh"

#include "vtkColorTransferFunction.h"
#include "vtkContourValues.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"

#include <cstdlib>

G4int G4VtkSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4VtkSceneHandler::G4VtkSceneHandler(G4VGraphicsSystem& system, const G4String& name)
  : G4VSceneHandler(system, fSceneIdCount++, name), polyhedronPipelineType(G4String("tensor"))
{}

void G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline&)" << G4endl;
#endif
  auto vc = MakeDefaultVisContext();

  if (fReadyForTransients)
    transientStore.AddPrimitive(polyline, vc);
  else
    store.AddPrimitive(polyline, vc);
}

void G4VtkSceneHandler::AddPrimitive(const G4Text& text)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Text& text)" << G4endl;
#endif

  auto vc = MakeDefaultVisContext();

  if (fReadyForTransients)
    transientStore.AddPrimitive(text, vc);
  else
    store.AddPrimitive(text, vc);
}

void G4VtkSceneHandler::AddPrimitive(const G4Circle& circle)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Circle& circle)" << G4endl;
#endif

  auto vc = MakeDefaultVisContext();
  G4VSceneHandler::MarkerSizeType sizeType;
  vc.fSize = GetMarkerSize(circle, sizeType);

  if (fReadyForTransients)
    transientStore.AddPrimitive(circle, vc);
  else
    store.AddPrimitive(circle, vc);
}

void G4VtkSceneHandler::AddPrimitive(const G4Square& square)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Square& square)" << G4endl;
#endif

  auto vc = MakeDefaultVisContext();
  G4VSceneHandler::MarkerSizeType sizeType;
  vc.fSize = GetMarkerSize(square, sizeType);

  if (fReadyForTransients)
    transientStore.AddPrimitive(square, vc);
  else
    store.AddPrimitive(square, vc);
}

void G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)" << G4endl;
#endif

  auto vc = MakeDefaultVisContext();
  auto visAtt = vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  auto colour = visAtt->GetColour();

  vc.fDrawingStyle = GetDrawingStyle(visAtt);
  vc.alpha = colour.GetAlpha();
  vc.red = colour.GetRed();
  vc.green = colour.GetGreen();
  vc.blue = colour.GetBlue();

  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel != nullptr) {
    vc.fDepth = pPVModel->GetCurrentDepth();
    vc.fDescription = pPVModel->GetCurrentDescription();
  }

  if (fReadyForTransients) {
    if (polyhedronPipelineType == "tensor")
      transientStore.AddPrimitiveTensorGlyph(polyhedron, vc);
    else if (polyhedronPipelineType == "append")
      transientStore.AddPrimitiveAppend(polyhedron, vc);
    else if (polyhedronPipelineType == "bake")
      transientStore.AddPrimitiveTransformBake(polyhedron, vc);
    else if (polyhedronPipelineType == "separate")
      transientStore.AddPrimitiveSeparate(polyhedron, vc);
  }
  else {
    if (polyhedronPipelineType == "tensor")
      store.AddPrimitiveTensorGlyph(polyhedron, vc);
    else if (polyhedronPipelineType == "append")
      store.AddPrimitiveAppend(polyhedron, vc);
    else if (polyhedronPipelineType == "bake")
      store.AddPrimitiveTransformBake(polyhedron, vc);
    else if (polyhedronPipelineType == "separate")
      store.AddPrimitiveSeparate(polyhedron, vc);
  }
}

void G4VtkSceneHandler::Modified()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::Modified()" << G4endl;
#endif

  store.Modified();
  transientStore.Modified();
}

void G4VtkSceneHandler::ClearStore()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::ClearStore()" << G4endl;
#endif
  store.Clear();
}

void G4VtkSceneHandler::ClearTransientStore()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::ClearTransientStore()" << G4endl;
#endif
  transientStore.Clear();
}

G4VtkVisContext G4VtkSceneHandler::MakeDefaultVisContext()
{
  auto vc = G4VtkVisContext(dynamic_cast<G4VtkViewer*>(fpViewer), fpVisAttribs, fProcessing2D,
                            fObjectTransformation);

  if (fpVisAttribs != nullptr) {
    G4Colour c = fpVisAttribs->GetColour();
    vc.red = c.GetRed();
    vc.green = c.GetGreen();
    vc.blue = c.GetBlue();
    vc.alpha = c.GetAlpha();
    vc.fDrawingStyle = fpViewer->GetViewParameters().GetDrawingStyle();
  }

  return vc;
}

void G4VtkSceneHandler::AddSolid(const G4Box& box)
{
  G4VSceneHandler::AddSolid(box);

  return;

  const G4VModel* pv_model = GetModel();
  if (pv_model == nullptr) {
    return;
  }

  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel == nullptr) {
    return;
  }

  //-- debug information
#ifdef G4VTKDEBUG
  G4VPhysicalVolume* pv = pPVModel->GetCurrentPV();
  G4LogicalVolume* lv = pv->GetLogicalVolume();
  G4cout << "name=" << box.GetName() << " volumeType=" << pv->VolumeType()
         << " pvName=" << pv->GetName() << " lvName=" << lv->GetName()
         << " multiplicity=" << pv->GetMultiplicity() << " isparametrised=" << pv->IsParameterised()
         << " isreplicated=" << pv->IsReplicated()
         << " parametrisation=" << pv->GetParameterisation()
         << G4endl

              G4Material* mat = pPVModel->GetCurrentMaterial();
  G4String name = mat->GetName();
  G4double dens = mat->GetDensity() / (g / cm3);
  G4int copyNo = pPVModel->GetCurrentPV()->GetCopyNo();
  G4int depth = pPVModel->GetCurrentDepth();
  G4cout << "    name    : " << box.GetName() << G4endl;
  G4cout << "    copy no.: " << copyNo << G4endl;
  G4cout << "    depth   : " << depth << G4endl;
  G4cout << "    density : " << dens << " [g/cm3]" << G4endl;
  G4cout << "    location: " << pPVModel->GetCurrentPV()->GetObjectTranslation() << G4endl;
  G4cout << "    Multiplicity        : " << pPVModel->GetCurrentPV()->GetMultiplicity() << G4endl;
  G4cout << "    Is replicated?      : " << pPVModel->GetCurrentPV()->IsReplicated() << G4endl;
  G4cout << "    Is parameterised?   : " << pPVModel->GetCurrentPV()->IsParameterised() << G4endl;
  G4cout << "    top phys. vol. name : " << pPVModel->GetTopPhysicalVolume()->GetName() << G4endl;
#endif
}

void G4VtkSceneHandler::AddCompound(const G4Mesh& mesh)
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddCompound> mesh type " << mesh.GetMeshType() << " "
         << fpViewer->GetViewParameters().GetSpecialMeshRenderingOption() << G4endl;
#endif

  if(fpViewer->GetViewParameters().GetSpecialMeshRenderingOption() == G4ViewParameters::meshAsDefault)
  {
    auto vc = MakeDefaultVisContext();

    if (fReadyForTransients)
      transientStore.AddCompound(mesh, vc);
    else
      store.AddCompound(mesh, vc);
  }
  else {
    StandardSpecialMeshRendering(mesh);
  }
}

void G4VtkSceneHandler::Print() {}

void G4VtkSceneHandler::SetPolyhedronPipeline(const G4String& type)
{
  polyhedronPipelineType = type;
}
