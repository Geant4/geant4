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
// Created by Stewart Boogert on 13/11/2021.
//

#include "G4VtkStore.hh"

#include "G4Circle.hh"
#include "G4Mesh.hh"
#include "G4Plane3D.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4Square.hh"
#include "G4Text.hh"
#include "G4VSceneHandler.hh"
#include "G4VisAttributes.hh"
#include "G4VtkUtility.hh"
#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"
// #include "G4VtkClipperClosedPipeline.hh"
#include "G4VtkClipClosedSurfacePipeline.hh"
#include "G4VtkClipOpenPipeline.hh"
#include "G4VtkCutterPipeline.hh"
#include "G4VtkImagePipeline.hh"
#include "G4VtkPolydataInstanceAppendPipeline.hh"
#include "G4VtkPolydataInstanceBakePipeline.hh"
#include "G4VtkPolydataInstanceTensorPipeline.hh"
#include "G4VtkPolydataPipeline.hh"
#include "G4VtkPolydataPolygonPipeline.hh"
#include "G4VtkPolydataPolyline2DPipeline.hh"
#include "G4VtkPolydataPolylinePipeline.hh"
#include "G4VtkPolydataSpherePipeline.hh"
#include "G4VtkText2DPipeline.hh"
#include "G4VtkTextPipeline.hh"
#include "G4VtkUnstructuredGridPipeline.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkBillboardTextActor3D.h"
#include "vtkLine.h"
#include "vtkMatrix3x3.h"
#include "vtkOBJReader.h"
#include "vtkPLYReader.h"
#include "vtkPolyDataReader.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkRegularPolygonSource.h"
#include "vtkSTLReader.h"
#include "vtkSmartPointer.h"
#include "vtkTextProperty.h"
#pragma GCC diagnostic pop

G4VtkStore::G4VtkStore(G4String nameIn) : name(nameIn) {}

G4VtkStore::~G4VtkStore()
{
  Clear();
}

void G4VtkStore::Modified()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkStore::Modified() " << name << G4endl;
#endif

  for (const auto& it : polylinePipeMap)
    it.second->Modified();

  for (const auto& it : polyline2DPipeMap)
    it.second->Modified();

  for (const auto& it : circlePipeMap)
    it.second->Modified();

  for (const auto& it : squarePipeMap)
    it.second->Modified();

  for (const auto& it : textPipeMap)
    it.second->Modified();

  for (const auto& it : text2DPipeMap)
    it.second->Modified();

  for (const auto& it : separatePipeMap)
    it.second->Modified();

  for (const auto& it : tensorGlyphPipeMap)
    it.second->Modified();

  for (const auto& it : appendPipeMap)
    it.second->Modified();

  for (const auto& it : bakePipeMap)
    it.second->Modified();
}

void G4VtkStore::Clear()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkStore::Clear() " << name << G4endl;
#endif

  for (const auto& v : polylinePipeMap) {
    v.second->Clear();
  }
  polylinePipeMap.clear();

  for (const auto& v : polyline2DPipeMap) {
    v.second->Clear();
  }
  polyline2DPipeMap.clear();

  for (const auto& v : circlePipeMap) {
    v.second->Clear();
  }
  circlePipeMap.clear();

  for (const auto& v : squarePipeMap) {
    v.second->Clear();
  }
  squarePipeMap.clear();

  for (const auto& v : textPipeMap) {
    v.second->Clear();
  }
  textPipeMap.clear();

  for (const auto& v : text2DPipeMap) {
    v.second->Clear();
  }
  text2DPipeMap.clear();

  for (const auto& v : separatePipeMap) {
    v.second->Clear();
  }
  separatePipeMap.clear();

  for (const auto& v : tensorGlyphPipeMap) {
    v.second->Clear();
  }
  tensorGlyphPipeMap.clear();

  for (const auto& v : appendPipeMap) {
    v.second->Clear();
  }
  appendPipeMap.clear();

  for (const auto& v : bakePipeMap) {
    v.second->Clear();
  }
  bakePipeMap.clear();
}

void G4VtkStore::ClearNonG4()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkStore::Clear() " << name << G4endl;
#endif

  for (const auto& v : imagePipeMap) {
    v.second->Clear();
  }
  imagePipeMap.clear();

  for (const auto& v : sideloadPipeMap) {
    v.second->Clear();
  }
  sideloadPipeMap.clear();
}

void G4VtkStore::Print()
{
  G4cout << "G4VtkStore::Print() " << name << std::endl;
  G4cout << "G4VtkStore::Print() polylinePipeMap size        > " << polylinePipeMap.size()
         << G4endl;
  G4cout << "G4VtkStore::Print() polyline2DPipeMap size      > " << polyline2DPipeMap.size()
         << G4endl;
  G4cout << "G4VtkStore::Print() circlePipeMap size          > " << circlePipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() squarePipeMap size          > " << squarePipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() textPipeMap size            > " << textPipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() text2DPipeMap size          > " << text2DPipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() separatePipeMap size        > " << separatePipeMap.size()
         << G4endl;
  G4cout << "G4VtkStore::Print() tensorGlyphPipeMap size     > " << tensorGlyphPipeMap.size()
         << G4endl;
  G4cout << "G4VtkStore::Print() transformAppendPipeMap size > " << appendPipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() bakePipeMap size            > " << bakePipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() ugridPipeMap size           > " << ugridPipeMap.size() << G4endl;
  G4cout << "G4VtkStore::Print() imagePipeMap size           > " << imagePipeMap.size() << G4endl;

  for (const auto& it : polylinePipeMap)
    it.second->Print();

  for (const auto& it : polyline2DPipeMap)
    it.second->Print();

  for (const auto& it : circlePipeMap)
    it.second->Print();

  for (const auto& it : squarePipeMap)
    it.second->Print();

  for (const auto& it : separatePipeMap)
    it.second->Print();

  for (const auto& it : tensorGlyphPipeMap)
    it.second->Print();

  for (const auto& it : appendPipeMap)
    it.second->Print();

  for (const auto& it : bakePipeMap)
    it.second->Print();

  for (const auto& it : ugridPipeMap) {
    G4cout << it.first << G4endl;
    it.second->Print();
  }

  for (const auto& it : imagePipeMap)
    it.second->Print();
}

void G4VtkStore::AddPrimitive(const G4Polyline& polyline, const G4VtkVisContext& vc)
{
  G4VSceneHandler::MarkerSizeType sizeType;
  if (vc.fProcessing2D) {
    sizeType = G4VSceneHandler::MarkerSizeType::screen;
  }
  else {
    sizeType = G4VSceneHandler::MarkerSizeType::world;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = vc.fViewer->GetApplicableVisAttributes(polyline.GetVisAttributes());

  if (sizeType == G4VSceneHandler::MarkerSizeType::world) {
    std::size_t hash = G4VtkPolydataPolylinePipeline::MakeHash(pVA);
    std::shared_ptr<G4VtkPolydataPolylinePipeline> pl;

    if (polylinePipeMap.find(hash) == polylinePipeMap.end()) {
      // Create new pipeline
      pl = std::make_shared<G4VtkPolydataPolylinePipeline>(G4String("name"), vc, pVA);

      // Add to map
      polylinePipeMap.insert(std::make_pair(hash, pl));
    }
    else {
      pl = polylinePipeMap[hash];
    }
    pl->SetPolydata(polyline);
  }

  else if (sizeType == G4VSceneHandler::MarkerSizeType::screen) {
    std::size_t hash = G4VtkPolydataPolyline2DPipeline::MakeHash(pVA);
    std::shared_ptr<G4VtkPolydataPolyline2DPipeline> pl;

    if (polyline2DPipeMap.find(hash) == polyline2DPipeMap.end()) {
      // Create new pipeline
      pl = std::make_shared<G4VtkPolydataPolyline2DPipeline>(G4String("name"), vc, pVA);

      // Add to map
      polyline2DPipeMap.insert(std::make_pair(hash, pl));
    }
    else {
      pl = polyline2DPipeMap[hash];
    }
    pl->SetPolydata(polyline);
  }
}

void G4VtkStore::AddPrimitive(const G4Text& text, const G4VtkVisContext& vc)
{
  G4VSceneHandler::MarkerSizeType sizeType;
  if (vc.fProcessing2D) {
    sizeType = G4VSceneHandler::MarkerSizeType::screen;
  }
  else {
    sizeType = G4VSceneHandler::MarkerSizeType::world;
  }

  auto pVA = vc.fViewer->GetApplicableVisAttributes(text.GetVisAttributes());

  switch (sizeType) {
    default:
    case (G4VSceneHandler::MarkerSizeType::screen): {
      std::size_t hash = G4VtkTextPipeline::MakeHash(text, vc, pVA);
      auto pl = std::make_shared<G4VtkText2DPipeline>(text, vc, pVA);
      text2DPipeMap.insert(std::make_pair(hash, pl));
      break;
    }
    case G4VSceneHandler::MarkerSizeType::world: {
      std::size_t hash = G4VtkTextPipeline::MakeHash(text, vc, pVA);
      auto pl = std::make_shared<G4VtkTextPipeline>(text, vc, pVA);
      textPipeMap.insert(std::make_pair(hash, pl));
      break;
    }
  }
}

void G4VtkStore::AddPrimitive(const G4Circle& circle, const G4VtkVisContext& vc)
{
  G4VSceneHandler::MarkerSizeType sizeType;

  if (vc.fProcessing2D) {
    sizeType = G4VSceneHandler::MarkerSizeType::screen;
  }
  else {
    sizeType = G4VSceneHandler::MarkerSizeType::world;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = vc.fViewer->GetApplicableVisAttributes(circle.GetVisAttributes());

  if (sizeType == G4VSceneHandler::MarkerSizeType::world) {
    std::size_t hash = G4VtkPolydataSpherePipeline::MakeHash(pVA);
    std::shared_ptr<G4VtkPolydataSpherePipeline> pl;

    if (circlePipeMap.find(hash) == circlePipeMap.end()) {
      // Create new pipeline
      pl = std::make_shared<G4VtkPolydataSpherePipeline>(G4String("name"), vc, pVA);

      // add to map
      circlePipeMap.insert(std::make_pair(hash, pl));
    }
    else {
      pl = circlePipeMap[hash];
    }

    // Data data point
    const CLHEP::HepRotation rot = vc.fTransform.getRotation();
    G4Point3D posPrime = rot * circle.GetPosition();
    pl->SetPolydataData(vc.fTransform.dx() + posPrime.x(), vc.fTransform.dy() + posPrime.y(),
                        vc.fTransform.dz() + posPrime.z());
  }
}

void G4VtkStore::AddPrimitive(const G4Square& square, const G4VtkVisContext& vc)
{
  G4VSceneHandler::MarkerSizeType sizeType;

  if (vc.fProcessing2D) {
    sizeType = G4VSceneHandler::MarkerSizeType::screen;
  }
  else {
    sizeType = G4VSceneHandler::MarkerSizeType::world;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = vc.fViewer->GetApplicableVisAttributes(square.GetVisAttributes());

  if (sizeType == G4VSceneHandler::MarkerSizeType::world) {
    std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);

    std::shared_ptr<G4VtkPolydataPolygonPipeline> pl;
    if (squarePipeMap.find(hash) == squarePipeMap.end()) {
      // Create new pipeline
      pl = std::make_shared<G4VtkPolydataPolygonPipeline>(G4String("name"), vc, pVA);

      // add to map
      squarePipeMap.insert(std::make_pair(hash, pl));
    }
    else {
      pl = squarePipeMap[hash];
    }

    // Data data point
    const CLHEP::HepRotation rot = vc.fTransform.getRotation();
    G4Point3D posPrime = rot * square.GetPosition();
    pl->SetPolydataData(vc.fTransform.dx() + posPrime.x(), vc.fTransform.dy() + posPrime.y(),
                        vc.fTransform.dz() + posPrime.z());
  }
}

void G4VtkStore::AddPrimitiveSeparate(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc)
{
  // Return if empty polyhedron
  if (polyhedron.GetNoFacets() == 0) {
    return;
  }

  auto hash = G4VtkPolydataPipeline::MakeHash(polyhedron, vc);

  auto pl = std::make_shared<G4VtkPolydataPipeline>(G4String("name"), vc);
  pl->SetPolydata(polyhedron);

  separatePipeMap.insert(std::make_pair(hash, pl));

  const G4VisAttributes* pVA =
    vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  G4Color colour = pVA->GetColour();
  G4Transform3D fInvObjTrans = vc.fTransform.inverse();

  pl->SetActorTransform(vc.fTransform.dx(), vc.fTransform.dy(), vc.fTransform.dz(),
                        fInvObjTrans.xx(), fInvObjTrans.xy(), fInvObjTrans.xz(), fInvObjTrans.yx(),
                        fInvObjTrans.yy(), fInvObjTrans.yz(), fInvObjTrans.zx(), fInvObjTrans.zy(),
                        fInvObjTrans.zz());
  pl->SetActorColour(colour.GetRed(), colour.GetGreen(), colour.GetBlue(), colour.GetAlpha());
}

void G4VtkStore::AddPrimitiveTensorGlyph(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc)
{
  // Return if empty polyhedron
  if (polyhedron.GetNoFacets() == 0) {
    return;
  }

  auto hash = G4VtkPolydataInstanceTensorPipeline::MakeHash(polyhedron, vc);
  std::shared_ptr<G4VtkPolydataInstanceTensorPipeline> pl;

  if (tensorGlyphPipeMap.find(hash) == tensorGlyphPipeMap.end()) {
    pl = std::make_shared<G4VtkPolydataInstanceTensorPipeline>(G4String("name"), vc);
    pl->SetPolydata(polyhedron);
    tensorGlyphPipeMap.insert(std::make_pair(hash, pl));
  }
  else {
    pl = tensorGlyphPipeMap[hash];
  }

  G4Transform3D fInvObjTrans = vc.fTransform.inverse();
  const G4VisAttributes* pVA =
    vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  G4Color colour = pVA->GetColour();

  pl->addInstance(vc.fTransform.dx(), vc.fTransform.dy(), vc.fTransform.dz(), fInvObjTrans.xx(),
                  fInvObjTrans.xy(), fInvObjTrans.xz(), fInvObjTrans.yx(), fInvObjTrans.yy(),
                  fInvObjTrans.yz(), fInvObjTrans.zx(), fInvObjTrans.zy(), fInvObjTrans.zz(),
                  colour.GetRed(), colour.GetGreen(), colour.GetBlue(), colour.GetAlpha(),
                  G4String("null"));
}

void G4VtkStore::AddPrimitiveAppend(const G4Polyhedron& polyhedron, const G4VtkVisContext& vc)
{
  // Empty polyhedron
  if (polyhedron.GetNoFacets() == 0) {
    return;
  }

  auto hash = G4VtkPolydataInstanceAppendPipeline::MakeHash(polyhedron, vc);
  std::shared_ptr<G4VtkPolydataInstanceAppendPipeline> pl;

  if (appendPipeMap.find(hash) == appendPipeMap.end()) {
    pl = std::make_shared<G4VtkPolydataInstanceAppendPipeline>(G4String("name"), vc);
    pl->SetPolydata(polyhedron);

    appendPipeMap.insert(std::make_pair(hash, pl));
  }
  else {
    pl = appendPipeMap[hash];
  }

  G4Transform3D fInvObjTrans = vc.fTransform.inverse();

  pl->addInstance(vc.fTransform.dx(), vc.fTransform.dy(), vc.fTransform.dz(), fInvObjTrans.xx(),
                  fInvObjTrans.xy(), fInvObjTrans.xz(), fInvObjTrans.yx(), fInvObjTrans.yy(),
                  fInvObjTrans.yz(), fInvObjTrans.zx(), fInvObjTrans.zy(), fInvObjTrans.zz(), 0, 0,
                  0, 0, G4String("null"));
};

void G4VtkStore::AddPrimitiveTransformBake(const G4Polyhedron& polyhedron,
                                           const G4VtkVisContext& vc)
{
  // Empty polyhedron
  if (polyhedron.GetNoFacets() == 0) {
    return;
  }

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());

  // Polyhedron colour
  G4Colour colour = pVA->GetColour();

  // Hash the vis (alpha) attributes and polyhedron
  std::size_t hash = G4VtkPolydataInstanceBakePipeline::MakeHash(polyhedron, vc);

  std::shared_ptr<G4VtkPolydataInstanceBakePipeline> pl;
  if (bakePipeMap.find(hash) == bakePipeMap.end()) {
    pl = std::make_shared<G4VtkPolydataInstanceBakePipeline>(G4String("none"), vc);

    bakePipeMap.insert(std::make_pair(hash, pl));
  }
  else {
    pl = bakePipeMap[hash];
  }

  G4Transform3D fInvObjTrans = vc.fTransform.inverse();
  pl->SetPolydata(polyhedron);

  pl->addInstance(vc.fTransform.dx(), vc.fTransform.dy(), vc.fTransform.dz(), fInvObjTrans.xx(),
                  fInvObjTrans.xy(), fInvObjTrans.xz(), fInvObjTrans.yx(), fInvObjTrans.yy(),
                  fInvObjTrans.yz(), fInvObjTrans.zx(), fInvObjTrans.zy(), fInvObjTrans.zz(),
                  colour.GetRed(), colour.GetGreen(), colour.GetBlue(), colour.GetAlpha(),
                  G4String("null"));
};

void G4VtkStore::AddCompound(const G4Mesh& mesh, const G4VtkVisContext& vc) {
  const auto& container = mesh.GetContainerVolume();
  auto& containerName = const_cast<G4String&>(container->GetName());

  auto pl = std::make_shared<G4VtkUnstructuredGridPipeline>(containerName,vc);
  pl->SetUnstructuredGridData(mesh);
  ugridPipeMap.insert(std::make_pair(G4String("test"), pl));
}

void G4VtkStore::UpdatePlanePipelines(G4String nameIn, G4String type, const G4Plane3D plane)
{
  if (type == "clipper") {
    UpdateClipper(nameIn, plane);
  }
  else if (type == "cutter") {
    UpdateCutter(nameIn, plane);
  }
}

void G4VtkStore::AddClipper(G4String nameIn, const G4Plane3D& plane)
{
  for (const auto& v : separatePipeMap) {
    auto clip_pl = new G4VtkClipClosedSurfacePipeline(nameIn, v.second->GetVtkVisContext(),
                                                      v.second->GetFinalFilter(), true);
    clip_pl->SetPlane(plane);
    v.second->AddChildPipeline(clip_pl);
  }

  for (const auto& v : tensorGlyphPipeMap) {
    auto clip_pl = new G4VtkClipClosedSurfacePipeline(nameIn, v.second->GetVtkVisContext(),
                                                      v.second->GetFinalFilter());
    clip_pl->SetPlane(plane);
    v.second->AddChildPipeline(clip_pl);
  }

  for (const auto& v : appendPipeMap) {
    auto clip_pl = new G4VtkClipClosedSurfacePipeline(nameIn, v.second->GetVtkVisContext(),
                                                      v.second->GetFinalFilter(), true);
    clip_pl->SetPlane(plane);
    v.second->AddChildPipeline(clip_pl);
  }

  for (const auto& v : bakePipeMap) {
    auto clip_pl = new G4VtkClipClosedSurfacePipeline(nameIn, v.second->GetVtkVisContext(),
                                                      v.second->GetFinalFilter(), true);
    clip_pl->SetPlane(plane);
    v.second->AddChildPipeline(clip_pl);
  }
}

void G4VtkStore::UpdateClipper(G4String nameIn, const G4Plane3D& plane)
{
  for (const auto& v : separatePipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkClipClosedSurfacePipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : tensorGlyphPipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkClipClosedSurfacePipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : appendPipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkClipClosedSurfacePipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : bakePipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkClipClosedSurfacePipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }
}

void G4VtkStore::RemoveClipper(G4String /*nameIn*/) {}

void G4VtkStore::AddCutter(G4String nameIn, const G4Plane3D& plane)
{
  for (const auto& v : separatePipeMap) {
    auto cut_pl = new G4VtkCutterPipeline(nameIn, v.second->GetVtkVisContext(),
                                          v.second->GetFinalFilter(), true);
    cut_pl->SetPlane(plane);
    v.second->AddChildPipeline(cut_pl);
  }

  for (const auto& v : tensorGlyphPipeMap) {
    auto cut_pl =
      new G4VtkCutterPipeline(nameIn, v.second->GetVtkVisContext(), v.second->GetFinalFilter());
    cut_pl->SetPlane(plane);
    v.second->AddChildPipeline(cut_pl);
  }

  for (const auto& v : appendPipeMap) {
    auto cut_pl =
      new G4VtkCutterPipeline(nameIn, v.second->GetVtkVisContext(), v.second->GetFinalFilter());
    cut_pl->SetPlane(plane);
    v.second->AddChildPipeline(cut_pl);
  }

  for (const auto& v : bakePipeMap) {
    auto cut_pl =
      new G4VtkCutterPipeline(nameIn, v.second->GetVtkVisContext(), v.second->GetFinalFilter());
    cut_pl->SetPlane(plane);
    v.second->AddChildPipeline(cut_pl);
  }
}

void G4VtkStore::UpdateCutter(G4String nameIn, const G4Plane3D& plane)
{
  for (const auto& v : separatePipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkCutterPipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : tensorGlyphPipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkCutterPipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : appendPipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkCutterPipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }

  for (const auto& v : bakePipeMap) {
    auto children = v.second->GetChildPipelines();
    for (auto c : children) {
      if (c->GetName() == nameIn) {
        auto cc = dynamic_cast<G4VtkCutterPipeline*>(c);
        cc->SetPlane(plane);
      }
    }
  }
}

void G4VtkStore::RemoveCutter(G4String /*nameIn*/) {}

void G4VtkStore::AddNonG4ObjectImage(const G4String& fileName, const G4VtkVisContext& vc)
{
  auto pl = std::make_shared<G4VtkImagePipeline>(G4String("name"), vc);
  pl->SetImage(fileName);
  imagePipeMap.insert(std::make_pair(fileName, pl));
}

void G4VtkStore::AddNonG4ObjectPolydata(const G4String fileName, const G4VtkVisContext& /*vc*/)
{
  if (fileName.find("obj") != G4String::npos) {
    vtkSmartPointer<vtkOBJReader> objReader;
  }
  else if (fileName.find("ply") != G4String::npos) {
    vtkSmartPointer<vtkPLYReader> plyReader;
  }
  else if (fileName.find("stl") != G4String::npos) {
    vtkSmartPointer<vtkSTLReader> stlReader;
  }
  else if (fileName.find("vtp") != G4String::npos && fileName.find("vtu") != G4String::npos) {
    vtkSmartPointer<vtkPolyDataReader> pldReader;
  }
}

void G4VtkStore::GetBounds(G4double maxBoundIn[6])
{
  G4double maxBound[6] = {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};

  for (const auto& v : separatePipeMap) {
    auto b = v.second->GetBounds();
    MaxBounds(maxBound, b);
  }

  for (const auto& v : tensorGlyphPipeMap) {
    auto b = v.second->GetBounds();
    MaxBounds(maxBound, b);
  }

  for (const auto& v : appendPipeMap) {
    auto b = v.second->GetBounds();
    MaxBounds(maxBound, b);
  }

  for (const auto& v : bakePipeMap) {
    auto b = v.second->GetBounds();
    MaxBounds(maxBound, b);
  }

  maxBoundIn[0] = maxBound[0];
  maxBoundIn[1] = maxBound[1];
  maxBoundIn[2] = maxBound[2];
  maxBoundIn[3] = maxBound[3];
  maxBoundIn[4] = maxBound[4];
  maxBoundIn[5] = maxBound[5];
}

void G4VtkStore::AddToRenderer(vtkRenderer* renderer)
{
  for (const auto& it : polylinePipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : polyline2DPipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : circlePipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : squarePipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : textPipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : text2DPipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : separatePipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : tensorGlyphPipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : appendPipeMap)
    renderer->AddActor(it.second->GetActor());
  for (const auto& it : bakePipeMap)
    renderer->AddActor(it.second->GetActor());
}
