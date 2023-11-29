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

#include "G4VtkPolydataInstanceAppendPipeline.hh"

#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkGeneralTransform.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkTransformPolyDataFilter.h"

std::size_t G4VtkPolydataInstanceAppendPipeline::MakeHash(const G4Polyhedron& polyhedron,
                                                          const G4VtkVisContext& vc)
{
  // Get view parameters that the user can force through the vis attributes, thereby over-riding the
  // current view parameter.
  const G4VisAttributes* pVA =
    vc.fViewer->GetApplicableVisAttributes(polyhedron.GetVisAttributes());
  G4Color colour = pVA->GetColour();

  // Hash the vis attributes
  std::size_t hash = std::hash<G4double>{}(colour.GetAlpha());
  std::size_t rhash = std::hash<G4double>{}(colour.GetRed());
  std::size_t ghash = std::hash<G4double>{}(colour.GetGreen());
  std::size_t bhash = std::hash<G4double>{}(colour.GetBlue());
  std::size_t phash = std::hash<G4Polyhedron>{}(polyhedron);
  std::size_t shash = std::hash<G4double>{}(vc.fDrawingStyle);

  std::hash_combine(hash, phash);
  std::hash_combine(hash, rhash);
  std::hash_combine(hash, bhash);
  std::hash_combine(hash, ghash);
  std::hash_combine(hash, shash);

  return hash;
}

G4VtkPolydataInstanceAppendPipeline::G4VtkPolydataInstanceAppendPipeline(G4String nameIn,
                                                                         const G4VtkVisContext& vc)
  : G4VtkPolydataInstancePipeline(nameIn, vc)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataInstanceAppendPipeline"));

  // append filter
  appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
  AddFilter(appendFilter);

  // set polydata mapper
  mapper->SetInputConnection(GetFinalFilter()->GetOutputPort());

  // set actor
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // colour parameters
  actor->GetProperty()->SetOpacity(vc.alpha);
  actor->GetProperty()->SetColor(vc.red, vc.green, vc.blue);

  // shading parameters
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.1);
  actor->GetProperty()->SetSpecularPower(1);

  // add to renderer
  vc.fViewer->renderer->AddActor(GetActor());
}

void G4VtkPolydataInstanceAppendPipeline::Print()
{
  G4cout << "G4VtkPolydataInstanceBakePipeline " << GetName() << G4endl;
  G4VtkPolydataPipeline::Print();
}

void G4VtkPolydataInstanceAppendPipeline::addInstance(G4double dx, G4double dy, G4double dz,
                                                      G4double r00, G4double r01, G4double r02,
                                                      G4double r10, G4double r11, G4double r12,
                                                      G4double r20, G4double r21, G4double r22,
                                                      G4double r, G4double g, G4double b,
                                                      G4double a, const G4String& name)
{
  // Add to base class
  G4VtkPolydataInstancePipeline::addInstance(dx, dy, dz, r00, r01, r02, r10, r11, r12, r20, r21,
                                             r22, r, g, b, a, name);

  // create transform
  auto transform = vtkSmartPointer<vtkGeneralTransform>::New();
  double transformArray[16] = {r00, r01, r02, dx, r10, r11, r12, dy,
                               r20, r21, r22, dz, 0.,  0.,  0.,  1.};
  transform->Concatenate(transformArray);
  transform->Update();

  // Create transform filter and add to local filter map
  vtkSmartPointer<vtkTransformPolyDataFilter> tf =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  tf->SetTransform(transform);
  tf->SetInputConnection(GetFilter(GetNumberOfFilters() - 2)->GetOutputPort());
  tf->Update();
  transformFilterMap[name] = tf;

  // Add transform filter to append filter
  GetFinalFilter()->AddInputConnection(tf->GetOutputPort());
}

void G4VtkPolydataInstanceAppendPipeline::removeInstance(const G4String& name)
{
  // Remove from base class
  G4VtkPolydataInstancePipeline::removeInstance(name);

  // Remove transform filter from append filter
  appendFilter->RemoveInputConnection(0, transformFilterMap[name]->GetOutputPort());

  // Remove from local filter map
  transformFilterMap.erase(name);
}
