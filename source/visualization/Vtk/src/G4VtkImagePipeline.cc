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

#include "G4VtkImagePipeline.hh"

#include "G4String.hh"
#include "G4ViewParameters.hh"
#include "G4VisAttributes.hh"
#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkImageActor.h"
#include "vtkImageAlgorithm.h"
#include "vtkImageMapper3D.h"
#include "vtkImageProperty.h"
#include "vtkImageReader2.h"
#include "vtkImageReader2Factory.h"
#include "vtkMatrix4x4.h"
#include "vtkSmartPointer.h"
#include "vtkImageFlip.h"

G4VtkImagePipeline::G4VtkImagePipeline(const G4String& nameIn, const G4VtkVisContext& vcIn)
  : G4VVtkPipeline(nameIn, "G4VtkImagePipeline", vcIn, false, vcIn.fViewer->renderer)
{
  actor = vtkSmartPointer<vtkImageActor>::New();

  actor->GetProperty()->SetOpacity(vc.alpha);

  auto transform = vtkSmartPointer<vtkMatrix4x4>::New();
  double transformArray[16] = {vc.fTransform.xx(),
                               vc.fTransform.xy(),
                               vc.fTransform.xz(),
                               vc.fTransform.dx(),
                               vc.fTransform.yx(),
                               vc.fTransform.yy(),
                               vc.fTransform.yz(),
                               vc.fTransform.dy(),
                               vc.fTransform.zx(),
                               vc.fTransform.zy(),
                               vc.fTransform.zz(),
                               vc.fTransform.dz(),
                               0.,
                               0.,
                               0.,
                               1.};
  transform->DeepCopy(transformArray);
  actor->SetUserMatrix(transform);
  actor->GetProperty()->SetOpacity(vc.alpha);

  vc.fViewer->renderer->AddActor(GetActor());
}

void G4VtkImagePipeline::SetImage(const G4String& fileName)
{
  vtkNew<vtkImageReader2Factory> readerFactory;
  vtkSmartPointer<vtkImageReader2> imageReader;
  imageReader.TakeReference(readerFactory->CreateImageReader2(fileName.c_str()));
  imageReader->SetFileName(fileName.c_str());
  imageReader->Update();

  vtkNew<vtkImageFlip> flip;
  flip->SetInputConnection(imageReader->GetOutputPort());
  flip->SetFilteredAxis(1);

  actor->GetMapper()->SetInputConnection(flip->GetOutputPort());
  AddFilter(imageReader);
}

void G4VtkImagePipeline::SetTransformation(const G4Transform3D& transformation)
{
  vtkSmartPointer<vtkMatrix4x4> t = vtkSmartPointer<vtkMatrix4x4>::New();
  t->SetElement(0, 0, transformation.xx());
  t->SetElement(0, 1, transformation.xy());
  t->SetElement(0, 2, transformation.xz());
  t->SetElement(0, 3, transformation.dx());

  t->SetElement(1, 0, transformation.yx());
  t->SetElement(1, 1, transformation.yy());
  t->SetElement(1, 2, transformation.yz());
  t->SetElement(1, 3, transformation.dy());

  t->SetElement(2, 0, transformation.zx());
  t->SetElement(2, 1, transformation.zy());
  t->SetElement(2, 2, transformation.zz());
  t->SetElement(2, 3, transformation.dz());

  t->SetElement(3, 0, 0);
  t->SetElement(3, 1, 0);
  t->SetElement(3, 2, 0);
  t->SetElement(3, 3, 1);
  actor->SetUserMatrix(t);
}

void G4VtkImagePipeline::Print()
{
  for (auto c : childPipelines)
    c->Print();
}

void G4VtkImagePipeline::Modified()
{
  actor->Update();

  for (auto c : childPipelines)
    c->Modified();
}

void G4VtkImagePipeline::Clear()
{
  if (actor != nullptr) renderer->RemoveActor(actor);

  for (auto c : childPipelines)
    c->Clear();
}

void G4VtkImagePipeline::Disable()
{
  GetActor()->SetVisibility(0);
}

void G4VtkImagePipeline::Enable()
{
  GetActor()->SetVisibility(1);
}