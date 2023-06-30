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

#include "G4VtkPolydataInstancePipeline.hh"

#include "vtkDoubleArray.h"
#include "vtkPoints.h"

G4VtkPolydataInstancePipeline::G4VtkPolydataInstancePipeline(G4String nameIn,
                                                             const G4VtkVisContext& vc)
  : G4VtkPolydataPipeline(nameIn, vc)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataInstancePipeline"));

  instanceColour = vtkSmartPointer<vtkDoubleArray>::New();
  instanceColour->SetName("colours");
  instanceColour->SetNumberOfComponents(4);

  instancePosition = vtkSmartPointer<vtkPoints>::New();

  instanceTransform = vtkSmartPointer<vtkDoubleArray>::New();
  instanceTransform->SetName("transform");
  instanceTransform->SetNumberOfComponents(9);
}

void G4VtkPolydataInstancePipeline::addInstance(G4double dx, G4double dy, G4double dz, G4double r00,
                                                G4double r01, G4double r02, G4double r10,
                                                G4double r11, G4double r12, G4double r20,
                                                G4double r21, G4double r22, const G4String& name)
{
  // add the instance without colour or alpha
  addInstance(dx, dy, dz, r00, r01, r02, r10, r11, r12, r20, r21, r22, 0, 0, 0, 0, name);
}

void G4VtkPolydataInstancePipeline::addInstance(G4double dx, G4double dy, G4double dz, G4double r00,
                                                G4double r01, G4double r02, G4double r10,
                                                G4double r11, G4double r12, G4double r20,
                                                G4double r21, G4double r22, G4double r, G4double g,
                                                G4double b, G4double a, const G4String& name)
{
  instanceColour->InsertNextTuple4(r, g, b, a);
  auto idp = instancePosition->InsertNextPoint(dx, dy, dz);
  instanceTransform->InsertNextTuple9(r00, r01, r02, r10, r11, r12, r20, r21, r22);
  instanceMap[name] = idp;
}

void G4VtkPolydataInstancePipeline::removeInstance(const G4String& /*name*/)
{
  // remove from instanceColour, instancePosition and instanceTransform
}