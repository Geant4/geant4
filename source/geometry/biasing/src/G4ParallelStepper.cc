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
// $Id: G4ParallelStepper.cc,v 1.12 2006/06/29 18:17:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelStepper.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelStepper.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelStepper::G4ParallelStepper() :
  G4VParallelStepper(),
  fPStep(0)
{}

G4ParallelStepper::~G4ParallelStepper()
{
  if (fPStep) {
    delete fPStep;
  }
}

G4ParallelStepper::G4ParallelStepper(const G4ParallelStepper &rhs)
  : G4VParallelStepper(), fPStep(new G4GeometryCellStep(rhs.GetPStep()))
{
  if (!fPStep) {
    Error("G4ParallelStepper:: new failed to create a G4GeometryCellStep!");
  }
}

G4ParallelStepper &G4ParallelStepper::operator=(const G4ParallelStepper &rhs)
{
  if (this != &rhs) {
    fPStep = new G4GeometryCellStep(rhs.GetPStep());
    if (!fPStep) {
      Error("operator=: new failed to create a G4GeometryCellStep!");
    }
  }
  return *this;
}

G4GeometryCellStep G4ParallelStepper::GetPStep() const {
  G4GeometryCellStep p = *fPStep;
  return p;
}


void G4ParallelStepper::Init(const G4GeometryCell &agCell)
{
  if (!fPStep) {
    fPStep = new G4GeometryCellStep(agCell, agCell);
      if (!fPStep) {
	Error("Init new failed to create a G4GeometryCellStep!");
      }

  }
  else {
    fPStep->SetPreGeometryCell(agCell);
    fPStep->SetPostGeometryCell(agCell);
    fPStep->SetCrossBoundary(false);
  }
}

void G4ParallelStepper::Update(const G4GeometryCell &agCell)
{
  if (!fPStep) {
    Error("fPStep == 0, Init not called?");
  }
  fPStep->SetPreGeometryCell(fPStep->GetPostGeometryCell());
  fPStep->SetPostGeometryCell(agCell);
  fPStep->SetCrossBoundary(true);
}

void G4ParallelStepper::UnSetCrossBoundary()
{
  if (!fPStep) {
    Error("fPStep == 0, Init not called?");
  }
  fPStep->SetCrossBoundary(false);
}

void G4ParallelStepper::Error(const G4String &m)
{
  G4cout << "ERROR: in G4ParallelStepper: " << m << G4endl;
  G4Exception("G4ParallelStepper::Error()",
              "FatalException", FatalException, m);
}
