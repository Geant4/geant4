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
// $Id: G4ParallelStepper.cc,v 1.6 2002-10-14 12:36:03 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelStepper.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelStepper.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelStepper::G4ParallelStepper()
 : fPStep(0)
{}

G4ParallelStepper::~G4ParallelStepper()
{
  if (fPStep) {
    delete fPStep;
  }
}

G4ParallelStepper::G4ParallelStepper(const G4ParallelStepper &rhs)
  :
  fPStep(new G4PStep(rhs.GetPStep()))
{
  if (!fPStep) {
    Error("G4ParallelStepper:: new failed to create a G4PStep!");
  }
}

G4ParallelStepper &G4ParallelStepper::operator=(const G4ParallelStepper &rhs)
{
  if (this != &rhs) {
    fPStep = new G4PStep(rhs.GetPStep());
    if (!fPStep) {
      Error("operator=: new failed to create a G4PStep!");
    }
  }
  return *this;
}

G4PStep G4ParallelStepper::GetPStep() const {
  G4PStep p = *fPStep;
  return p;
}


void G4ParallelStepper::Init(const G4GeometryCell &agCell)
{
  if (!fPStep) {
    fPStep = new G4PStep(agCell, agCell);
      if (!fPStep) {
	Error("Init new failed to create a G4PStep!");
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
  G4std::G4cout << "ERROR: in G4ParallelStepper::" << m << G4endl;
  G4std::G4Exception("Program aborted.");
}
