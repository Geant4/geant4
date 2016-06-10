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
// $Id: G4VTrajectoryModel.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay May 2006
//
// See header for details
//
#include "G4VTrajectoryModel.hh"
#include "G4VisTrajContext.hh"

#include <assert.h>

G4VTrajectoryModel::G4VTrajectoryModel(const G4String& name, 
				       G4VisTrajContext* context)
  :fName(name)
  ,fVerbose(false)
  ,fpContext(context) 
{
  // Create context object if none is provided. Model will
  // then use default G4VisTrajContext configuration
  if (0 == fpContext) fpContext = new G4VisTrajContext();
}

G4VTrajectoryModel::~G4VTrajectoryModel()
{
  delete fpContext;
}

const G4VisTrajContext&
G4VTrajectoryModel::GetContext() const 
{
  // Expect context to exist
  assert (0 != fpContext);
  return *fpContext;
}

G4String 
G4VTrajectoryModel::Name() const 
{
  return fName;
}

void
G4VTrajectoryModel::SetVerbose(const G4bool& verbose)
{
  fVerbose = verbose;
}

G4bool
G4VTrajectoryModel::GetVerbose() const
{
  return fVerbose;
}
