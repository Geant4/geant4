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
// $Id: G4VTrajectoryModel.cc,v 1.2 2006-05-05 01:46:33 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
