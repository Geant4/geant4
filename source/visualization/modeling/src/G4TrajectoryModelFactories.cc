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
// $Id: G4TrajectoryModelFactories.cc,v 1.1 2005-11-21 05:44:44 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005

#include "G4TrajectoryModelFactories.hh"
#include "G4ModelCommandsDrawByCharge.hh"
#include "G4TrajectoryDrawByCharge.hh"

G4TrajectoryDrawByChargeFactory::G4TrajectoryDrawByChargeFactory()
  :G4VModelFactory<G4VTrajectoryModel>("drawByCharge") 
{}

G4TrajectoryDrawByChargeFactory::~G4TrajectoryDrawByChargeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByChargeFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryDrawByCharge* model = new G4TrajectoryDrawByCharge(name);
  
  // Create associated messengers
  Messengers messengers;

  messengers.push_back(new G4ModelCommandDrawByChargeSet(model, placement));
  messengers.push_back(new G4ModelCommandDrawByChargeSetRGBA(model, placement));
  
  return ModelAndMessengers(model, messengers);
}
