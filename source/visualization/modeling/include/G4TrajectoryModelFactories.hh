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
// $Id: G4TrajectoryModelFactories.hh,v 1.2 2005/11/23 05:19:23 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Trajectory model factories creating trajectory models
// and associated messengers
// Class Description - End:

#ifndef G4TRAJECTORYMODELFACTORIES_HH
#define G4TRAJECTORYMODELFACTORIES_HH

#include "G4VModelFactory.hh"
#include "G4VTrajectoryModel.hh"

namespace {
  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair<G4VTrajectoryModel*, Messengers > ModelAndMessengers;
}

class G4TrajectoryDrawByChargeFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByChargeFactory();

  virtual ~G4TrajectoryDrawByChargeFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

class G4TrajectoryDrawByParticleIDFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByParticleIDFactory();

  virtual ~G4TrajectoryDrawByParticleIDFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

#endif


