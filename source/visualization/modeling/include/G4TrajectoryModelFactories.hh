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
// $Id: G4TrajectoryModelFactories.hh 98766 2016-08-09 14:17:17Z gcosmo $
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

class G4TrajectoryDrawByAttributeFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByAttributeFactory();

  virtual ~G4TrajectoryDrawByAttributeFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

class G4TrajectoryDrawByChargeFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByChargeFactory();

  virtual ~G4TrajectoryDrawByChargeFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

class G4TrajectoryGenericDrawerFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryGenericDrawerFactory();

  virtual ~G4TrajectoryGenericDrawerFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

class G4TrajectoryDrawByParticleIDFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByParticleIDFactory();

  virtual ~G4TrajectoryDrawByParticleIDFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

class G4TrajectoryDrawByOriginVolumeFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByOriginVolumeFactory();

  virtual ~G4TrajectoryDrawByOriginVolumeFactory();

  ModelAndMessengers Create(const G4String& placement, const G4String& name);

};

class G4TrajectoryDrawByEncounteredVolumeFactory : public G4VModelFactory<G4VTrajectoryModel> {

public: // With description

  G4TrajectoryDrawByEncounteredVolumeFactory();

  virtual ~G4TrajectoryDrawByEncounteredVolumeFactory();

  ModelAndMessengers Create(const G4String& placement, const G4String& name);

};

#endif


