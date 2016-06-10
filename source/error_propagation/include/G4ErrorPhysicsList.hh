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
// $Id: G4ErrorPhysicsList.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
// Class Description:
//
//  Default physics list for GEANT4e (should not be overridden, unless by
//  experts). No multiple scattering and no production of secondaries.
//  The energy loss process is G4eMuIonisation or G4EnergyLossForExtrapolator
//  (depending on the value of the enviromental variable G4EELOSSEXTRAP)
//  It also defines the geant4e processes to limit the step:
//  G4eMagneticFieldLimitProcess, G4eStepLimitProcess.

// History:
// - Created:   P. Arce
// ---------------------------------------------------------------------

#ifndef G4ErrorPhysicsList_hh
#define G4ErrorPhysicsList_hh

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class G4ErrorPhysicsList: public G4VUserPhysicsList
{
 public:  // with description

  G4ErrorPhysicsList();
  virtual ~G4ErrorPhysicsList();
  
 protected:

  virtual void ConstructParticle();
    // constructs gamma, e+/-, mu+/- and stable hadrons

  virtual void ConstructProcess();
    // construct physical processes

  virtual void SetCuts();  
    // SetCutsWithDefault

  virtual void ConstructEM();
    // constructs electromagnetic processes
};

#endif
