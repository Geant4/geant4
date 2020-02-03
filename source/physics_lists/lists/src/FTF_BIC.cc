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
//---------------------------------------------------------------------------
//
// ClassName:   FTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
// 19.06.2008 G.Folger: don't use chips quasielastic in FTF
// 04.06.2010 G.Folger: Use new ctor for builders
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 16.10.2012 A.Ribon: Use new default stopping and ion physics
// 06.08.2019 A.Ribon: Use the new G4StoppingPhysicsFritiofWithBinaryCascade
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysicsFritiofWithBinaryCascade.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "FTF_BIC.hh"
#include "G4HadronPhysicsFTF_BIC.hh"

#include "G4WarnPLStatus.hh"

FTF_BIC::FTF_BIC(G4int ver)
{
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;

  G4cout << "<<< Geant4 Physics List simulation engine: FTF_BIC"<<G4endl;
  G4cout <<G4endl;

  defaultCutValue = 0.7*CLHEP::mm;  
  SetVerboseLevel(ver);

  G4WarnPLStatus exp;
  exp.Experimental("FTF_BIC");

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
   RegisterPhysics( new G4HadronElasticPhysics(ver) );

   // Hadron Physics
  RegisterPhysics(  new G4HadronPhysicsFTF_BIC(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysicsFritiofWithBinaryCascade(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));
  
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));

  // Neutron cross sections
  //AR-31Oct2012 RegisterPhysics( new G4NeutronCrossSectionXS(ver) );

}


