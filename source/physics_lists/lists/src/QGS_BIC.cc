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
// ClassName:   QGS_BIC
//
// Author: 2007 Gunter Folger
//     created from QGSP_BIC  by H.P.Wellisch
//
// Modified:
// 04.06.2010 G.Folger: Use new ctor for builders
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 16.10.2012 A.Ribon: Use new default stopping and ion physics
// 06.08.2019 A.Ribon: Use the new G4StoppingPhysicsFritiofWithBinaryCascade
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysicsFritiofWithBinaryCascade.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "QGS_BIC.hh"
#include "G4HadronPhysicsQGS_BIC.hh"

#include "G4WarnPLStatus.hh"

QGS_BIC::QGS_BIC(G4int ver)
{
  if(ver > 0) {
    G4cout << "<<< Geant4 Physics List simulation engine: QGS_BIC"<<G4endl;
    G4cout <<G4endl;
    G4WarnPLStatus exp;
    exp.Experimental("QGS_BIC");
  }
  defaultCutValue = 0.7*CLHEP::mm;  
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver) );

   // Hadron Physics
  RegisterPhysics(  new G4HadronPhysicsQGS_BIC(ver) );

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysicsFritiofWithBinaryCascade(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver) );
  
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver) );

  // Neutron cross sections
  //AR-31Oct2012 RegisterPhysics( new G4NeutronCrossSectionXS(ver) );

}


