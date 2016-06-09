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
// File name:     RadmonApplicationPhysicsSetup.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.cc,v 1.4.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-03 $
//

// Include files
#include "RadmonApplicationPhysicsSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonSubPhysicsListWithLabelFactory.hh"

#include "RadmonPhysicsElectronEEDL.hh"
#include "RadmonPhysicsElectronStandard.hh"
#include "RadmonPhysicsPhotonEPDL.hh"
#include "RadmonPhysicsPhotonStandard.hh"
#include "RadmonPhysicsPositronStandard.hh"
#include "RadmonPhysicsMuonStandard.hh"
#include "RadmonPhysicsTauStandard.hh"
#include "RadmonPhysicsNuclear.hh"
#include "RadmonPhysicsDecay.hh"
#include "RadmonPhysicsNeutronBinary.hh"
#include "RadmonPhysicsNeutronBertini.hh"
#include "RadmonPhysicsHadronsBinary.hh"
#include "RadmonPhysicsHadronsBertini.hh"
#include "RadmonPhysicsICRUIonization.hh"
#include "RadmonPhysicsProductionCuts.hh"
#include "RadmonPhysicsParticles.hh"

#define DECLARE_SUBPHYSICS_LIST(name)           subPhysicsList=new name();                                                               \
                                                if (subPhysicsList==0)                                                                   \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendSubPhysicsListWithLabel(subPhysicsList)

G4bool                                          RadmonApplicationPhysicsSetup :: CreateSubPhysicsList(RadmonSubPhysicsListWithLabelFactory * factory)
{
 RadmonVSubPhysicsListWithLabel * subPhysicsList;
 
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsElectronEEDL);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsElectronStandard);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsPhotonEPDL);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsPhotonStandard);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsPositronStandard);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsMuonStandard);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsTauStandard);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsNuclear);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsDecay);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsNeutronBinary);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsNeutronBertini);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsHadronsBinary);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsHadronsBertini);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsICRUIonization);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsProductionCuts);
 DECLARE_SUBPHYSICS_LIST(RadmonPhysicsParticles);

 return true;
}
