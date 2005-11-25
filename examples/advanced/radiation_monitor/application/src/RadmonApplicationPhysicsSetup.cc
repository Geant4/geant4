//
// File name:     RadmonApplicationPhysicsSetup.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.cc,v 1.4 2005-11-25 11:54:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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

 return true;
}
