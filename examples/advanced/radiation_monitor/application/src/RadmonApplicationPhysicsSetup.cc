//
// File name:     RadmonApplicationPhysicsSetup.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.cc,v 1.2 2005-11-10 08:17:16 capra Exp $
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

 return true;
}
