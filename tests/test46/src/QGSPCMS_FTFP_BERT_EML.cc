#include "QGSPCMS_FTFP_BERT_EML.hh"
#include "CMSEmStandardPhysics92.hh"
#include "G4EmStandardPhysics_option1.hh"

#include "G4SystemOfUnits.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4DataQuestionaire.hh"
#include "HadronPhysicsQGSP_FTFP_BERT.hh"

QGSPCMS_FTFP_BERT_EML::QGSPCMS_FTFP_BERT_EML()
{
  int ver = 1;

  G4DataQuestionaire it(photon);
  defaultCutValue = 0.7*mm;  
  SetVerboseLevel(ver);
  
  G4cout << "You are using " << "QGSP_FTFP_BERT_EML "
	 << G4endl;

  // EM Physics
  //RegisterPhysics( new CMSEmStandardPhysics92(ver));
  RegisterPhysics( new G4EmStandardPhysics_option1(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver));

  // Hadron Physics
  //G4bool quasiElastic=true;
  //RegisterPhysics( new HadronPhysicsQGSPCMS_FTFP_BERT("hadron",quasiElastic)); 
  RegisterPhysics( new HadronPhysicsQGSP_FTFP_BERT(ver)); 
  
  // Stopping Physics
  RegisterPhysics( new G4QStoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));
}

QGSPCMS_FTFP_BERT_EML::~QGSPCMS_FTFP_BERT_EML()
{}

void QGSPCMS_FTFP_BERT_EML::SetCuts()
{
  SetCutsWithDefault();   
}
