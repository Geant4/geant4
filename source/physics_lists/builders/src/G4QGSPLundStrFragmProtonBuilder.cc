//
//---------------------------------------------------------------------------
//
// ClassName:   QGSPStrFragmLundProtonBuilder
//
// Author: Julia Yarba, FNAL/CD (2014)
//
//
//----------------------------------------------------------------------------
//
#include "G4QGSPLundStrFragmProtonBuilder.hh"
//
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4BGGNucleonInelasticXS.hh"

G4QGSPLundStrFragmProtonBuilder::
G4QGSPLundStrFragmProtonBuilder(G4bool quasiElastic) 
 {
   theMin = 100*GeV;
   theModel = new G4TheoFSGenerator("QGSP");

   theStringModel = new G4QGSModel< G4QGSParticipants >;
   theStringDecay = new G4ExcitedStringDecay(theStrFragm = new G4LundStringFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);

   theCascade = new G4GeneratorPrecompoundInterface;

   theModel->SetTransport(theCascade);
   theModel->SetHighEnergyGenerator(theStringModel);
   if (quasiElastic)
   {
      theQuasiElastic=new G4QuasiElasticChannel;
      theModel->SetQuasiElasticChannel(theQuasiElastic);
   } else 
   {  theQuasiElastic=0;}  

 }

void G4QGSPLundStrFragmProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
 {
   aP->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(100*TeV);
   aP->RegisterMe(theModel);
 }

void G4QGSPLundStrFragmProtonBuilder::
Build(G4HadronElasticProcess * )
 {
 }

G4QGSPLundStrFragmProtonBuilder::~G4QGSPLundStrFragmProtonBuilder() 
 {
   if ( theQuasiElastic ) delete theQuasiElastic;
   delete theStringDecay;
   delete theStringModel;
   delete theStrFragm;
 }

 
