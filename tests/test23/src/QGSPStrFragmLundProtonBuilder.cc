//
//---------------------------------------------------------------------------
//
// ClassName:   QGSPStrFragmLundProtonBuilder
//
// Author: Julia Yarba, FNAL/CD (2013)
//
//
//----------------------------------------------------------------------------
//
//specific no NuMI-X
#include "QGSPStrFragmLundProtonBuilder.hh"
//
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4BGGNucleonInelasticXS.hh"

QGSPStrFragmLundProtonBuilder::
QGSPStrFragmLundProtonBuilder(G4bool quasiElastic) 
 {
   theMin = 100*GeV;
   theModel = new G4TheoFSGenerator("QGSP");

   theStringModel = new G4QGSModel< G4QGSParticipants >;
   theStringDecay = new G4ExcitedStringDecay(theStrFragm = new G4LundStringFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);

   theCascade = new G4GeneratorPrecompoundInterface;
   thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
   theCascade->SetDeExcitation(thePreEquilib);  

   theModel->SetTransport(theCascade);
   theModel->SetHighEnergyGenerator(theStringModel);
   if (quasiElastic)
   {
      theQuasiElastic=new G4QuasiElasticChannel;
      theModel->SetQuasiElasticChannel(theQuasiElastic);
   } else 
   {  theQuasiElastic=0;}  

 }

void QGSPStrFragmLundProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
 {
   aP->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(100*TeV);
   aP->RegisterMe(theModel);
 }

void QGSPStrFragmLundProtonBuilder::
Build(G4HadronElasticProcess * )
 {
 }

QGSPStrFragmLundProtonBuilder::~QGSPStrFragmLundProtonBuilder() 
 {
   delete thePreEquilib;
   delete theCascade;
   if ( theQuasiElastic ) delete theQuasiElastic;
   delete theStringDecay;
   delete theStringModel;
   delete theModel;
   delete theStrFragm;
   //delete theHandler;
 }

 
