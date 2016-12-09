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
// $Id: G4HadronInelasticQBBC.cc 99977 2016-10-13 07:26:42Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 2 October 2009 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQBBC.hh"

#include "G4SystemOfUnits.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronicInteraction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4PiNuclearCrossSection.hh"

#include "G4QGSBuilder.hh"
#include "G4FTFBuilder.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ChipsHyperonInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4NeutronRadCapture.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4HadronicInteractionRegistry.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronInelasticQBBC);

G4ThreadLocal G4ComponentAntiNuclNuclearXS* G4HadronInelasticQBBC::theAntiNuclXS = nullptr;
G4ThreadLocal G4ComponentGGHadronNucleusXsc* G4HadronInelasticQBBC::theKaonXS = nullptr;

G4HadronInelasticQBBC::G4HadronInelasticQBBC(G4int ver) 
  : G4VHadronPhysics("hInelasticQBBC"),verbose(ver)
{
  theAntiNuclXS = nullptr;
  theKaonXS = nullptr;
}

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String&, G4int ver, 
    G4bool, G4bool,G4bool, G4bool, G4bool) : G4HadronInelasticQBBC(ver)
{}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{
  delete theAntiNuclXS; theAntiNuclXS=nullptr;
  delete theKaonXS;     theKaonXS=nullptr;
}

void G4HadronInelasticQBBC::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### HadronInelasticQBBC Construct Process " << G4endl;
  }

  G4double emax = 100.*TeV;

  //G4cout << "G4HadronInelasticQBBC::ConstructProcess new PRECO"<< G4endl;

  // PreCompound and Evaporation models are instantiated here
  G4PreCompoundModel* thePreCompound = nullptr;
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  thePreCompound = static_cast<G4PreCompoundModel*>(p);
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel(); }
 
  // configure models
  //G4HadronicInteraction* theQGSP = 
  //  BuildModel(new G4QGSBuilder("QGSP",thePreCompound,true,false),12.5*GeV,emax);
  G4HadronicInteraction* theFTFP = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),3.0*GeV,emax);
  G4HadronicInteraction* theFTFP1 = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),3.0*GeV,emax);
  G4HadronicInteraction* theFTFP2 = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),0.0,emax);

  G4CascadeInterface* casc = new G4CascadeInterface();
  casc->usePreCompoundDeexcitation();
  G4HadronicInteraction* theBERT = NewModel(casc,1.0*GeV,5.0*GeV);

  casc = new G4CascadeInterface();
  casc->usePreCompoundDeexcitation();
  G4HadronicInteraction* theBERT1 = NewModel(casc,0.0*GeV,5.0*GeV);

  //G4cout << "G4HadronInelasticQBBC::ConstructProcess new Binary"<< G4endl;
  G4BinaryCascade* bic = new G4BinaryCascade(thePreCompound);
  G4HadronicInteraction* theBIC = NewModel(bic,0.0,1.5*GeV);

  // cross sections
  theAntiNuclXS = new G4ComponentAntiNuclNuclearXS();
  G4CrossSectionInelastic* anucxs = 
    new G4CrossSectionInelastic(theAntiNuclXS);
  theKaonXS = new G4ComponentGGHadronNucleusXsc();
  G4CrossSectionInelastic* kaonxs = new G4CrossSectionInelastic(theKaonXS);

  // loop over particles
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ) {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String pname = particle->GetParticleName();
    if(verbose > 1) { 
      G4cout << "### HadronInelasticQBBC:  " << pname << G4endl; 
    }

    //
    // model and X-section configuration per particle type
    //
    if(pname == "proton") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet(new G4BGGNucleonInelasticXS(particle));
      
      //hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
      hp->RegisterMe(theBERT);
      hp->RegisterMe(theBIC);

    } else if(pname == "neutron") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet((G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name()));
      //hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
       
      G4HadronicProcess* capture = FindCaptureProcess();
      capture->AddDataSet((G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name()));
      hp->RegisterMe(theBERT);
      hp->RegisterMe(theBIC);
      capture->RegisterMe(new G4NeutronRadCapture());

    } else if(pname == "pi-" || pname == "pi+") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->AddDataSet(new G4BGGPionInelasticXS(particle));
      //hp->AddDataSet(new G4CrossSectionPairGG((G4PiNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4PiNuclearCrossSection::Default_Name()), 91*GeV));
      //hp->RegisterMe(theQGSP);
      hp->RegisterMe(theFTFP);
      hp->RegisterMe(theBERT1);

    } else if(pname == "kaon-" ) {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->RegisterMe(theFTFP1);
      hp->RegisterMe(theBERT1);
      hp->AddDataSet(kaonxs);

    } else if(pname == "kaon+" ) {
        G4HadronicProcess* hp = FindInelasticProcess(particle);
        hp->RegisterMe(theFTFP1);
        hp->RegisterMe(theBERT1);
        hp->AddDataSet(kaonxs);

    } else if(pname == "kaon0S"    ||
              pname == "kaon0L") {
        G4HadronicProcess* hp = FindInelasticProcess(particle);
        hp->RegisterMe(theFTFP1);
        hp->RegisterMe(theBERT1);
        hp->AddDataSet(kaonxs);
        
    } else if(pname == "lambda"    ||
              pname == "omega-"    ||
	      pname == "sigma-"    || 
	      pname == "sigma+"    || 
	      pname == "sigma0"    || 
	      pname == "xi-"       || 
	      pname == "xi0") {
      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->RegisterMe(theFTFP1);
      hp->RegisterMe(theBERT1);
      hp->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsHyperonInelasticXS::Default_Name()));

    } else if(pname == "anti_alpha"   ||
	      pname == "anti_deuteron"||
              pname == "anti_He3"     ||
	      pname == "anti_proton"  || 
              pname == "anti_triton"  ||  
	      pname == "anti_lambda"  ||
              pname == "anti_neutron" ||
              pname == "anti_omega-"  || 
              pname == "anti_sigma-"  || 
              pname == "anti_sigma+"  || 
              pname == "anti_xi-"     || 
              pname == "anti_xi0"     
	      ) {

      G4HadronicProcess* hp = FindInelasticProcess(particle);
      hp->RegisterMe(theFTFP2);
      hp->AddDataSet(anucxs);

    } 
  }
}
