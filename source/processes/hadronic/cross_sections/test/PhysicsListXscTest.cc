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
// $Id: PhysicsListXscTest.cc,v 1.2 2007-03-15 16:00:44 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test of G4CrossSectionDataStore and G4CrossSectionDataSet classes
//
// History:
//
// 14.03.07 V.Grichine: creation with PLs and NIST elements/materials, write in file
//



#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4HadronicProcess.hh"
#include "G4HadronElasticProcess.hh"
//#include "G4HadronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4UHadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"

#include "G4CrossSectionDataStore.hh"
#include "G4HadronCrossSections.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4HadronFissionDataSet.hh"
#include "G4HadronCaptureDataSet.hh"
#include "G4HadronInelasticDataSet.hh"

// #include "G4HadronCrossSectionPlugin.hh"

#include "G4GlauberGribovCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4TripathiCrossSection.hh"

#include "G4VQCrossSection.hh"
#include "G4QElasticCrossSection.hh"


// #include "G4ProcessVectorTypeIndex.hh"


#include "G4VUserPhysicsList.hh"
#include "FTFC.hh"
#include "FTFP.hh"

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronInelasticQLHEP.hh"

#include "HadronPhysicsFTFC.hh"
#include "HadronPhysicsFTFP.hh"

#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_BERT_HP.hh"
#include "HadronPhysicsLHEP_BIC.hh"
#include "HadronPhysicsLHEP_BIC_HP.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_HP.hh"
#include "HadronPhysicsLHEP_LEAD.hh"
#include "HadronPhysicsLHEP_LEAD_HP.hh"
#include "HadronPhysicsLHEP_PRECO.hh"
#include "HadronPhysicsLHEP_PRECO_HP.hh"

#include "HadronPhysicsQGSC_EFLOW.hh"
#include "HadronPhysicsQGSC_EMV.hh"
#include "HadronPhysicsQGSC.hh"
#include "HadronPhysicsQGSC_LEAD.hh"
#include "HadronPhysicsQGSC_LEAD_HP.hh"

#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"
#include "HadronPhysicsQGSP_EMV.hh"
#include "HadronPhysicsQGSP_EMX.hh"
#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_HP.hh"

#include "LBE.hh"

#include "LHEP_BERT.hh"
#include "LHEP_BERT_HP.hh"
#include "LHEP_BIC.hh"
#include "LHEP_BIC_HP.hh"
#include "LHEP_EMV.hh"
#include "LHEP.hh"
#include "LHEP_HP.hh"
#include "LHEP_LEAD.hh"
#include "LHEP_LEAD_HP.hh"
#include "LHEP_PRECO.hh"
#include "LHEP_PRECO_HP.hh"

#include "QBBC.hh"

#include "QGSC_EFLOW.hh"
#include "QGSC_EMV.hh"
#include "QGSC.hh"
#include "QGSC_LEAD.hh"
#include "QGSC_LEAD_HP.hh"

#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_EMV.hh"
#include "QGSP_EMX.hh"
#include "QGSP.hh"
#include "QGSP_HP.hh"
#include "QGSP_QEL.hh"






#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif


int main()
{

  // Physics list selection 

  G4VUserPhysicsList* pList; 


  // G4String mylist; 
  // mylist = "ftfc";

  // if (mylist == "ftfc") pList = new  FTFC();
  // else ;


  // pList = new  FTFP();
  // pList = new  LHEP();
  // pList = new  QGSP();
  pList = new  QGSC();

  pList->ConstructParticle();
  // pList->ConstructProcess();
  pList->Construct();



  // Element definition

  G4cout << " 1 copper" << G4endl;
  G4cout << " 2 iron" << G4endl;
  G4cout << " 3 lead" << G4endl;
  G4cout << " 4 hydrogen" << G4endl;
  G4cout << " 5 oxigen" << G4endl;
  G4cout << " 6 carbon" << G4endl;
  G4cout << " 7 nitrogen" << G4endl;
  G4cout << " 8 argon" << G4endl;
  G4cout << " 9 silicon" << G4endl;
  G4cout << "10 tugnsten" << G4endl;
  G4cout << "11 cadmium" << G4endl;
  G4cout << "12 uranium" << G4endl;
  G4cout << "13 berillium" << G4endl;
  G4cout << "14 aluminium" << G4endl;
  G4cout << "15 helium" << G4endl;
  G4int choice;
  // G4cin >> choice;
  choice = 1;



  G4Element*     theElement;
  G4Material*    theMaterial;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  switch (choice) 
  {
    case 1:

      theElement  = man->FindOrBuildElement("Cu");
      theMaterial = man->FindOrBuildMaterial("G4_Cu");
      break;

    case 2:

      theElement  = man->FindOrBuildElement("Fe");
      theMaterial = man->FindOrBuildMaterial("G4_Fe");
      break;

    case 3:

      theElement  = man->FindOrBuildElement("Pb");
      theMaterial = man->FindOrBuildMaterial("G4_Pb");
      break;

    case 4:

      theElement  = man->FindOrBuildElement("H");
      theMaterial = man->FindOrBuildMaterial("G4_H");
      break;

    case 5:

      theElement  = man->FindOrBuildElement("O");
      theMaterial = man->FindOrBuildMaterial("G4_O");
      break;

    case 6:

      theElement  = man->FindOrBuildElement("C");
      theMaterial = man->FindOrBuildMaterial("G4_C");
      break;

    case 7:

      theElement  = man->FindOrBuildElement("N");
      theMaterial = man->FindOrBuildMaterial("G4_N");
      break;

    case 8:

      theElement  = man->FindOrBuildElement("Ar");
      theMaterial = man->FindOrBuildMaterial("G4_Ar");
      break;

    case 9:

      theElement  = man->FindOrBuildElement("Si");
      theMaterial = man->FindOrBuildMaterial("G4_Si");
      break;

    case 10:

      theElement  = man->FindOrBuildElement("W");
      theMaterial = man->FindOrBuildMaterial("G4_W");
      break;

    case 11:

      theElement  = man->FindOrBuildElement("Cd");
      theMaterial = man->FindOrBuildMaterial("G4_Cd");
      break; 
 
    case 12:

      theElement  = man->FindOrBuildElement("U");
      theMaterial = man->FindOrBuildMaterial("G4_U");
      break;
  
    case 13:

      theElement  = man->FindOrBuildElement("Be");
      theMaterial = man->FindOrBuildMaterial("G4_Be");
      break; 
 
    case 14:

      theElement  = man->FindOrBuildElement("Al");
      theMaterial = man->FindOrBuildMaterial("G4_Al");
      break;
  
    case 15:

      theElement  = man->FindOrBuildElement("He");
      theMaterial = man->FindOrBuildMaterial("G4_He");
      break;  
  }
   //   G4cout << "Dumping element info:" << G4endl;
   //   theElement->DumpInfo();
   //   G4cout << "Dumping material info:" << G4endl;
   //   theMaterial->DumpInfo();

// Particle definition

  G4cout << " 1 proton" << G4endl;
  G4cout << " 2 pion+" << G4endl;
  G4cout << " 3 neutron" << G4endl;
  G4cout << " 4 kaon+" << G4endl;
  G4cout << " 5 kaon0short" << G4endl;
  G4cout << " 6 pion-" << G4endl;
  //  G4cin >> choice;
  choice = 3;

  G4ParticleDefinition* theParticleDefinition;
  G4VProcess* theProcess;

  switch (choice) 
  {
    case 1:

      theParticleDefinition = G4Proton::ProtonDefinition();
     
      break;

    case 2:

      theParticleDefinition = G4PionPlus::PionPlusDefinition();
    
      break;

    case 3:

      theParticleDefinition = G4Neutron::NeutronDefinition();
      
      break;

    case 4:

      theParticleDefinition = G4KaonPlus::KaonPlusDefinition();
     
      break;

    case 5:

      theParticleDefinition = G4KaonZeroShort::KaonZeroShortDefinition();
      
      break;

    case 6:

      theParticleDefinition = G4PionMinus::PionMinusDefinition();
    
      break;
  }
   //   G4cout << "Dumping particle info:" << G4endl;
   //   theParticleDefinition->DumpTable();

  // const G4ParticleDefinition partDef = *theParticleDefinition; 

  G4int i, iMax = 70, iProcess, iMaxProcess, iElastic=0, iInelastic=0;
  G4double kinEnergy;
  G4String processName, pn; //  = "hElastic";
  // G4String processName = "NeutronInelastic";



  // G4ProcessVectorTypeIndex typ = 1;

  iMaxProcess = theParticleDefinition->GetProcessManager()->GetProcessListLength();

// Make a dynamic particle too
  G4DynamicParticle* theDynamicParticle;

// Process definition

  G4cout << " 1 elastic" << G4endl;
  G4cout << " 2 fission" << G4endl;
  G4cout << " 3 capture" << G4endl;
  G4cout << " 4 inelastic" << G4endl;

  // G4cin >> choice;

  choice = 1;

  G4double sig = 0;   

  switch (choice) 
  {
      case 1:

        for (i = 0; i < iMaxProcess; i++) 
        {
          pn = (*theParticleDefinition->GetProcessManager()->
		   GetPostStepProcessVector(typeDoIt))[i]->GetProcessName();

          G4cout << i << "\t" << processName << G4endl;

          if( pn  == "hElastic" || pn == "CHIPSElasticScattering" ) 
          {
            iElastic = i;
            theProcess = (*theParticleDefinition->GetProcessManager()->
		     GetPostStepProcessVector(typeDoIt))[i];
             processName = (*theParticleDefinition->GetProcessManager()->
		   GetPostStepProcessVector(typeDoIt))[i]->GetProcessName();
          }    
        }
        G4cout<<"iElastic = "<<iElastic <<G4endl;
        G4cout<<"processName = "<<processName <<G4endl;
	iProcess = iElastic;

      break;

      case 2:

      break;


      case 3:
        break;

      case 4:
	iProcess = iInelastic;
      break;
  }
  G4HadronicProcess* hadProcess = dynamic_cast<G4HadronicProcess*>(theProcess);

  std::ofstream writef("plxsc.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  G4VQCrossSection* qElastic;
  G4int pPDG, iz, N;
  G4bool boolChips=false;

  if( processName != "CHIPSElasticScattering") 
  {
    hadProcess->BuildPhysicsTable(*theParticleDefinition);
  }
  else
  {
    qElastic = G4QElasticCrossSection::GetPointer();

    pPDG = theParticleDefinition->GetPDGEncoding();
    iz   = G4int(theElement->GetZ()); 
    N    = G4int(theElement->GetN()+0.5);
    N   -= iz;
    if( N < 0 ) N  = 0;
    boolChips = true;                              
  }
  iMax = 90;
    
  writef <<iMax<< G4endl; 

  kinEnergy = 0.01*GeV;

  for(i = 0; i < iMax; i++)
  {      
    theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);

    if( processName != "CHIPSElasticScattering")
    {
      sig = hadProcess->GetMicroscopicCrossSection(theDynamicParticle,
                                                 theElement,273*kelvin);
      G4cout << kinEnergy/GeV <<" GeV, \t" << "PL elastic" <<" \t"
           << sig/millibarn << " mb" << G4endl;
    }
    else
    {
      if(i > 0) boolChips = false;    
      G4double momentum = theDynamicParticle->GetTotalMomentum();
      sig = qElastic->GetCrossSection(boolChips,momentum,iz,N,pPDG);
      G4cout << "Q-elastic          " <<" \t"<< sig/millibarn << " mb" << G4endl;
    }

    writef << kinEnergy/GeV <<"\t"<< sig/millibarn << G4endl;

    kinEnergy *= 1.138;
    delete theDynamicParticle;
  }
  G4cout << theProcess->GetProcessName() << " cross section for " << 
           theParticleDefinition->GetParticleName() <<
          " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;
                         

  return 1;
} // end of main
