//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CrossSectionDataTest.cc,v 1.7 2006-06-07 13:53:12 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test of G4CrossSectionDataStore and G4CrossSectionDataSet classes
//
// History:
//
// 29.05.06 V.Grichine: NIST elements/materials, write in file
// F.W. Jones, TRIUMF, 22-JAN-98
//                     19-MAY-98
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

#include "G4HadronElasticProcess.hh"
//#include "G4HadronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

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

#include "G4HadronCrossSectionPlugin.hh"

#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif


int main()
{
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
  G4cout << "11 uranium" << G4endl;
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

      theElement  = man->FindOrBuildElement("U");
      theMaterial = man->FindOrBuildMaterial("G4_U");
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
      theProcess = new G4ProtonInelasticProcess("Inelastic");
      break;

    case 2:

      theParticleDefinition = G4PionPlus::PionPlusDefinition();
      theProcess = new G4PionPlusInelasticProcess("Inelastic");
      break;

    case 3:

      theParticleDefinition = G4Neutron::NeutronDefinition();
      theProcess = new G4NeutronInelasticProcess("Inelastic");
      break;

    case 4:

      theParticleDefinition = G4KaonPlus::KaonPlusDefinition();
      theProcess = new G4KaonPlusInelasticProcess("Inelastic");
      break;

    case 5:

      theParticleDefinition = G4KaonZeroShort::KaonZeroShortDefinition();
      theProcess = new G4KaonZeroSInelasticProcess("Inelastic");
      break;

    case 6:

      theParticleDefinition = G4PionMinus::PionMinusDefinition();
      break;
  }
   //   G4cout << "Dumping particle info:" << G4endl;
   //   theParticleDefinition->DumpTable();


  G4int i, iMax = 70;
  G4double kinEnergy;
 

 // G4cout << "Kinetic energy in GeV: "<<G4endl;
  // G4cin >> kinEnergy;

// Make a dynamic particle too
  G4DynamicParticle* theDynamicParticle;

// Process definition

  G4cout << " 1 elastic" << G4endl;
  G4cout << " 2 fission" << G4endl;
  G4cout << " 3 capture" << G4endl;
  G4cout << " 4 inelastic" << G4endl;

  // G4cin >> choice;
  choice = 4;

  /////////////////////////////

  G4CrossSectionDataStore  theCrossSectionDataStore;

  G4HadronElasticDataSet   theElasticDataSet;
  G4HadronCaptureDataSet   theCaptureDataSet;
  G4HadronInelasticDataSet theInelasticDataSet;


  G4HadronFissionDataSet   theFissionDataSet;
  theFissionDataSet.SetVerboseLevel(2);

  G4HadronCrossSectionPlugin theCrossSectionPlugin;

  // theCrossSectionPlugin.SetVerboseLevel(2);

  G4double sig = 0;   // , mfp = 0;

  switch (choice) 
  {
      case 1:

        theProcess = new G4HadronElasticProcess("ELASTIC");

      // ((G4HadronElasticProcess*)theProcess)->
      //   SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theElasticDataSet);

        theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      // mfp = ((G4HadronElasticProcess*)theProcess)->
      //         GetMeanFreePathBasic(theDynamicParticle, theMaterial);
      break;

      case 2:

        theProcess = new G4HadronFissionProcess("Fission");

      // ((G4HadronFissionProcess*)theProcess)->
      //    SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theFissionDataSet);

      //      theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      // mfp = ((G4HadronFissionProcess*)theProcess)->
      //         GetMeanFreePathBasic(theDynamicParticle, theMaterial);
      break;


      case 3:

        theProcess = new G4HadronCaptureProcess("Capture");

      // ((G4HadronCaptureProcess*)theProcess)->
      //    SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theCaptureDataSet);

        theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      //mfp = ((G4HadronCaptureProcess*)theProcess)->
      //       GetMeanFreePathBasic(theDynamicParticle, theMaterial);
        break;

      case 4:
      //      theProcess = new G4HadronInelasticProcess("Inelastic",
      //                                                theParticleDefinition);

      // ((G4HadronInelasticProcess*)theProcess)->
      //    SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theInelasticDataSet);

        theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      //      mfp = ((G4HadronInelasticProcess*)theProcess)->
      //               GetMeanFreePathBasic(theDynamicParticle, theMaterial);
      break;
  }

     //   sig = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,
     //                                                  theElement);

  std::ofstream writef("txs.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  kinEnergy = 500.*GeV; //1.

  // for(i = 0; i < iMax; i++)
  {
   
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);
     sig = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,
                                                    theElement, 273*kelvin);

     // G4cout << theProcess->GetProcessName() << " cross section for " << 
     //       theParticleDefinition->GetParticleName() <<
     //      " on " << theElement->GetName() << G4endl;
     G4cout << kinEnergy/GeV <<" GeV, \t"<< sig/millibarn << " mb" << G4endl;
     // G4cout << "Mean free path = " << mfp << " mm" << G4endl;

     // writef << kinEnergy/GeV <<"\t"<< sig/millibarn << G4endl;

     kinEnergy *= 1.1;
     delete theDynamicParticle;
  }
  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << theProcess->GetProcessName() << " cross section for " << 
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;

  return 1;
} // end of main
