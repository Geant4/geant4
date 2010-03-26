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
// $Id: G4CrossSectionDataTest.cc,v 1.20 2010-03-26 11:11:25 grichine Exp $
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

#include "G4NeutronHPInelasticData.hh"

// #include "G4HadronCrossSectionPlugin.hh"

#include "G4GlauberGribovCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4TripathiCrossSection.hh"

// #include "G4VQCrossSection.hh"
// #include "G4QElasticCrossSection.hh"
// #include "G4QuasiFreeRatios.hh"



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
  G4cout << "11 cadmium" << G4endl;
  G4cout << "12 uranium" << G4endl;
  G4cout << "13 berillium" << G4endl;
  G4cout << "14 aluminium" << G4endl;
  G4cout << "15 helium" << G4endl;
  G4int choice;
  // G4cin >> choice;
  choice = 6;



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
  choice = 1;

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
      theProcess = new G4PionMinusInelasticProcess("Inelastic");
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

  //  G4HadronCrossSectionPlugin theCrossSectionPlugin;

  // theCrossSectionPlugin.SetVerboseLevel(2);

  G4double sig = 0;   // , mfp = 0;

  switch (choice) 
  {
      case 1:

        theProcess = new G4HadronElasticProcess("ELASTIC");

      // ((G4HadronElasticProcess*)theProcess)->
      //   SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theElasticDataSet);

        // theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

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

        // theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      //mfp = ((G4HadronCaptureProcess*)theProcess)->
      //       GetMeanFreePathBasic(theDynamicParticle, theMaterial);
        break;

      case 4:
      //      theProcess = new G4HadronInelasticProcess("Inelastic",
      //                                                theParticleDefinition);

      // ((G4HadronInelasticProcess*)theProcess)->
      //    SetCrossSectionDataStore(&theCrossSectionDataStore);

        theCrossSectionDataStore.AddDataSet(&theInelasticDataSet);

        // theCrossSectionDataStore.AddDataSet(&theCrossSectionPlugin);

      //      mfp = ((G4HadronInelasticProcess*)theProcess)->
      //               GetMeanFreePathBasic(theDynamicParticle, theMaterial);
      break;
  }

     //   sig = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,
     //                                                  theElement);

  std::ofstream writef("g4txs.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );
  G4double ratio = 1.;

  /*

  G4ProtonInelasticCrossSection hpwPrIn;
  G4NeutronInelasticCrossSection hpwNeIn;
  G4PiNuclearCrossSection barIn;
  G4NucleonNuclearCrossSection barNucIn;

  G4NeutronHPInelasticData  nhpXscData;

  std::pair<G4double,G4double> chipsRat;

  // G4VQCrossSection*  qElastic = G4QElasticCrossSection::GetPointer();
  // G4QuasiFreeRatios* qRatio = G4QuasiFreeRatios::GetPointer();

  G4int pPDG = theParticleDefinition->GetPDGEncoding();
  G4int iz   = G4int(theElement->GetZ()); 
  G4int N    = G4int(theElement->GetN()+0.5);
        N   -= iz;
  if( N < 0 ) N  = 0;
  G4bool boolChips = true; 

                                  
  kinEnergy = 0.01*MeV;

  G4double barashXsc, hpwXsc, geishaXsc, neutronhpXsc;

  iMax = 90;
    
  writef <<iMax<< G4endl; 
  for(i = 0; i < iMax; i++)
  {
   
   theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);
    geishaXsc = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // G4cout << "GHEISHA" <<" \t"<< sig/millibarn << " mb" << G4endl;


    G4double momentum = theDynamicParticle->GetTotalMomentum();
    // if(i > 0) boolChips = false;
    // sig = qElastic->GetCrossSection(boolChips,momentum,iz,N,pPDG);

    // chipsRat = qRatio->GetRatios(momentum,pPDG,iz,N);
    // ratio = chipsRat.first*chipsRat.second; // !!!

    // G4cout << "Q-elastic          " <<" \t"<< sig/millibarn << " mb" << G4endl;


    // sig = hpwPrIn.GetCrossSection(theDynamicParticle, theElement, 273*kelvin);
    hpwXsc = hpwNeIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // G4cout << kinEnergy/GeV <<" GeV, \t"<< sig/millibarn << " mb" << G4endl;

    barashXsc = barNucIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    neutronhpXsc = nhpXscData.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // sig = barNucIn.GetTotalXsc();
    // sig = barNucIn.GetElasticXsc();

    // sig = barIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // sig = barIn.GetTotalXsc();
    // sig = barIn.GetElasticXsc();

    // G4cout << theProcess->GetProcessName() << " cross section for " << 
    //        theParticleDefinition->GetParticleName() <<
    //    " on " << theElement->GetName() << G4endl;
    // G4cout << kinEnergy/GeV <<" GeV, \t"<< sig/millibarn << " mb" << G4endl;
     //  G4cout << "Mean free path = " << mfp << " mm" << G4endl;

G4cout << kinEnergy/MeV <<" MeV, \t" << geishaXsc/millibarn << " mb \t"
       << hpwXsc/millibarn << " mb \t"<< barashXsc/millibarn << " mb \t"
       << neutronhpXsc/millibarn << " mb \t" << G4endl;      
writef  << kinEnergy/MeV <<" MeV, \t" << geishaXsc/millibarn << " mb \t"
       << hpwXsc/millibarn << " mb \t"<< barashXsc/millibarn << " mb \t"
       << neutronhpXsc/millibarn << " mb \t" << G4endl;      


    // writef << kinEnergy/GeV <<"\t"<< sig/millibarn << G4endl;


    // G4cout << kinEnergy/GeV <<" GeV, \t"<< "chips qe/in = "<< ratio << G4endl;
    //  writef << kinEnergy/GeV <<"\t"<< ratio <<G4endl;

    kinEnergy *= 1.138;
    delete theDynamicParticle;
  }
  G4cout<<"iz = "<<iz<<";  N = "<<N<<"; sum = "<<iz+N<<G4endl;                         
  
*/  
  // Check Glauber-Gribov model
  G4cout<<"Check Glauber-Gribov model"<<G4endl;

  std::ofstream writegg("ggtxs.dat", std::ios::out ) ;
  writegg.setf( std::ios::scientific, std::ios::floatfield );

  G4GlauberGribovCrossSection ggXsc;
  G4double ggTotXsc, ggElaXsc, ggIneXsc, ggProdXsc, ggDifXsc;

  
                  
  kinEnergy = 100.*GeV;
  iMax = 80;
  // writegg <<iMax<< G4endl; 
  writef <<iMax<< G4endl; 
  for(i = 0; i < iMax; i++)
  {
   
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);
     ggTotXsc = ggXsc.GetCrossSection(theDynamicParticle,
                                 theElement, 273*kelvin);

     ggElaXsc = ggXsc.GetElasticGlauberGribovXsc();
     // ggElaXsc = ggXsc.GetHadronNucleaonXscPDG(theDynamicParticle, theElement);

     ggIneXsc  = ggXsc.GetInelasticGlauberGribovXsc();
     ggProdXsc = ggXsc.GetProductionGlauberGribovXsc();
     ggDifXsc  = ggXsc.GetDiffractionGlauberGribovXsc();

     // if(ggIneXsc != 0.) ratio = 1 - ggProdXsc/ggIneXsc;
     // if(ratio < 0.) ratio = 0.;
     // if(ggIneXsc != 0.) ratio = ggDifXsc/ggIneXsc;

     // ratio = ggXsc.GetRatioSD(theDynamicParticle, theElement->GetN(), theElement->GetZ() );
     ratio = ggXsc.GetRatioQE(theDynamicParticle, theElement->GetN(), theElement->GetZ() );

     // G4cout << kinEnergy/GeV <<" GeV, total Xsc = "<<  ggTotXsc/millibarn 
     // << " mb; elastic Xsc = "
     // <<  ggElaXsc/millibarn << " mb; inelastic Xsc =" <<  ggIneXsc/millibarn << G4endl;
    
     // G4cout << kinEnergy/GeV <<" GeV, in Xsc = "<<  ggIneXsc/millibarn << " mb; prod Xsc = "
     //       <<  ggProdXsc/millibarn << " mb; ratio = 1 - prod/in = " <<  ratio << G4endl;
    
      G4cout << kinEnergy/GeV <<" GeV, in-Xsc = "<<  ggIneXsc/millibarn << " mb; sdif-Xsc = "
            <<  ggDifXsc/millibarn << " mb; ratio = dif/in = " <<  ratio << G4endl;
    

  // writegg << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
     // writef << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
  writef << kinEnergy/GeV <<"\t"<< ggIneXsc/millibarn <<"\t"<< ggProdXsc/millibarn << G4endl;
  //  writef << kinEnergy/GeV <<"\t"<< ratio <<G4endl;

     kinEnergy *= 1.145;
     delete theDynamicParticle;
  }
              
      

 
  /*
  G4double At = 1.;

  for(i = 0; i < iMax; i++)
  {
    G4double nuclRadius = ggXsc.GetNucleusRadius(At); 
    G4cout << "At = "<<At<<"      R = "<<nuclRadius/fermi<<";   " 
           << ggXsc.GetRadiusConst()*std::pow(At, 1./3.)/fermi<<G4endl;
    At *= 1.1;   
  }
  */

  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << theProcess->GetProcessName() << " cross section for " << 
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;


  return 1;
} // end of main
