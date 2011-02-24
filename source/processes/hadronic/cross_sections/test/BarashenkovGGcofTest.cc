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
// $Id: BarashenkovGGcofTest.cc,v 1.4 2008-09-18 15:58:52 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test for normalisation coefficient between Barashenkov and GG models
//
// History:
//
// 10.09.08 V.Grichine: NIST elements/materials, write in file
//
//
//
//
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

#include "G4VQCrossSection.hh"
#include "G4QElasticCrossSection.hh"
#include "G4QuasiFreeRatios.hh"



#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif


int main()
{

  G4int choice;
  G4int i, iMax, iz, N;

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
      theProcess = new G4PionMinusInelasticProcess("Inelastic");
      break;
  }
   //   G4cout << "Dumping particle info:" << G4endl;
   //   theParticleDefinition->DumpTable();


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


  G4ProtonInelasticCrossSection hpwPrIn;
  G4NeutronInelasticCrossSection hpwNeIn;


  G4PiNuclearCrossSection barPiIn;
  G4NucleonNuclearCrossSection barNucIn;

  G4NeutronHPInelasticData  nhpXscData;

  G4double ratio = 1.;
  std::pair<G4double,G4double> chipsRat;

  G4VQCrossSection*  qElastic = G4QElasticCrossSection::GetPointer();
  G4QuasiFreeRatios* qRatio = G4QuasiFreeRatios::GetPointer();

  G4int pPDG = theParticleDefinition->GetPDGEncoding();


  G4Element*     theElement;
  G4Material*    theMaterial;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
                                  
  kinEnergy = 90*GeV;

  G4double barashInXsc, barashTotXsc, barashElXsc, hpwXsc, geishaXsc, neutronhpXsc;

  std::ofstream writef("g4txs.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream writeCorTot("corBarTot.dat", std::ios::out ) ;
  writeCorTot.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream writeCorIn("corBarIn.dat", std::ios::out ) ;
  writeCorIn.setf( std::ios::scientific, std::ios::floatfield );


  std::ofstream writegg("ggtxs.dat", std::ios::out ) ;
  writegg.setf( std::ios::scientific, std::ios::floatfield );

  G4GlauberGribovCrossSection ggXsc;
  G4double ggTotXsc, ggElaXsc, ggIneXsc, ggProdXsc,ggDifXsc;
  G4double ratioTot, ratioIn, ratioEl;

  iMax = 93;
    
  writef <<iMax-2<< G4endl; 

  for(i = 2; i < iMax; i++)
  {
   
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);

    choice = i;

    switch (choice) 
    {
      case 1:

      theElement  = man->FindOrBuildElement("H");
      theMaterial = man->FindOrBuildMaterial("G4_H");
      break;

      case 2:

      theElement  = man->FindOrBuildElement("He");
      theMaterial = man->FindOrBuildMaterial("G4_He");
      break; 
 
      case 3:

      theElement  = man->FindOrBuildElement("Li");
      theMaterial = man->FindOrBuildMaterial("G4_Li");
      break; 
 
      case 4:

      theElement  = man->FindOrBuildElement("Be");
      theMaterial = man->FindOrBuildMaterial("G4_Be");
      break; 

      case 5:

      theElement  = man->FindOrBuildElement("B");
      theMaterial = man->FindOrBuildMaterial("G4_B");
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

      theElement  = man->FindOrBuildElement("O");
      theMaterial = man->FindOrBuildMaterial("G4_O");
      break;

      case 9:

      theElement  = man->FindOrBuildElement("F");
      theMaterial = man->FindOrBuildMaterial("G4_F");
      break;

      case 10:

      theElement  = man->FindOrBuildElement("Ne");
      theMaterial = man->FindOrBuildMaterial("G4_Ne");
      break;

      case 11:

      theElement  = man->FindOrBuildElement("Na");
      theMaterial = man->FindOrBuildMaterial("G4_Na");
      break;

      case 12:

      theElement  = man->FindOrBuildElement("Mg");
      theMaterial = man->FindOrBuildMaterial("G4_Mg");
      break;


      case 13:

      theElement  = man->FindOrBuildElement("Al");
      theMaterial = man->FindOrBuildMaterial("G4_Al");
      break;
  
      case 14:

      theElement  = man->FindOrBuildElement("Si");
      theMaterial = man->FindOrBuildMaterial("G4_Si");
      break;

      case 15:

      theElement  = man->FindOrBuildElement("P");
      theMaterial = man->FindOrBuildMaterial("G4_P");
      break;

      case 16:

      theElement  = man->FindOrBuildElement("S");
      theMaterial = man->FindOrBuildMaterial("G4_S");
      break;

      case 17:

      theElement  = man->FindOrBuildElement("Cl");
      theMaterial = man->FindOrBuildMaterial("G4_Cl");
      break;

      case 18:

      theElement  = man->FindOrBuildElement("Ar");
      theMaterial = man->FindOrBuildMaterial("G4_Ar");
      break;

      case 19:

      theElement  = man->FindOrBuildElement("K");
      theMaterial = man->FindOrBuildMaterial("G4_K");
      break;

      case 20:

      theElement  = man->FindOrBuildElement("Ca");
      theMaterial = man->FindOrBuildMaterial("G4_Ca");
      break;

      case 21:

      theElement  = man->FindOrBuildElement("Sc");
      theMaterial = man->FindOrBuildMaterial("G4_Sc");
      break;

      case 22:

      theElement  = man->FindOrBuildElement("Ti");
      theMaterial = man->FindOrBuildMaterial("G4_Ti");
      break;

      case 23:

      theElement  = man->FindOrBuildElement("V");
      theMaterial = man->FindOrBuildMaterial("G4_V");
      break;

      case 24:

      theElement  = man->FindOrBuildElement("Cr");
      theMaterial = man->FindOrBuildMaterial("G4_Cr");
      break;

      case 25:

      theElement  = man->FindOrBuildElement("Mn");
      theMaterial = man->FindOrBuildMaterial("G4_Mn");
      break;


      case 26:

      theElement  = man->FindOrBuildElement("Fe");
      theMaterial = man->FindOrBuildMaterial("G4_Fe");
      break;

      case 27:

      theElement  = man->FindOrBuildElement("Co");
      theMaterial = man->FindOrBuildMaterial("G4_Co");
      break;

      case 28:

      theElement  = man->FindOrBuildElement("Ni");
      theMaterial = man->FindOrBuildMaterial("G4_Ni");
      break;


      case 29:

      theElement  = man->FindOrBuildElement("Cu");
      theMaterial = man->FindOrBuildMaterial("G4_Cu");
      break;

      case 30:

      theElement  = man->FindOrBuildElement("Zn");
      theMaterial = man->FindOrBuildMaterial("G4_Zn");
      break;

      case 31:

      theElement  = man->FindOrBuildElement("Ga");
      theMaterial = man->FindOrBuildMaterial("G4_Ga");
      break;

      case 32:

      theElement  = man->FindOrBuildElement("Ge");
      theMaterial = man->FindOrBuildMaterial("G4_Ge");
      break;

      case 33:

      theElement  = man->FindOrBuildElement("As");
      theMaterial = man->FindOrBuildMaterial("G4_As");
      break;

      case 34:

      theElement  = man->FindOrBuildElement("Se");
      theMaterial = man->FindOrBuildMaterial("G4_Se");
      break;

      case 35:

      theElement  = man->FindOrBuildElement("Br");
      theMaterial = man->FindOrBuildMaterial("G4_Br");
      break;

    case 36:

      theElement  = man->FindOrBuildElement("Kr");
      theMaterial = man->FindOrBuildMaterial("G4_Kr");
      break;

      case 37:

      theElement  = man->FindOrBuildElement("Rb");
      theMaterial = man->FindOrBuildMaterial("G4_Rb");
      break;

      case 38:

      theElement  = man->FindOrBuildElement("Sr");
      theMaterial = man->FindOrBuildMaterial("G4_Sr");
      break;

      case 39:

      theElement  = man->FindOrBuildElement("Y");
      theMaterial = man->FindOrBuildMaterial("G4_Y");
      break;

      case 40:

      theElement  = man->FindOrBuildElement("Zr");
      theMaterial = man->FindOrBuildMaterial("G4_Zr");
      break;

      case 41:

      theElement  = man->FindOrBuildElement("Nb");
      theMaterial = man->FindOrBuildMaterial("G4_Nb");
      break;

      case 42:

      theElement  = man->FindOrBuildElement("Mo");
      theMaterial = man->FindOrBuildMaterial("G4_Mo");
      break;

      case 43:

      theElement  = man->FindOrBuildElement("Tc");
      theMaterial = man->FindOrBuildMaterial("G4_Tc");
      break;

      case 44:

      theElement  = man->FindOrBuildElement("Ru");
      theMaterial = man->FindOrBuildMaterial("G4_Ru");
      break;

      case 45:

      theElement  = man->FindOrBuildElement("Rh");
      theMaterial = man->FindOrBuildMaterial("G4_Rh");
      break;

    case 46:

      theElement  = man->FindOrBuildElement("Pd");
      theMaterial = man->FindOrBuildMaterial("G4_Pd");
      break;

      case 47:

      theElement  = man->FindOrBuildElement("Ag");
      theMaterial = man->FindOrBuildMaterial("G4_Ag");
      break;

      case 48:

      theElement  = man->FindOrBuildElement("Cd");
      theMaterial = man->FindOrBuildMaterial("G4_Cd");
      break; 
 
      case 49:

      theElement  = man->FindOrBuildElement("In");
      theMaterial = man->FindOrBuildMaterial("G4_In");
      break;

      case 50:

      theElement  = man->FindOrBuildElement("Sn");
      theMaterial = man->FindOrBuildMaterial("G4_Sn");
      break;

      case 51:

      theElement  = man->FindOrBuildElement("Sb");
      theMaterial = man->FindOrBuildMaterial("G4_Sb");
      break;
 
      case 52:

      theElement  = man->FindOrBuildElement("Te");
      theMaterial = man->FindOrBuildMaterial("G4_Te");
      break;

      case 53:

      theElement  = man->FindOrBuildElement("I");
      theMaterial = man->FindOrBuildMaterial("G4_I");
      break;

      case 54:

      theElement  = man->FindOrBuildElement("Xe");
      theMaterial = man->FindOrBuildMaterial("G4_Xe");
      break;

      case 55:

      theElement  = man->FindOrBuildElement("Cs");
      theMaterial = man->FindOrBuildMaterial("G4_Cs");
      break;

      case 56:

      theElement  = man->FindOrBuildElement("Ba");
      theMaterial = man->FindOrBuildMaterial("G4_Ba");
      break;

      case 57:

      theElement  = man->FindOrBuildElement("La");
      theMaterial = man->FindOrBuildMaterial("G4_La");
      break;

      case 58:

      theElement  = man->FindOrBuildElement("Ce");
      theMaterial = man->FindOrBuildMaterial("G4_Ce");
      break;

      case 59:

      theElement  = man->FindOrBuildElement("Pr");
      theMaterial = man->FindOrBuildMaterial("G4_Pr");
      break;

      case 60:

      theElement  = man->FindOrBuildElement("Nd");
      theMaterial = man->FindOrBuildMaterial("G4_Nd");
      break;

      case 61:

      theElement  = man->FindOrBuildElement("Pm");
      theMaterial = man->FindOrBuildMaterial("G4_Pm");
      break;

      case 62:

      theElement  = man->FindOrBuildElement("Sm");
      theMaterial = man->FindOrBuildMaterial("G4_Sm");
      break;

      case 63:

      theElement  = man->FindOrBuildElement("Eu");
      theMaterial = man->FindOrBuildMaterial("G4_Eu");
      break;

      case 64:

      theElement  = man->FindOrBuildElement("Gd");
      theMaterial = man->FindOrBuildMaterial("G4_Gd");
      break;

    case 65:

      theElement  = man->FindOrBuildElement("Tb");
      theMaterial = man->FindOrBuildMaterial("G4_Tb");
      break;

      case 66:

      theElement  = man->FindOrBuildElement("Dy");
      theMaterial = man->FindOrBuildMaterial("G4_Dy");
      break;

      case 67:

      theElement  = man->FindOrBuildElement("Ho");
      theMaterial = man->FindOrBuildMaterial("G4_Ho");
      break;

      case 68:

      theElement  = man->FindOrBuildElement("Er");
      theMaterial = man->FindOrBuildMaterial("G4_Er");
      break;

      case 69:

      theElement  = man->FindOrBuildElement("Tm");
      theMaterial = man->FindOrBuildMaterial("G4_Tm");
      break;

      case 70:

      theElement  = man->FindOrBuildElement("Yb");
      theMaterial = man->FindOrBuildMaterial("G4_Yb");
      break;

      case 71:

      theElement  = man->FindOrBuildElement("Lu");
      theMaterial = man->FindOrBuildMaterial("G4_Lu");
      break;

      case 72:

      theElement  = man->FindOrBuildElement("Hf");
      theMaterial = man->FindOrBuildMaterial("G4_Hf");
      break;

      case 73:

      theElement  = man->FindOrBuildElement("Ta");
      theMaterial = man->FindOrBuildMaterial("G4_Ta");
      break;

      case 74:

      theElement  = man->FindOrBuildElement("W");
      theMaterial = man->FindOrBuildMaterial("G4_W");
      break;

      case 75:

      theElement  = man->FindOrBuildElement("Re");
      theMaterial = man->FindOrBuildMaterial("G4_Re");
      break;

      case 76:

      theElement  = man->FindOrBuildElement("Os");
      theMaterial = man->FindOrBuildMaterial("G4_Os");
      break;

      case 77:

      theElement  = man->FindOrBuildElement("Ir");
      theMaterial = man->FindOrBuildMaterial("G4_Ir");
      break;

      case 78:

      theElement  = man->FindOrBuildElement("Pt");
      theMaterial = man->FindOrBuildMaterial("G4_Pt");
      break;

      case 79:

      theElement  = man->FindOrBuildElement("Au");
      theMaterial = man->FindOrBuildMaterial("G4_Au");
      break;

      case 80:

      theElement  = man->FindOrBuildElement("Hg");
      theMaterial = man->FindOrBuildMaterial("G4_Hg");
      break;

      case 81:

      theElement  = man->FindOrBuildElement("Tl");
      theMaterial = man->FindOrBuildMaterial("G4_Tl");
      break;

      case 82:

      theElement  = man->FindOrBuildElement("Pb");
      theMaterial = man->FindOrBuildMaterial("G4_Pb");
      break;

      case 83:

      theElement  = man->FindOrBuildElement("Bi");
      theMaterial = man->FindOrBuildMaterial("G4_Bi");
      break;

      case 84:

      theElement  = man->FindOrBuildElement("Po");
      theMaterial = man->FindOrBuildMaterial("G4_Po");
      break;

      case 85:

      theElement  = man->FindOrBuildElement("At");
      theMaterial = man->FindOrBuildMaterial("G4_At");
      break;

      case 86:

      theElement  = man->FindOrBuildElement("Rn");
      theMaterial = man->FindOrBuildMaterial("G4_Rn");
      break;

      case 87:

      theElement  = man->FindOrBuildElement("Fr");
      theMaterial = man->FindOrBuildMaterial("G4_Fr");
      break;

      case 88:

      theElement  = man->FindOrBuildElement("Ra");
      theMaterial = man->FindOrBuildMaterial("G4_Ra");
      break;

      case 89:

      theElement  = man->FindOrBuildElement("Ac");
      theMaterial = man->FindOrBuildMaterial("G4_Ac");
      break;

      case 90:

      theElement  = man->FindOrBuildElement("Th");
      theMaterial = man->FindOrBuildMaterial("G4_Th");
      break;

      case 91:

      theElement  = man->FindOrBuildElement("Pa");
      theMaterial = man->FindOrBuildMaterial("G4_Pa");
      break;

      case 92:

      theElement  = man->FindOrBuildElement("U");
      theMaterial = man->FindOrBuildMaterial("G4_U");
      break;
  
    }

    iz   = G4int(theElement->GetZ()); 
    N    = G4int(theElement->GetN()+0.5);
        N   -= iz;
    if( N < 0 ) N  = 0;
    G4bool boolChips = true; 




   //   G4cout << "Dumping element info:" << G4endl;
   //   theElement->DumpInfo();
   //   G4cout << "Dumping material info:" << G4endl;
   //   theMaterial->DumpInfo();

     //   sig = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,
     //                                                  theElement);

// geishaXsc = theCrossSectionDataStore.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // G4cout << "GHEISHA" <<" \t"<< sig/millibarn << " mb" << G4endl;


    // G4double momentum = theDynamicParticle->GetTotalMomentum();
    // if(i > 0) boolChips = false;
    // sig = qElastic->GetCrossSection(boolChips,momentum,iz,N,pPDG);

    // chipsRat = qRatio->GetRatios(momentum,pPDG,iz,N);
    // ratio = chipsRat.first*chipsRat.second; // !!!

    // G4cout << "Q-elastic          " <<" \t"<< sig/millibarn << " mb" << G4endl;


    // sig = hpwPrIn.GetCrossSection(theDynamicParticle, theElement, 273*kelvin);
    // hpwXsc = hpwNeIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // G4cout << kinEnergy/GeV <<" GeV, \t"<< sig/millibarn << " mb" << G4endl;

    barashInXsc  = barNucIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    barashTotXsc = barNucIn.GetTotalXsc();
    barashElXsc  = barNucIn.GetElasticXsc();

    // neutronhpXsc = nhpXscData.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);

    // barashInXsc  = barPiIn.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);
    // barashTotXsc = barPiIn.GetTotalXsc();
    // barashElXsc  = barPiIn.GetElasticXsc();

     ggTotXsc  = ggXsc.GetCrossSection(theDynamicParticle,theElement, 273*kelvin);

     ggIneXsc  = ggXsc.GetInelasticGlauberGribovXsc();

     ggElaXsc  = ggXsc.GetElasticGlauberGribovXsc();

     // ggElaXsc = ggXsc.GetHadronNucleaonXscPDG(theDynamicParticle, theElement);


     ratioTot = barashTotXsc/ggTotXsc;
     ratioIn  = barashInXsc/ggIneXsc;
     ratioEl  = barashElXsc/ggElaXsc;


    // G4cout << theProcess->GetProcessName() << " cross section for " << 
    //        theParticleDefinition->GetParticleName() <<
    //    " on " << theElement->GetName() << G4endl;
    // G4cout << kinEnergy/GeV <<" GeV, \t"<< sig/millibarn << " mb" << G4endl;
     //  G4cout << "Mean free path = " << mfp << " mm" << G4endl;

    G4cout <<i<<"    ratioTot = " << ratioTot << " ; ratioIn = "
       << ratioIn << " ; ratioEl = "<< ratioEl << G4endl;      
    writef <<i<<"    ratioTot = " << ratioTot << " ; ratioIn = "
       << ratioIn << " ; ratioEl = "<< ratioEl << G4endl;

    writeCorTot<<ratioTot<<", ";
    writeCorIn<<ratioIn<<", ";

    // delete theElement;
    // delete theMaterial;
      

  }
  writeCorTot<<G4endl;
  writeCorIn<<G4endl;

  //  G4cout<<"iz = "<<iz<<";  N = "<<N<<"; sum = "<<iz+N<<G4endl;                         
  
  /*  
  // Check Glauber-Gribov model
  G4cout<<"Check Glauber-Gribov model"<<G4endl;
  
                  
  kinEnergy = 0.1*GeV;
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
  // writef << kinEnergy/GeV <<"\t"<< ggIneXsc/millibarn <<"\t"<< ggProdXsc/millibarn << G4endl;
     writef << kinEnergy/GeV <<"\t"<< ratio <<G4endl;

     kinEnergy *= 1.145;
     delete theDynamicParticle;
  }
              
  */    

 
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

  G4cout<<"kin energy = "<<kinEnergy/GeV<<"in GeV"<<G4endl;
  G4cout << theProcess->GetProcessName() << " cross section for " << 
            theParticleDefinition->GetParticleName()  << G4endl;

  writef<<"kin energy = "<<kinEnergy/GeV<<"in GeV"<<G4endl;
  writef << theProcess->GetProcessName() << " cross section for " << 
            theParticleDefinition->GetParticleName()  << G4endl;


  return 1;
} // end of main
