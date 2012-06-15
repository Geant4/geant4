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
// $Id: G4GGNNCrossSectionDataTest.cc,v 1.5 2010-06-09 08:29:47 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test of G4GGNuclNuclCrossSection class
//
// History:
//
// 24.11.08 V.Grichine: NIST elements/materials, write in file
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
#include "G4ParticleTable.hh"
#include "G4DecayPhysics.hh"
#include "G4DynamicParticle.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"

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
#include "G4GGNuclNuclCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"



#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonsKoxCrossSection.hh"
#include "G4IonsSihverCrossSection.hh"

#include "G4VQCrossSection.hh"
#include "G4QElasticCrossSection.hh"
#include "G4QuasiFreeRatios.hh"



#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif


int main()
{
  // Element definition

  G4cout << " 1 hydrogen" << G4endl;
  G4cout << " 2 helium" << G4endl;
  G4cout << " 4 berillium" << G4endl;
  G4cout << " 6 carbon" << G4endl;
  G4cout << " 7 nitrogen" << G4endl;
  G4cout << " 8 oxigen" << G4endl;
  G4cout << "13 aluminium" << G4endl;
  G4cout << "14 silicon" << G4endl;
  G4cout << "18 argon" << G4endl;
  G4cout << "26 iron" << G4endl;
  G4cout << "29 copper" << G4endl;
  G4cout << "39 ytrium" << G4endl;
  G4cout << "48 cadmium" << G4endl;
  G4cout << "74 tugnsten" << G4endl;
  G4cout << "77 iridium" << G4endl;
  G4cout << "82 lead" << G4endl;
  G4cout << "92 uranium" << G4endl;

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

      theElement  = man->FindOrBuildElement("H");
      theMaterial = man->FindOrBuildMaterial("G4_H");
      break;

    case 2:

      theElement  = man->FindOrBuildElement("He");
      theMaterial = man->FindOrBuildMaterial("G4_He");
      break;

    case 4:

      theElement  = man->FindOrBuildElement("Be");
      theMaterial = man->FindOrBuildMaterial("G4_Be");
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

    case 13:

      theElement  = man->FindOrBuildElement("Al");
      theMaterial = man->FindOrBuildMaterial("G4_Al");
      break;

    case 14:

      theElement  = man->FindOrBuildElement("Si");
      theMaterial = man->FindOrBuildMaterial("G4_Si");
      break;

    case 18:

      theElement  = man->FindOrBuildElement("Ar");
      theMaterial = man->FindOrBuildMaterial("G4_Ar");
      break;

    case 26:

      theElement  = man->FindOrBuildElement("Fe");
      theMaterial = man->FindOrBuildMaterial("G4_Fe");
      break;

    case 29:

      theElement  = man->FindOrBuildElement("Cu");
      theMaterial = man->FindOrBuildMaterial("G4_Cu");
      break;

    case 39:

      theElement  = man->FindOrBuildElement("Y");
      theMaterial = man->FindOrBuildMaterial("G4_Y");
      break;

    case 48:

      theElement  = man->FindOrBuildElement("Cd");
      theMaterial = man->FindOrBuildMaterial("G4_Cd");
      break;


    case 74:

      theElement  = man->FindOrBuildElement("W");
      theMaterial = man->FindOrBuildMaterial("G4_W");
      break;

    case 77:

      theElement  = man->FindOrBuildElement("Ir");
      theMaterial = man->FindOrBuildMaterial("G4_Ir");
      break;

    case 82:

      theElement  = man->FindOrBuildElement("Pb");
      theMaterial = man->FindOrBuildMaterial("G4_Pb");
      break;

    case 92:

      theElement  = man->FindOrBuildElement("U");
      theMaterial = man->FindOrBuildMaterial("G4_U");
      break;
  }


  G4StateManager* g4State = G4StateManager::GetStateManager();
  if (! g4State->SetNewState(G4State_Init) );

  // C
  G4int Z = 6;
  G4int A = 12;
  // Ne
  // G4int Z = 10;
  // G4int A = 20;

  G4DecayPhysics decays;
  decays.ConstructParticle();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4ParticleDefinition* theParticleDefinition = 
                        partTable->FindIon(Z, A, 0, Z);


  G4int i, iMax = 70;
  G4double kinEnergy;
 

 // G4cout << "Kinetic energy in GeV: "<<G4endl;
  // G4cin >> kinEnergy;

// Make a dynamic particle too

  // G4DynamicParticle* aDynamicParticle = 0;


  G4GGNuclNuclCrossSection ggXsc;

  G4IonsKoxCrossSection    koxXsc;
  G4IonsShenCrossSection   shenXsc;
  G4IonsSihverCrossSection sihverXsc;
  G4TripathiCrossSection   triXsc;

  G4double ggTotXsc, ggElaXsc, ggIneXsc, ggProdXsc;
  // G4double ggDifXsc, ratio;
  G4double koxInXsc, shenInXsc, sihverInXsc, triInXsc;


  iMax = 100;

  // Check Glauber-Gribov model
  G4cout<<"Check Glauber-Gribov model"<<G4endl;

  std::ofstream writegg("ggtxs.dat", std::ios::out ) ;
  writegg.setf( std::ios::scientific, std::ios::floatfield );

  
                  
  kinEnergy = 1.*MeV;
  iMax = 100;
  writegg <<iMax<< G4endl; 
 

  for(i = 0; i < iMax; i++)
  {
   
    G4DynamicParticle* aDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);

     ggIneXsc = ggXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);

     ggElaXsc = ggXsc.GetElasticGlauberGribovXsc();
     // ggElaXsc = ggXsc.GetHadronNucleaonXscPDG(theDynamicParticle, theElement);

     ggTotXsc  = ggXsc.GetTotalGlauberGribovXsc();
     ggProdXsc = ggXsc.GetProductionGlauberGribovXsc();
     // ggDifXsc  = ggXsc.GetDiffractionGlauberGribovXsc();

     // if(ggIneXsc != 0.) ratio = 1 - ggProdXsc/ggIneXsc;
     // if(ratio < 0.) ratio = 0.;
     // if(ggIneXsc != 0.) ratio = ggDifXsc/ggIneXsc;

     // ratio = ggXsc.GetRatioSD(theDynamicParticle, theElement->GetN(), theElement->GetZ() );
     // ratio = ggXsc.GetRatioQE(theDynamicParticle, theElement->GetN(), theElement->GetZ() );

     G4cout << kinEnergy/A <<" MeV/n,   tot = "<<  ggTotXsc/millibarn 
     << " mb;   el = "
     <<  ggElaXsc/millibarn << " mb;   in = " <<  ggIneXsc/millibarn<< " mb;   prod = " 
     <<  ggProdXsc/millibarn << G4endl;

     writegg << kinEnergy/A <<"  "<<  ggTotXsc/millibarn // <<  "  "<<  ggElaXsc/millibarn 
             <<"  " <<  ggIneXsc/millibarn <<"  " <<  ggProdXsc/millibarn
        <<"  " <<  ggElaXsc/millibarn 
             << G4endl;
    
     // G4cout << kinEnergy/GeV <<" GeV, in Xsc = "<<  ggIneXsc/millibarn << " mb; prod Xsc = "
     //       <<  ggProdXsc/millibarn << " mb; ratio = 1 - prod/in = " <<  ratio << G4endl;
    
     // G4cout << kinEnergy/GeV <<" GeV, in-Xsc = "<<  ggIneXsc/millibarn << " mb; sdif-Xsc = "
     //        <<  ggDifXsc/millibarn << " mb; ratio = dif/in = " <<  ratio << G4endl;
    

  // writegg << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
     // writegg << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
  // writegg << kinEnergy/GeV <<"\t"<< ggIneXsc/millibarn <<"\t"<< ggProdXsc/millibarn << G4endl;
     //     writegg << kinEnergy/GeV <<"\t"<< ratio <<G4endl;

     kinEnergy *= 1.145;
     delete aDynamicParticle;
  }
              
      


  // Check Nucleus-Nucleus inelastic Xsc models
  G4cout<<"Check Nucleus-Nucleus inelastic Xsc models"<<G4endl;
 
  std::ofstream writenn("nntxs.dat", std::ios::out ) ;
  writenn.setf( std::ios::scientific, std::ios::floatfield );

  
                  
  kinEnergy = 1.*MeV;
  iMax = 120;
  writenn <<iMax<< G4endl; 
 

  for(i = 0; i < iMax; i++)
  {
   
    G4DynamicParticle* aDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);

     ggIneXsc    = ggXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);
     koxInXsc    = koxXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);
     shenInXsc   = shenXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);
     sihverInXsc = sihverXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);
     triInXsc    = triXsc.GetCrossSection(aDynamicParticle,
                                 theElement, 273*kelvin);

     G4cout << kinEnergy/A <<" MeV/n,   gg = "<<  ggIneXsc/millibarn 
     << " mb;   kox = "
     <<  koxInXsc/millibarn << " mb;   shen = " <<  shenInXsc/millibarn<< " mb;   sihver = " 
     <<  sihverInXsc/millibarn << G4endl;

     writenn << kinEnergy/A <<"  "<<  ggIneXsc/millibarn <<  "  "<<  triInXsc/millibarn 
             <<"  " <<  shenInXsc/millibarn <<"  " <<  sihverInXsc/millibarn << G4endl;
    
     // G4cout << kinEnergy/GeV <<" GeV, in Xsc = "<<  ggIneXsc/millibarn << " mb; prod Xsc = "
     //       <<  ggProdXsc/millibarn << " mb; ratio = 1 - prod/in = " <<  ratio << G4endl;
    
     // G4cout << kinEnergy/GeV <<" GeV, in-Xsc = "<<  ggIneXsc/millibarn << " mb; sdif-Xsc = "
     //        <<  ggDifXsc/millibarn << " mb; ratio = dif/in = " <<  ratio << G4endl;
    

  // writenn << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
     // writenn << kinEnergy/GeV <<"\t"<< ggTotXsc/millibarn <<"\t"<< ggIneXsc/millibarn << G4endl;
  // writenn << kinEnergy/GeV <<"\t"<< ggIneXsc/millibarn <<"\t"<< ggProdXsc/millibarn << G4endl;
     //     writenn << kinEnergy/GeV <<"\t"<< ratio <<G4endl;

     kinEnergy *= 1.145;
     delete aDynamicParticle;
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
  G4cout << " cross section for " << 
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;


  return 1;
} // end of main
