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
// Unit test for coherent elastic models
//
//  18.05.07 V. Grichine
//
//

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
#include <complex>


#include "G4Element.hh"
#include "G4NistManager.hh"


#include "G4ParticleDefinition.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleMomentum.hh"


#include "G4MscRadiation.hh"


using namespace std;


G4double SimpleMF(G4double x)
{
  G4double order = 6.*x;
  return 1. - std::exp(-order);
}


G4double SimpleMPsi(G4double x)
{
  G4double order = 4.*x;
  return 1. - std::exp(-order);
}


G4double SimpleMG(G4double x)
{  
  return 3.*SimpleMPsi(x) - 2.*SimpleMF(x);
}



G4double MediumMF(G4double x)
{
  G4double order = 6.*x;
  order *= 1. + (3. - pi)*x;

  return 1. - std::exp(-order);
}


G4double MediumMPsi(G4double x)
{
  G4double order = 4.*x + 8.*x*x;
  return 1. - std::exp(-order);
}


G4double MediumMG(G4double x)
{  
  return 3.*MediumMPsi(x) - 2.*MediumMF(x);
}


G4double ComplexMF(G4double x)
{
  G4double order = 6.*x;
  order *= 1. + (3. - pi)*x;
  order -= x*x*x/(0.623+0.796*x+0.658*x*x);

  return 1. - std::exp(-order);
}


G4double ComplexMPsi(G4double x)
{
  G4double order = 4.*x + 8.*x*x/(1.+ 3.96*x + 4.97*x*x - 0.05*x*x*x + 7.5*x*x*x*x);
  return 1. - std::exp(-order);
}



G4double ComplexMG(G4double x)
{  
  return 3.*ComplexMPsi(x) - 2.*ComplexMF(x);
}



int main()
{

  G4int i, j, k, iMax;
  G4double x;
  G4double expXrad=0., g4Xrad;

  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  // iMax = 200;
  iMax = 0;
  // writef<<iMax<<G4endl;
  for( i = 1; i <= iMax; i++ )
  {
    x = 0.01*i;
    // G4cout<<x<<"\t"<<SimpleMF(x)<<"\t"<<MediumMF(x)<<"\t"<<ComplexMF(x)<<G4endl;
    // writef<<x<<"\t"<<SimpleMF(x)<<"\t"<<MediumMF(x)<<"\t"<<ComplexMF(x)<<G4endl;
  
    // G4cout<<x<<"\t"<<SimpleMPsi(x)<<"\t"<<MediumMPsi(x)<<"\t"<<ComplexMPsi(x)<<G4endl;
    // writef<<x<<"\t"<<SimpleMPsi(x)<<"\t"<<MediumMPsi(x)<<"\t"<<ComplexMPsi(x)<<G4endl;

    // G4cout<<x<<"\t"<<SimpleMG(x)<<"\t"<<MediumMG(x)<<"\t"<<ComplexMG(x)<<G4endl;
    // writef<<x<<"\t"<<SimpleMG(x)<<"\t"<<MediumMG(x)<<"\t"<<ComplexMG(x)<<G4endl;
  }

  G4Element*     theElement;
  G4Material*    theMaterial;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

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
  G4cout << "48 cadmium" << G4endl;
  G4cout << "74 tugnsten" << G4endl;
  G4cout << "77 iridium" << G4endl;
  G4cout << "82 lead" << G4endl;
  G4cout << "92 uranium" << G4endl;
  G4int choice;
  // G4cin >> choice;
  choice = 13;


  switch (choice)
  {
    case 1:

      theElement  = man->FindOrBuildElement("H");
      theMaterial = man->FindOrBuildMaterial("G4_H");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 2:

      theElement  = man->FindOrBuildElement("He");
      theMaterial = man->FindOrBuildMaterial("G4_He");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 4:

      theElement  = man->FindOrBuildElement("Be");
      theMaterial = man->FindOrBuildMaterial("G4_Be");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 6:

      theElement  = man->FindOrBuildElement("C");
      theMaterial = man->FindOrBuildMaterial("G4_C");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 19.6*cm;
      break;

    case 7:

      theElement  = man->FindOrBuildElement("N");
      theMaterial = man->FindOrBuildMaterial("G4_N");
      g4Xrad = theMaterial->GetRadlen();
      break;


    case 8:

      theElement  = man->FindOrBuildElement("O");
      theMaterial = man->FindOrBuildMaterial("G4_O");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 13:

      theElement  = man->FindOrBuildElement("Al");
      theMaterial = man->FindOrBuildMaterial("G4_Al");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 8.9*cm;
      break;

    case 14:

      theElement  = man->FindOrBuildElement("Si");
      theMaterial = man->FindOrBuildMaterial("G4_Si");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 18:

      theElement  = man->FindOrBuildElement("Ar");
      theMaterial = man->FindOrBuildMaterial("G4_Ar");
      break;

    case 26:

      theElement  = man->FindOrBuildElement("Fe");
      theMaterial = man->FindOrBuildMaterial("G4_Fe");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 1.76*cm;
      break;

    case 29:

      theElement  = man->FindOrBuildElement("Cu");
      theMaterial = man->FindOrBuildMaterial("G4_Cu");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 48:

      theElement  = man->FindOrBuildElement("Cd");
      theMaterial = man->FindOrBuildMaterial("G4_Cd");
      g4Xrad = theMaterial->GetRadlen();
      break;


    case 74:

      theElement  = man->FindOrBuildElement("W");
      theMaterial = man->FindOrBuildMaterial("G4_W");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 0.35*cm;
      break;

    case 77:

      theElement  = man->FindOrBuildElement("Ir");
      theMaterial = man->FindOrBuildMaterial("G4_Ir");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 82:

      theElement  = man->FindOrBuildElement("Pb");
      theMaterial = man->FindOrBuildMaterial("G4_Pb");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 0.56*cm;
      break;

    case 92:

      theElement  = man->FindOrBuildElement("U");
      theMaterial = man->FindOrBuildMaterial("G4_U");
      g4Xrad = theMaterial->GetRadlen();
      expXrad = 0.35*cm;
      break;
  }

// Particle definition

  G4cout << " 1 electron" << G4endl;
  G4cout << " 2 proton" << G4endl;
  G4cout << " 3 pion+" << G4endl;
  G4cout << " 4 pion-" << G4endl;
  G4cout << " 4 muon+" << G4endl;
  G4cout << " 5 muon-" << G4endl;

  //  G4cin >> choice;
  choice = 1;

  G4ParticleDefinition* theParticleDefinition;


  switch (choice)
  {
    case 1:

      theParticleDefinition = G4Electron::ElectronDefinition();  
      break;

    case 2:

      theParticleDefinition = G4Proton::ProtonDefinition();
      break;

    case 3:

      theParticleDefinition = G4PionPlus::PionPlusDefinition(); 
      break;

    case 4:

      theParticleDefinition = G4PionMinus::PionMinusDefinition();
      break;
 
    case 5:
     
      theParticleDefinition = G4MuonPlus::MuonPlusDefinition(); 
      break;

    case 6:
     
      theParticleDefinition = G4MuonMinus::MuonMinusDefinition();
      break;

  }

  G4double energyMscXR, kinEnergy;


  // kinEnergy = 8*GeV;   // 25.0*GeV;

  kinEnergy = 25.0*GeV;


  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              kinEnergy);

  G4double m1 = theParticleDefinition->GetPDGMass();
  G4double plab = theDynamicParticle->GetTotalMomentum();
  G4cout <<"lab momentum, plab = "<<plab/GeV<<" GeV"<<G4endl;
  G4double plabLowLimit = 20.0*MeV;

  G4int Z   = G4int(theElement->GetZ());
  G4int A    = G4int(theElement->GetN()+0.5);

  G4double step = 4.10*mm;

  step = expXrad;

  G4double m2 = man->GetAtomicMassAmu(Z)*GeV;
  // G4double m2 = man->GetAtomicMass( Z, A);
  G4cout <<" target mass, m2 = "<<m2/GeV<<" GeV"<<G4endl<<G4endl;
  G4cout <<"step = "<<step<<" mm; g4Xrad = "<<g4Xrad<<" mm; expXrad = "
         <<expXrad<<"  mm"<<G4endl<<G4endl;


  // G4MscRadiation* mscRad = new G4MscRadiation(theMaterial, step);
  G4MscRadiation* mscRad = new G4MscRadiation(theMaterial, step, theDynamicParticle);

  mscRad->SetVerboseLevel(0);

  G4double numberMscXR, numberMGYXR, numberE146XR, absorption, lincofXR, radLength; 

  G4double sRe, sIm, mcRe, mcIm, tmRe, tmIm, tmc, sf;

  G4complex ms, mc, tm;

  G4double s, phi, xi, G, psi, tmef;

  iMax = 50;

  writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {
    energyMscXR = std::exp(i*0.2)*0.1*MeV;
    lincofXR = mscRad->GetPlateLinearPhotoAbs(energyMscXR);
    absorption = (1. - std::exp(-lincofXR*step))/lincofXR;
 
    /*    
   
    numberMscXR  = mscRad->CalculateMscDiffdNdx(theDynamicParticle,energyMscXR);
    numberMGYXR  = mscRad->CalculateMscMigdalDiffdNdx(theDynamicParticle,energyMscXR);
    numberE146XR = mscRad->CalculateMscE146DiffdNdx(energyMscXR);

    radLength = mscRad->GetRadLength();

    numberMscXR  *= energyMscXR; // *expXrad;
    numberMGYXR  *= energyMscXR; // *expXrad;
    numberE146XR *= energyMscXR; // *expXrad;

    // numberMscXR *= expXrad;
    // numberMscXR *= expXrad;

    numberMscXR  *= absorption;
    numberMGYXR  *= absorption;
    numberE146XR *= absorption;

    G4cout <<"effStep = "<<absorption<<" mm; eXR  = "
           <<energyMscXR/MeV<<" MeV; MscXR =  "
           <<numberMscXR<<" "<<"; MGYXR =  "
           <<numberMGYXR<<" "<<"; E146XR =  "
           <<numberE146XR<<" "<<G4endl;
    writef <<energyMscXR/MeV<<"\t"<<numberMscXR<<"\t"<<numberMGYXR<<"\t"<<numberE146XR<<G4endl;

    
    

    mscRad->CalculateCorrectionTMGY(energyMscXR);
    tm = mscRad->GetCorrectionTMGY();
    ms = mscRad->CalculateMigdalS(energyMscXR);
    mc = mscRad->CalculateCorrectionMsc(energyMscXR);
    tmc = real(mc/tm);
    

    // sf = mscRad->SupressionFunction(kinEnergy, energyMscXR);

     mscRad->CalcLPMFunctions(kinEnergy, energyMscXR);

     s = mscRad->GetMigdalS();     
     phi = mscRad->GetMigdalPhi();     
     xi = mscRad->GetMigdalXi();     
     psi = mscRad->GetMigdalPsi();     
     G = mscRad->GetMigdalG(); 
     tmef = mscRad->GetTMEffect();    
         
    G4cout <<energyMscXR/MeV<<"\t"<<real(tm)<<"\t"<<imag(tm)
           <<"\t"<<real(ms)<<"\t"<<imag(ms)
           <<"\t"<<real(mc)<<"\t"<<imag(mc)<<"\t"<<sf<<G4endl; 
    writef<<energyMscXR/MeV<<"\t"<<real(tm)<<"\t"<<real(mc)<<"\t"<<tmc<<"\t"<<sf<<G4endl;
    
   

    // G4cout <<energyMscXR/MeV<<"\t"<<sf<<G4endl;
    // G4cout <<energyMscXR/MeV<<"\t"<<s<<G4endl;
    // writef <<energyMscXR/MeV<<"\t"<<s<<G4endl;
    // writef <<energyMscXR/MeV<<"\t"<<real(mc)<<G4endl;
    // writef <<energyMscXR/MeV<<"\t"<<s<<G4endl;


    G4cout <<energyMscXR/MeV<<"\t"<<phi<<"\t"<<xi<<"\t"<<psi<<"\t"<<G<<"\t"<<tmef<<G4endl;
    writef <<energyMscXR/MeV<<"\t"<<phi<<"\t"<<xi<<"\t"<<psi<<"\t"<<G<<"\t"<<tmef<<G4endl;

*/ 



  }


  G4double eTkin, lambda;
  iMax = 50;
  
  // writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {

    eTkin = std::exp(i*0.5)*0.001*MeV;

    G4DynamicParticle*  aDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              eTkin);

    G4MscRadiation* aMscRad = new G4MscRadiation(theMaterial, step, aDynamicParticle);

    lambda = aMscRad->CalculateMscE146Lambda();

    G4cout <<eTkin/GeV<<"\t"<<lambda<<G4endl;
    writef <<eTkin/GeV<<"\t"<<lambda<<G4endl;


    delete aMscRad;
    delete aDynamicParticle;
  }




  return 1;
} // end of main
