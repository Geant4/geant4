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
// Unit test for nucleus-nucleus diffuse elastic models
//
//  17.06.10 V. Grichine
//
//

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

#include "G4DiffuseElastic.hh"
#include "G4NuclNuclDiffuseElastic.hh"


#include "G4UHadronElasticProcess.hh"

#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

#include  "G4ParticleTable.hh"
#include  "G4IonTable.hh"
#include  "G4StateManager.hh"
#include  "G4DecayPhysics.hh"


using namespace std;


int main()
{

  G4int k; // i, j, iMax;
  // G4double x;

  G4DiffuseElastic* diffelastic = new G4DiffuseElastic();
  G4NuclNuclDiffuseElastic* nndiffelastic = new G4NuclNuclDiffuseElastic();

  /*
  std::ofstream writeb("bessel.dat", std::ios::out ) ;
  writeb.setf( std::ios::scientific, std::ios::floatfield );

  for( i = 0; i < 51; i++) // J0,1 test for the table from Hanbook of Special Functions p.390
  {
    x = 0.1*i;
    G4cout<<x<<"\t"<<diffelastic->BesselJzero(x)<<"\t"<<diffelastic->BesselJone(x)<<G4endl;
    writeb<<x<<"\t"<<diffelastic->BesselJzero(x)<<"\t"<<diffelastic->BesselJone(x)<<G4endl;
  }
  */

  // Element definition

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
  G4cout << "82 lead" << G4endl;
  G4cout << "92 uranium" << G4endl;
  G4int choice;
  // G4cin >> choice;

  choice = 8;


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

    case 48:

      theElement  = man->FindOrBuildElement("Cd");
      theMaterial = man->FindOrBuildMaterial("G4_Cd");
      break;


    case 74:

      theElement  = man->FindOrBuildElement("W");
      theMaterial = man->FindOrBuildMaterial("G4_W");
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

// Particle definition

  G4cout << " 1 proton" << G4endl;
  G4cout << " 2 neutron" << G4endl;
  G4cout << " 3 pion+" << G4endl;
  G4cout << " 4 pion-" << G4endl;
  G4cout << " 5 kaon+" << G4endl;
  G4cout << " 6 kaon0short" << G4endl;
  G4cout << " 7 generic ion" << G4endl;

  //  G4cin >> choice;
  choice = 7;

  G4ParticleDefinition* theParticleDefinition;

  //  G4NucleonNuclearCrossSection* barash = new G4NucleonNuclearCrossSection();

  // G4PiNuclearCrossSection* barash = new G4PiNuclearCrossSection();

  G4StateManager* g4State = G4StateManager::GetStateManager();
  if(! g4State->SetNewState(G4State_Init) ); 

  // C

  G4int Z1 = 6;
  G4int A1 = 12;

  G4DecayPhysics decays;
  decays.ConstructParticle();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  switch (choice)
  {
    case 1:

      theParticleDefinition = G4Proton::ProtonDefinition();
      // barash = new G4NucleonNuclearCrossSection();
      break;

    case 2:

      theParticleDefinition = G4Neutron::NeutronDefinition();
      // barash = new G4NucleonNuclearCrossSection();     
      break;

    case 3:

      theParticleDefinition = G4PionPlus::PionPlusDefinition();
      // barash = new G4PiNuclearCrossSection();
    
      break;

    case 4:

      theParticleDefinition = G4PionMinus::PionMinusDefinition();
      // barash = new G4PiNuclearCrossSection();
      
      break;
 
    case 5:

      theParticleDefinition = G4KaonPlus::KaonPlusDefinition();
     
      break;

    case 6:

      theParticleDefinition = G4KaonZeroShort::KaonZeroShortDefinition();
     
      break;

    case 7:

      theParticleDefinition = partTable->FindIon(Z1, A1, 0, Z1);
     
      break;

  }

  // G4double momentum = 9.92*GeV;

  // G4double pMass    = theParticleDefinition->GetPDGMass();

  // G4double thetaMax  = 15.*degree; 

  // G4double kinEnergy = std::sqrt(momentum*momentum + pMass*pMass) - pMass;

  G4double kinEnergy = 168.*MeV;

  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              kinEnergy);
  // G4double m1 = theParticleDefinition->GetPDGMass();
  G4double plab = theDynamicParticle->GetTotalMomentum();
  G4cout <<"lab momentum, plab = "<<plab/GeV<<" GeV"<<G4endl;
  G4double plabLowLimit = 20.0*MeV;

  G4int Z   = G4int(theElement->GetZ());
  G4int A    = G4int(theElement->GetN()+0.5);



  G4double m2 = man->GetAtomicMassAmu(Z)*GeV;
  // G4double m2 = man->GetAtomicMass( Z, A);
  G4cout <<" target mass, m2 = "<<m2/GeV<<" GeV"<<G4endl;

  G4LorentzVector lv1 = theDynamicParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,m2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1   = lv1.vect();
  G4double      ptot = p1.mag();
  G4cout <<"cms momentum, ptot = "<<ptot/GeV<<" GeV"<<G4endl;

  // Choose generator

  G4bool swave = false;

  // S-wave for very low energy

  if(plab < plabLowLimit) swave = true;

  // Angle sampling


  G4double dData, thetaLab[200];
  G4double distrDif[200], distrXsc[200];

  G4double thetaCMS;

  G4double rad, rad2;
  
  const G4int kAngle = 100; // number of angles
  G4double thetaMin = 15.*degree;
  dData = 30.*degree/kAngle;

  nndiffelastic->SetCofLambda(1.01);  // 1.01, 
  nndiffelastic->SetCofDelta(0.095);   // 0.092
  nndiffelastic->SetCofAlpha(0.095);   // 0.092
  nndiffelastic->SetCofPhase(1.03);   
  nndiffelastic->SetCofFar(1.0);    
  // nndiffelastic->SetEtaRatio(1.40);   // 1


  nndiffelastic->InitParameters(theParticleDefinition, ptot, Z, A);

  // nndiffelastic->InitParametersGla(theDynamicParticle, ptot, Z, A);
  // nndiffelastic->SetMaxL(10);   // 10
  // nndiffelastic->SetEtaRatio(7.50);   // 1

  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  writef<<kAngle<<G4endl;

  for( k = 0; k < kAngle; k++) 
  {
    thetaLab[k] = k*dData + thetaMin;
    thetaCMS    = thetaLab[k];
    
    distrDif[k] = diffelastic->GetDiffuseElasticSumXsc(theParticleDefinition, thetaCMS, ptot, A, Z); 

    // distrXsc[k] = nndiffelastic->AmplitudeMod2(thetaCMS);
    // distrXsc[k] = nndiffelastic->AmplitudeSimMod2(thetaCMS);
    distrXsc[k] = nndiffelastic->CoulombAmplitudeMod2(thetaCMS)*nndiffelastic->GetElCoulRatioSim(thetaCMS);
    // distrXsc[k] = nndiffelastic->AmplitudeGlaMod2(thetaCMS);
    // distrXsc[k] = nndiffelastic->AmplitudeGGMod2(thetaCMS);

    rad = diffelastic->GetNuclearRadius();
    rad2 = rad*rad;

    distrDif[k] /= rad2;
    distrXsc[k] /= rad2;

    G4cout <<thetaLab[k]/degree<<"\t"<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;
    writef <<thetaLab[k]/degree<<"\t"<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;
  
  }




  G4cout <<G4endl<< " elastic cross section for " <<
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;

  return 1;
} // end of main
