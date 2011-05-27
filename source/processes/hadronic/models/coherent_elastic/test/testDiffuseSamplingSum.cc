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
// Unit test for all elastic models vs. experimental data
//
//  16.01.09 V. Grichine
//
//

#include "G4ios.hh"
#include "G4Timer.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

#include "G4DiffuseElastic.hh"
#include "G4HadronElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4QElasticCrossSection.hh"
#include "G4VQCrossSection.hh"
#include "G4ElasticHadrNucleusHE.hh"


#include "G4UHadronElasticProcess.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticle.hh"

#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include  "G4IonTable.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

#include "G4VCrossSectionDataSet.hh"


G4int verboseLevel = 0; 

G4double SampleThetaLab( G4DynamicParticle*  theDynamicParticle, G4double tmass, G4double t  




		// const G4HadProjectile* theDynamicParticle, G4double tmass, G4double A
              )
{

  const G4ParticleDefinition* theParticle = theDynamicParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  // G4double plab = theDynamicParticle->GetTotalMomentum();
  G4LorentzVector lv1 = theDynamicParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,tmass);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot    = p1.mag();

  G4double tmax    = 4.0*ptot*ptot;

  // NaN finder

  if(!(t < 0.0 || t >= 0.0)) 
  {
    if (verboseLevel > 0) 
    {
      G4cout << "G4DiffuseElastic:WARNING: "  
	// << " mom(GeV)= " << plab/GeV 
             << " S-wave will be sampled" 
	     << G4endl; 
    }
    t = G4UniformRand()*tmax; 
  }
  if(verboseLevel>1)
  {
    G4cout <<" t= " << t << " tmax= " << tmax 
	   << " ptot= " << ptot << G4endl;
  }
  // Sampling of angles in CM system

  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  if( cost >= 1.0 ) 
  {
    cost = 1.0;
    sint = 0.0;
  }
  else if( cost <= -1.0) 
  {
    cost = -1.0;
    sint =  0.0;
  }
  else  
  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }    
  if (verboseLevel>1) 
  {
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;

  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4ThreeVector np1 = nlv1.vect();

    // G4double theta = std::acos( np1.z()/np1.mag() );  // degree;

  G4double theta = np1.theta();

  return theta;
}


using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

int main()
{
  G4Timer timer;
  G4int i, k, iMax, iMod;
 
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

  choice = 2;


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
  G4cout << " 4 kaon+" << G4endl;
  G4cout << " 5 kaon0short" << G4endl;

  //  G4cin >> choice;
  choice = 1;

  G4ParticleDefinition* theParticleDefinition;

  // G4NucleonNuclearCrossSection* barash = new G4NucleonNuclearCrossSection();

  //   G4PiNuclearCrossSection* barash = new G4PiNuclearCrossSection();

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

  }

  // Initialisation of Geant4 hadron elastic models


  G4DiffuseElastic*       diffelastic = new G4DiffuseElastic(theParticleDefinition);
  G4ElasticHadrNucleusHE* hElastic    = new G4ElasticHadrNucleusHE();
  G4HadronElastic*        gElastic    = new G4HadronElastic(); // has Gheisha call



  // Physics data
  /*
  G4double momentum = 9.92*GeV;
  G4double pMass    = theParticleDefinition->GetPDGMass();
  G4double kinEnergy = std::sqrt(momentum*momentum + pMass*pMass) - pMass;
  */

  // G4double kinEnergy = 1.*GeV;
  G4double kinEnergy = 301.*GeV;

  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              kinEnergy);

  // G4HadProjectile* projectile = new G4HadProjectile(*theDynamicParticle); 

  G4double m1 = theParticleDefinition->GetPDGMass();
  G4double plab = theDynamicParticle->GetTotalMomentum();
  G4cout <<"lab momentum, plab = "<<plab/GeV<<" GeV"<<G4endl;
  G4double plabLowLimit = 20.0*MeV;

  G4int Z   = G4int(theElement->GetZ());
  G4int A    = G4int(theElement->GetN()+0.5);
  G4int N = A - Z;
  if(N < 0) N = 0;
  G4int projPDG = theParticleDefinition->GetPDGEncoding();


  /*
  G4ParticleDefinition * theTargetDef = 0;

  if      (Z == 1 && A == 1) theTargetDef = G4Proton::Proton();
  else if (Z == 1 && A == 2) theTargetDef = G4Deuteron::Deuteron();
  else if (Z == 1 && A == 3) theTargetDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theTargetDef = G4He3::He3();
  else if (Z == 2 && A == 4) theTargetDef = G4Alpha::Alpha();
  else                       theTargetDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);
 
  G4double m2 = theTargetDef->GetPDGMass();
  */


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
  G4double      tmax = 4.0*ptot*ptot;

  G4double      tDif, tGla, tChi, tGhe;



  // Choose generator
  G4bool swave = false;

  // S-wave for very low energy
  if(plab < plabLowLimit) swave = true;

  // Angle sampling

  G4int  numberOfExpPoints = 0;

  G4double g2 = GeV*GeV; 

  G4double sData, dData, thetaLab[200],  tData[200];
  G4double distrDif[200], distrGla[200],distrChi[200],distrGhe[200];

  std::ifstream simRead;
  // protons
  // simRead.open("pPbT1GeV.dat");
  // simRead.open("pSiT1GeV.dat");
  // simRead.open("pPb9p92GeVc.dat");
  // simRead.open("pCT1GeV.dat");
  // simRead.open("pOT1GeV.dat");
  // simRead.open("pHe4T45GeV.dat");
  simRead.open("pHe4T301GeV.dat");
  // pions
  // simRead.open("pipPb9p92GeVc.exp");
  // simRead.open("pimPb9p92GeVc.exp");
  // simRead.open("pipC9p92GeVc.exp");

  simRead>>numberOfExpPoints;
  G4cout<<"numberOfExpPoints = "<<numberOfExpPoints<<G4endl;

  for( i = 0; i < numberOfExpPoints; i++ )
  {
    tData[i] = 0.0;

    // simRead>>pData>>sData>>dData;
    // tData[i] = pData*GeV*GeV;

    // simRead>>tData[i]>>sData;
    simRead>>tData[i]>>sData>>dData;
    // simRead>>sData>>tData[i]>>dData;
    tData[i] *= g2;
    // tData[i] *= degree;
    // tData[i] *= mrad;
  }
  simRead.close();

  const G4int kAngle = 100; // numberOfExpPoints;

  // dData = ( tData[numberOfExpPoints-1] - tData[0] )/kAngle;

  // G4double thetaMin = 3.*degree;
  // G4double thetaMin = 3.*mrad;
  G4double thetaMin = 0.;

  dData = tData[numberOfExpPoints-1]/kAngle;

  for( k = 0; k < kAngle; k++) 
  {
    // thetaLab[k] = tData[0] + k*dData;

    thetaLab[k] = k*dData + thetaMin;
    distrDif[k] = 0;
    distrGla[k] = 0;
    distrChi[k] = 0;
    distrGhe[k] = 0;
  }
  std::ofstream writes("sigma.dat", std::ios::out ) ;
  writes.setf( std::ios::scientific, std::ios::floatfield );

  G4double thetaLabGla, thetaLabChi, thetaLabDif, thetaLabGhe;
  // G4double thetaCmsDif; // thetaCmsGla, thetaCmsChi,  thetaCmsGhe;
  
  G4VQCrossSection*       qCManager   = G4QElasticCrossSection::GetPointer();
  
  G4double cs = qCManager->GetCrossSection(false,plab,Z,N,projPDG);
  cs *= 1.;
  
  // Sampling loop up to iMax



  iMax = 1000000;   // numberOfExpPoints;
  // iMax = 1;   // numberOfExpPoints;
  iMod = iMax/10;
  writes << iMax  << G4endl;

  G4cout <<"Start Sampling ... "<<G4endl;

  timer.Start();

  // ptot = 0.1*GeV + 10*GeV*G4UniformRand();

  for( i = 0; i < iMax; i++)
  {
    // ptot = 0.1*GeV + 10*GeV*G4UniformRand();
           
    tDif = diffelastic->SampleTableT(theParticleDefinition, ptot, Z, A);
    // thetaLabDif =  SampleThetaLab( theDynamicParticle, m2, tDif );
    thetaLabDif = tDif;
    // tDif = diffelastic->SampleT(theParticleDefinition, ptot, A);

    // thetaCmsDif = diffelastic->SampleTableThetaCMS( theParticleDefinition, ptot, Z, A); 
    // thetaLabDif = diffelastic->ThetaCMStoThetaLab(theDynamicParticle, m2, thetaCmsDif);
    
    for( k = 0; k < kAngle; k++)
    {
      if( thetaLabDif <= thetaLab[k] )
      {
        distrDif[k] += 1;
        break;
      }
    }
       
    tGla = hElastic->SampleT(theParticleDefinition,plab,Z,A);
    // thetaLabGla = SampleThetaLab( theDynamicParticle, m2, tGla ); 
    thetaLabGla = tGla;
    for( k = 0; k < kAngle; k++)
    {
      if( thetaLabGla <= thetaLab[k] )
      {
        distrGla[k] += 1;
        break;
      }
    }
           
    tChi = qCManager->GetExchangeT(Z,N,projPDG);
    // thetaLabChi = SampleThetaLab( theDynamicParticle, m2, tChi );
    thetaLabChi = tChi;
    for( k = 0; k < kAngle; k++)
    {
      if( thetaLabChi <= thetaLab[k] )
      {
        distrChi[k] += 1;
        break;
      }
    }
        
    tGhe = g2*gElastic->SampleT(tmax/g2,m1,m2,G4double(A));
    // thetaLabGhe = SampleThetaLab( theDynamicParticle, m2, tGhe ); 
    thetaLabGhe = tGhe;
    for( k = 0; k < kAngle; k++)
    {
      if( thetaLabGhe <= thetaLab[k] )
      {
        distrGhe[k] += 1;
        break;
      }
    }
        
    if (i%iMod == 0) G4cout <<"done = "<<100.*G4double(i)/G4double(iMax)<<" %"<<G4endl;
  }
  timer.Stop();
  G4cout.precision(16);
  G4cout<<"statistics time = "<<timer.GetUserElapsed()<<" s"<<G4endl;
  G4cout.precision(6);

  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  // distrDif[0] = distrDif[1];



  writef <<kAngle<<G4endl;

  for( k = 0; k < kAngle; k++) 
  {
    // angleDistr[k] *= sig/iMax/(2*pi)/std::sin(k*thetaMax/kAngle);
    // angleDistr[k] *= steradian/millibarn;
    // G4cout <<k*thetaMax/kAngle/degree<<"\t"<<"\t"<<angleDistr[k]<<G4endl;
    // writef <<k*thetaMax/kAngle/degree<<"\t"<<angleDistr[k]<<G4endl;
    /*
    distrDif[k] /= 2*pi*std::sin(thetaLab[k]+0.001);
    distrGla[k] /= 2*pi*std::sin(thetaLab[k]+0.001);
    distrChi[k] /= 2*pi*std::sin(thetaLab[k]+0.001);
    distrGhe[k] /= 2*pi*std::sin(thetaLab[k]+0.001);
    */

    // G4cout <<thetaLab[k]/degree<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
    //     <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;
    //  writef <<thetaLab[k]/degree<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
    //       <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;

    // G4cout <<thetaLab[k]/mrad<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
    //          <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;
    //  writef <<thetaLab[k]/mrad<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
    //     <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;


    G4cout <<thetaLab[k]/g2<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
        <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;
     writef <<thetaLab[k]/g2<<"\t"<<distrDif[k]<<"\t"<<distrGla[k]
        <<"\t"<<distrChi[k]<<"\t"<<distrGhe[k]<<G4endl;
  
  }

  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << " elastic cross section for T = " <<kinEnergy/GeV<<" GeV   "<<
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;

  return 1;
} // end of main
