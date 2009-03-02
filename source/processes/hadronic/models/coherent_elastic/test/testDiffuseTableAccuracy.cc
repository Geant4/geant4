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
#include "G4Integrator.hh"

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

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

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

  // Physics data

  // G4double momentum = 9.92*GeV;
  // G4double pMass    = theParticleDefinition->GetPDGMass();
  // G4double kinEnergy = std::sqrt(momentum*momentum + pMass*pMass) - pMass;

  // G4double kinEnergy = 1.*GeV;
  G4double kinEnergy = 45.*GeV;

  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              kinEnergy);

  // G4HadProjectile* projectile = new G4HadProjectile(*theDynamicParticle); 

  // G4double m1 = theParticleDefinition->GetPDGMass();
  G4double plab = theDynamicParticle->GetTotalMomentum();
  G4cout <<"lab momentum, plab = "<<plab/GeV<<" GeV"<<G4endl;


  G4int Z   = G4int(theElement->GetZ());
  G4int A    = G4int(theElement->GetN()+0.5);
  G4int N = A - Z;
  if(N < 0) N = 0;
 


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
  // G4double      tmax = 4.0*ptot*ptot;


  // Angle sampling

  G4int  numberOfExpPoints = 0;

  G4double sData, dData, thetaLab[200],  tData[200];
  G4double distrDif[200],distrXsc[200];

  G4double thetaLabDif, thetaCmsDif, thetaCmsDif2, tDif;

  G4double g2 = GeV*GeV;

  // read exp data and set angle simulation limit

  std::ifstream simRead;
  // protons
  // simRead.open("pPbT1GeV.dat");
  // simRead.open("pSiT1GeV.dat");
  // simRead.open("pPb9p92GeVc.dat");
  // simRead.open("pCT1GeV.dat");
  // simRead.open("pOT1GeV.dat");
  simRead.open("pHe4T45GeV.dat");
  // simRead.open("pHe4T301GeV.dat");
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
  G4double thetaMin = 0.;

  dData = tData[numberOfExpPoints-1]/kAngle;

  for( k = 0; k < kAngle; k++) 
  {
    // thetaLab[k] = tData[0] + k*dData;

    thetaLab[k] = k*dData + thetaMin;
    distrDif[k] = 0;

    thetaLabDif = thetaLab[k];

    // thetaCmsDif = diffelastic->ThetaLabToThetaCMS(theDynamicParticle,m2,thetaLabDif);
    thetaCmsDif = std::sqrt(thetaLabDif/ptot/ptot);

    distrXsc[k] = diffelastic->GetDiffuseElasticSumXsc( theParticleDefinition, thetaCmsDif, ptot, A, Z);
    distrXsc[k] /= diffelastic->GetNuclearRadius()*diffelastic->GetNuclearRadius();

  }



  std::ofstream writes("testTable.dat", std::ios::out ) ;
  writes.setf( std::ios::scientific, std::ios::floatfield );


  // Sampling loop up to iMax



  // iMax = 1000000;   // numberOfExpPoints;
  iMax = 1;   // numberOfExpPoints;
  iMod = iMax/10;
 
  /*
  G4cout <<"Start table testing ... "<<G4endl;

  timer.Start();

  for( i = 0; i < iMax; i++)
  {
    
    diffelastic->TestAngleTable(theParticleDefinition, ptot, Z, A);

    
    // if (i%iMod == 0) G4cout <<"done = "<<100.*G4double(i)/G4double(iMax)<<" %"<<G4endl;
  }
  timer.Stop();
  G4cout.precision(16);
  G4cout<<"statistics time = "<<timer.GetUserElapsed()<<" s"<<G4endl;
  G4cout.precision(6);

  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << " elastic cross section for " <<
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;

  // Cross-section test

  G4double sigma[200], sumsigma[200],coulsigma[200], samplesigma[200];
  G4double dif, mean;

  // writes<<numberOfExpPoints<< G4endl;

  for( i = 0; i < numberOfExpPoints; i++)
  {

    thetaLabDif = tData[i];

    thetaCmsDif = diffelastic->ThetaLabToThetaCMS(theDynamicParticle,m2,thetaLabDif);

    thetaCmsDif2 = thetaCmsDif*thetaCmsDif;

    sigma[i] = diffelastic->GetDiffuseElasticXsc( theParticleDefinition, thetaCmsDif, ptot, A);
    sumsigma[i] = diffelastic->GetDiffuseElasticSumXsc( theParticleDefinition, thetaCmsDif, ptot, A, Z);
    samplesigma[i] = diffelastic->GetDiffElasticSumProbA( thetaCmsDif2);
    // samplesigma[i] = diffelastic->GetDiffElasticSumProb( thetaCmsDif);
    samplesigma[i] *= diffelastic->GetNuclearRadius()*diffelastic->GetNuclearRadius();

    dif = std::abs(sumsigma[i]-samplesigma[i]);
    mean = 0.5*(sumsigma[i]+samplesigma[i]);

    coulsigma[i] = diffelastic->GetCoulombElasticXsc( theParticleDefinition, thetaCmsDif, ptot, Z);

     G4cout << thetaLabDif/degree << "\t" << "\t" << sigma[i]/millibarn << "\t\t" 
            << sumsigma[i]/millibarn<< "\t\t" << samplesigma[i]/millibarn << "\t\t" <<dif/mean<< G4endl;
     // writes << thetaLabDif/degree << "\t" << sigma[i]/millibarn << "\t" 
     //        << sumsigma[i]/millibarn<< "\t" << samplesigma[i]/millibarn <<  "\t" <<dif/mean<<G4endl;

  }
  */

  iMax = 1000000;   // numberOfExpPoints;
  // iMax = 1;   // numberOfExpPoints;
  iMod = iMax/10;
  // writes << iMax  << G4endl;

  G4cout <<"Start Sampling ... "<<G4endl;

  timer.Start();

  for( i = 0; i < iMax; i++)
  {
    
    tDif = diffelastic->SampleTableT(theParticleDefinition, ptot, Z, A);
    // thetaLabDif =  SampleThetaLab( theDynamicParticle, m2, tDif );

    thetaLabDif =  tDif;

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
    if (i%iMod == 0) G4cout <<"done = "<<100.*G4double(i)/G4double(iMax)<<" %"<<G4endl;
  }
  timer.Stop();
  G4cout.precision(16);
  G4cout<<"statistics time = "<<timer.GetUserElapsed()<<" s"<<G4endl;
  G4cout.precision(6);

  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );




  distrDif[0] = distrDif[1];
  distrXsc[0] = distrXsc[1];

  writef <<kAngle<<G4endl;

  for( k = 0; k < kAngle; k++) 
  {
    // angleDistr[k] *= sig/iMax/(2*pi)/std::sin(k*thetaMax/kAngle);
    // angleDistr[k] *= steradian/millibarn;
    // G4cout <<k*thetaMax/kAngle/degree<<"\t"<<"\t"<<angleDistr[k]<<G4endl;
    // writef <<k*thetaMax/kAngle/degree<<"\t"<<angleDistr[k]<<G4endl;

    // distrDif[k] /= 2*pi*std::sin(thetaLab[k]+0.001);

    // G4cout <<thetaLab[k]/degree<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;
    // writef <<thetaLab[k]/degree<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;

    // G4cout <<thetaLab[k]/mrad<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;
    // writef <<thetaLab[k]/mrad<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;

    G4cout <<thetaLab[k]/g2<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;
    writef <<thetaLab[k]/g2<<"\t"<<distrDif[k]<<"\t"<<distrXsc[k]<<G4endl;

  
  
  }


  for( k = 0; k < kAngle; k++) 
  {

    thetaLabDif = thetaLab[k];

    thetaCmsDif = diffelastic->ThetaLabToThetaCMS(theDynamicParticle,m2,thetaLabDif);

    tDif = 2*ptot*ptot*( 1. - std::cos(thetaCmsDif) );

    thetaLabDif =  SampleThetaLab( theDynamicParticle, m2, tDif );

    G4cout << thetaLab[k]/degree << "\t" << thetaCmsDif/degree<< "\t" << thetaLabDif/degree << G4endl;

  }



  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << " elastic cross section for " <<
            theParticleDefinition->GetParticleName() <<
           " on " << theElement->GetName() << G4endl;
  G4cout <<"with atomic weight = "<<theElement->GetN() << G4endl;

  return 1;
} // end of main
