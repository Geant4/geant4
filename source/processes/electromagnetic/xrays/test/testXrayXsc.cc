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
// Unit test for X-ray xsc models
//
//  08.03.11 V. Grichine
//
//

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4Timer.hh"
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
#include "G4Gamma.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleMomentum.hh"

#include "G4VEmModel.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4HeatedKleinNishinaCompton.hh"
#include "G4PairProductionRelModel.hh"
#include "G4BetheHeitlerModel.hh"

// Photo electric effect

#include "G4PEEffectModel.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4Penelope08PhotoElectricModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

// Rayleigh (coherent elastic) scattering

#include "G4LivermoreRayleighModel.hh"
#include "G4PenelopeRayleighModel.hh"

///////////////////////////////////////////////

using namespace std;

const G4double cofA = 2.*pi2*Bohr_radius*Bohr_radius;

const G4double cofR = 8.*pi*classic_electr_radius*classic_electr_radius/3.;

//////////////////////////////////////////////////

G4double CalculateRayleighFF( G4double Z, G4double energy, G4double theta )
{
  G4double k    = energy/hbarc;
  G4double k2   = k*k;
  G4double z2p3 = std::pow(Z,2./3.);
  G4double a    = 1. + cofA*z2p3*k2*( 1. - std::cos(theta) );
  G4double a2   = a*a;
  G4double ff   = Z/a2;
  return ff;
}

//////////////////////////////////////////////////

G4double CalculateRayleighFF( G4double Z, G4double x)
{
  G4double z2p3 = std::pow(Z,2./3.);
  G4double x2   = x*x;
  G4double a    = 1. + cofA*z2p3*x2*2.;
  G4double a2   = a*a;
  G4double ff   = Z/a2;
  return ff;
}

////////////////////////////////////////////////////

G4double CalculateRayleighXsc( G4double Z, G4double energy )
{
  const G4double p0 =  1.78076e+00;  // 1.77457e+00;
  const G4double p1 = -6.0911e-02;    // -1.78171e-02;
  // const G4double p2 =  4.60444e-04;
  G4double    alpha = p0 + p1*std::log(Z); //  + p2*Z*Z;

  G4double k    = energy/hbarc;
  G4double k2   = std::pow(k, alpha); // k*k; 1.69 for Z=6, 1.5 for Z=82
  G4double z2p3 = std::pow( Z, 0.66667); // 2/3~0.66667
  G4double a    = cofA*z2p3*k2; 
  G4double b    = 1. + 2.*a;
  G4double b2   = b*b;
  G4double b3   = b*b2;
  G4double xsc  = cofR*Z*Z/b3;
           xsc *= a*a + (1. + a)*(1. + a);  
  return   xsc;   
}

G4double CalculateCosTheta(G4double c)
{
  G4double cosTheta, delta, cofA, signc = 1., a = c, power = 1./3.;
 
  if( c < 0.)
  {
    signc = -1.;
    a = -c;
  }
  delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*std::pow(delta, power);

  cosTheta = cofA - 1./cofA;
  return cosTheta; 
}



/////////////////////////////////////////////////////
//
//

int main()
{
  G4Timer timer;
  G4int i, j, k, iMax;
  G4double x;
  G4double expXrad=0., g4Xrad;

  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  /*
  iMax = 100;
  writef<<iMax<<G4endl;
  for( i = 0; i < iMax; i++ )
  {
    x = -4. + 8.*G4double(i)/G4double(iMax);
    g4Xrad = CalculateCosTheta(x);
    G4cout<<x<<"\t"<<g4Xrad<<G4endl;
    writef<<x<<"\t"<<g4Xrad<<G4endl;
  }
  */


  G4Element*     theElement;
  G4Material*    theMaterial;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  G4int choice;
  // G4cin >> choice;
  choice = 29;

  G4cout << " 1 hydrogen" << G4endl;
  G4cout << " 2 helium" << G4endl;
  G4cout << " 4 berillium" << G4endl;
  G4cout << " 6 carbon" << G4endl;
  G4cout << " 7 nitrogen" << G4endl;
  G4cout << " 8 oxigen" << G4endl;
  G4cout << "10 neon" << G4endl;
  G4cout << "13 aluminium" << G4endl;
  G4cout << "14 silicon" << G4endl;
  G4cout << "18 argon" << G4endl;
  G4cout << "26 iron" << G4endl;
  G4cout << "29 copper" << G4endl;
  G4cout << "48 cadmium" << G4endl;
  G4cout << "54 xenon" << G4endl;
  G4cout << "74 tugnsten" << G4endl;
  G4cout << "77 iridium" << G4endl;
  G4cout << "79 gold" << G4endl;
  G4cout << "82 lead" << G4endl;
  G4cout << "92 uranium" << G4endl;
  G4cout << "200 water" << G4endl;
  G4cout << "201 CO2" << G4endl;

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

    case 10:

      theElement  = man->FindOrBuildElement("Ne");
      theMaterial = man->FindOrBuildMaterial("G4_Ne");
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


    case 54:

      theElement  = man->FindOrBuildElement("Xe");
      theMaterial = man->FindOrBuildMaterial("G4_Xe");
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

    case 79:

      theElement  = man->FindOrBuildElement("Au");
      theMaterial = man->FindOrBuildMaterial("G4_Au");
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

    case 200:
     
      theMaterial = man->FindOrBuildMaterial("G4_WATER");
      g4Xrad = theMaterial->GetRadlen();
      break;

    case 201:

      theMaterial = man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
      // g4Xrad = theMaterial->GetRadlen();
      break;



  }

// Particle definition


  G4cout << " 0 gamma" << G4endl;
  G4cout << " 1 electron" << G4endl;
  G4cout << " 2 proton" << G4endl;
  G4cout << " 3 pion+" << G4endl;
  G4cout << " 4 pion-" << G4endl;
  G4cout << " 4 muon+" << G4endl;
  G4cout << " 5 muon-" << G4endl;

  //  G4cin >> choice;
  choice = 0;

  G4ParticleDefinition* theParticleDefinition;


  switch (choice)
  {
    case 0:

      theParticleDefinition = G4Gamma::GammaDefinition();  
      break;

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

  G4double energyXray, xsc=0., xsccomp, xscray, xscpair, xscphoto, xscref, delta;



  const G4DataVector cuts(1,0.01);


  G4VEmModel* penray = new G4PenelopeRayleighModel();
  G4VEmModel* livray = new G4LivermoreRayleighModel();
              livray->Initialise(theParticleDefinition,cuts);


  G4VEmModel* comp = new G4KleinNishinaCompton();
  // G4VEmModel* compheated = new G4HeatedKleinNishinaCompton();

  G4VEmModel* photo      = new G4PEEffectModel();
  G4VEmModel* fluphoto   = new G4PEEffectFluoModel();

  G4VEmModel* penphoto   = new G4PenelopePhotoElectricModel();
              penphoto->Initialise(theParticleDefinition,cuts);
  G4VEmModel* pen08photo = new G4Penelope08PhotoElectricModel();
              pen08photo->Initialise(theParticleDefinition,cuts);
  G4VEmModel* livphoto   = new G4LivermorePhotoElectricModel();
              livphoto->Initialise(theParticleDefinition,cuts);



  G4VEmModel* pair = new G4PairProductionRelModel(theParticleDefinition,"pp");
  // G4VEmModel* bhpair = new G4BetheHeitlerModel(theParticleDefinition,"bhpp");





  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(0.,0.,1.),
                                              energyXray);

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


  energyXray = 0.02*MeV;
  //xsc = pair->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );


  /*  

  timer.Start();
  
  iMax = 1000000; // 140; // 30; // 130;
  // writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {
    // energyXray = std::exp(i*0.08)*0.01*keV; // for iMax=140
    // energyXray = std::exp(i*0.11)*0.01*keV; // for iMax=130
    // energyXray = std::exp(i*0.41)*0.01*keV; // for iMax=30

    // xscray    = CalculateRayleighXsc(theElement->GetZ(),energyXray);

    // xscray = penray->CrossSectionPerVolume(theMaterial,theParticleDefinition,energyXray)
    //          /theMaterial->GetTotNbOfAtomsPerVolume();

    // xscray = livray->CrossSectionPerVolume(theMaterial,theParticleDefinition,energyXray)
    //      /theMaterial->GetTotNbOfAtomsPerVolume();



    // xscphoto = photo->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = fluphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = penphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    xscphoto = pen08photo->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = livphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );


    // xsccomp = comp->ComputeCrossSectionPerAtom(theParticleDefinition,energyXray,G4double(Z ) );
    // xscpair = pair->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscpair = bhpair->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );

    xsc = xscphoto; // xscray + xsccomp; // + xscphoto + xscpair;

    // G4cout<<energyXray/eV<<"\t\t"<<xsc/barn<<G4endl;
    // writef<<energyXray/eV<<"\t\t"<<xsc/barn<<G4endl;
  }
  timer.Stop();  
  G4cout<<"calculation time = "<<timer.GetUserElapsed()<<" s"<<G4endl;  
  */

  std::ifstream simRead;
  simRead.open("pe-cs-29.dat");

  simRead>>iMax;

  G4cout<<"iMax = "<<iMax<<G4endl;
  writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {
    simRead >> energyXray >> xscref;
    energyXray *= MeV;
    xscref     *= barn;

    // xscphoto = photo->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = fluphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = penphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    // xscphoto = pen08photo->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );
    xscphoto = livphoto->ComputeCrossSectionPerAtom(theParticleDefinition, energyXray, G4double(Z ) );

    delta = std::fabs(xscphoto-xscref)/xscref;

    
    G4cout<<energyXray/eV<<"\t\t"<<delta<<G4endl;
    writef<<energyXray/eV<<"\t\t"<<delta<<G4endl;
  }
  simRead.close();



  return 1;
} // end of main
