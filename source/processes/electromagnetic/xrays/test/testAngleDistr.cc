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
// Unit test for relativistic EM angle distribution  models
//
//  02.05.11 V. Grichine
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
#include "G4ModifiedTsai.hh"


///////////////////////////////////////////////

using namespace std;

const G4double cofA = 2.*pi2*Bohr_radius*Bohr_radius;

const G4double cofR = 8.*pi*classic_electr_radius*classic_electr_radius/3.;


///////////////////////////////////////////////////////////////
//
// c = 8 - 4*r, where r is the uniform random (0,1)

G4double CalculateCosTheta(G4double c)
{
  G4double cosTheta, delta, cofA, signc = 1., a = c, power = 1./3.;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*std::pow(delta, power);

  cosTheta = cofA - 1./cofA;

  return cosTheta; 
}

///////////////////////////////////////////////////////////////
//
// 1 + cos^2(theta) random distribution in the projectile rest frame, 3rd equation algorithm

G4double RandCosThetaDipole()
{
  G4double c, cosTheta, delta, cofA, signc = 1., a, power = 1./3.;

  c = 4. - 8.*G4UniformRand();
  a = c;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*std::pow(delta, power);

  cosTheta = cofA - 1./cofA;

  return cosTheta; 
}


///////////////////////////////////////////////////////////////
//
// 1 + cos^2(theta) random distribution in the projectile rest frame, 3rd equation algorithm

G4double RandCosThetaDipPen()
{
  G4double x, cosTheta, signX, modX, power = 1./3.;

  if( G4UniformRand() > 0.25) 
  {
    cosTheta = 2.*G4UniformRand()-1.;
  }
  else
  {
    x     = 2.*G4UniformRand()-1.;

    if ( x < 0. ) 
    {
      modX = -x;
      signX = -1.;
    }
    else
    {
      modX = -x;
      signX = -1.;
    }
    cosTheta = signX*std::pow(modX,power);
  }
  return cosTheta; 
}

//////////////////////////////////////////////////////
//
// Bust the cosTheta to the laboratory frame

G4double BustCosTheta(G4double cosTheta, G4double beta)
{
  G4double result;
  result  = cosTheta + beta;
  result /= 1. + cosTheta*beta;
  return result;
}

////////////////////////////////////////////////////////////////////////////////////////
//
// vmg: clean-up of G4ModifiedTsai::PolarAngle method

G4double ModTsaiCosTheta( const G4double initial_energy ) 
{
  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution ( Rev Mod Phys 49,421(1977) )

  G4double gamma = 1. + initial_energy/electron_mass_c2;   
  G4double uMax  = gamma*pi;

  const G4double a1     = 0.625;
  const G4double a2     = 1.875;
  const G4double border = 0.25;
  G4double u, theta;

  do
  {
    u = - std::log( G4UniformRand()*G4UniformRand() );

    if ( border > G4UniformRand() ) u /= a1; 
    else                            u /= a2; 
    
  } while( u > uMax );

  theta = u/gamma;

  return std::cos(theta);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


int main()
{
  G4Timer timer;
  G4int i, iMax, ix;
  G4double cosTheta;
  G4double g4Xrad, expXrad, eTkin;

  std::ofstream writef("emangle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );


  G4Element*     theElement;
  G4Material*    theMaterial;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  G4int choice;
  // G4cin >> choice;
  choice = 82;
  /*  
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
  */
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
  G4double dipMy[100], dipPen[100], tsai[100];

  for( i = 0; i < 100; i++ )
  {
    dipMy[i]  = 0.;
    dipPen[i] = 0.;
    tsai[i]   = 0.;
  }
  // G4ModifiedTsai* modTsai = new G4ModifiedTsai();

  eTkin = 1.*MeV;

  G4double gamma = 1. + eTkin/electron_mass_c2;

  G4double beta = std::sqrt(1. - 1./gamma/gamma);

  iMax = 1000000;

  timer.Start();

  for( i = 0; i < iMax; i++ )
  {
    cosTheta = RandCosThetaDipole();
    cosTheta = BustCosTheta(cosTheta,beta);
    ix = G4int( 50.*(cosTheta + 1. ) + 0.5 );
    if( ix >=0 && ix < 200 ) dipMy[ix] += 1.; 


    cosTheta = RandCosThetaDipPen();
    cosTheta = BustCosTheta(cosTheta,beta);
    ix = G4int( 50.*(cosTheta + 1. ) + 0.5 );
    if( ix >=0 && ix < 200 ) dipPen[ix] += 1.; 

    cosTheta = ModTsaiCosTheta(eTkin); // already in lab frame
    ix = G4int( 50.*(cosTheta + 1. ) + 0.5 );
    if( ix >=0 && ix < 200 ) tsai[ix] += 1.; 
  }
  timer.Stop();  
  G4cout<<"calculation time = "<<timer.GetUserElapsed()<<" s"<<G4endl; 

  iMax = 100;

  writef << iMax << G4endl;

  for( i = 0; i < iMax; i++ )
  {
    cosTheta = -1. + 0.005 + i/50.;

    G4cout<<cosTheta<<"\t"<<dipMy[i]<<"\t"<<dipPen[i]<<"\t"<<tsai[i]<<G4endl;

    writef <<cosTheta<<"\t"<<dipMy[i]<<"\t"<<dipPen[i]<<"\t"<<tsai[i]<<G4endl;
  }

  return 1;
} // end of main
