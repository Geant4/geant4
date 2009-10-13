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
#include "G4RandomTools.hh"
#include "G4Integrator.hh"
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



#include "G4WeMoSoftMscModel.hh"




using namespace std;


// Z-distribution for integral case I

G4double GetDeltaZDistribution(G4double v)
{
  G4double u, r;
  if ( v > 0.) u = std::exp(-1./v);
  else         u = 0.;

  if( v <= 3.)
  {
    r  = 2.*u*( 1 - 3.*std::pow(u, 8.) + 5.*std::pow(u, 24.) );
    r /= v*std::sqrt(pi*v);
  }
  else
  {
    r = 0.25*pi*std::exp(-pi*pi*v/16.);
  }
  return r;
}

// ln( G-function(x) )


G4double GammaLog(G4double xx)
{
  static G4double cof[6] = { 76.18009172947146,     -86.50532032941677,
                             24.01409824083091,      -1.231739572450155,
                              0.1208650973866179e-2, -0.5395239384953e-5  } ;
  register G4int j;
  G4double x = xx - 1.0 ;
  G4double tmp = x + 5.5 ;
  tmp -= (x + 0.5) * std::log(tmp) ;
  G4double ser = 1.000000000190015 ;

  for ( j = 0; j <= 5; j++ )
  {
    x += 1.0 ;
    ser += cof[j]/x ;
  }
  return -tmp + std::log(2.5066282746310005*ser) ;
}





////////////////////////////////////////////////////////////////////////////////
//
//

int main()
{

  G4int i, j, k, iMax, iStat;
  G4double x, xm, mm, a, B[200], sG[200], Bint[200], sGint[200],g, G[200], Gint[200], delta;
  G4double expXrad=0., g4Xrad, sumB=0., sumG=0.;

  for( i = 0; i < 200; i++)
  {
    B[i]     = 0.;
    Bint[i]  = 0.;
    G[i]     = 0.;
    Gint[i]  = 0.;
    sG[i]    = 0.;
    sGint[i] = 0.;
  }
  std::ofstream writef("angle.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  delta = 0.1;
  xm = 1.8;
  a  = 16./pi/pi;
  g  = GammaLog(a);
  g  = std::exp(g);
  G4cout<<"xm = "<<xm<<"; a = "<<a<<";   g(a) = "<<g<<G4endl;
  G4cout<<G4endl;



   iMax = 100;
  // iMax = 0;
  // writef<<iMax<<G4endl;

  for( i = 1; i <= iMax; i++ )
  {
    sG[i] = 0.;
    x     = delta*i;
    B[i]     = GetDeltaZDistribution(x);
    sumB += B[i]*delta;

    mm = a*x/xm;
    G[i] = std::pow(mm,a)*std::exp(-mm)/x/g;

    sumG += G[i]*delta;
    Bint[i] = 1. - sumB;
    Gint[i] = 1. - sumG;

    // G4cout<<x<<"\t"<<B[i]<<"\t"<<G[i]<<G4endl;
  }
  G4cout<<"sumB = "<<sumB<<"; sumG = "<<sumG<<G4endl;

  writef<<iMax<<G4endl;

  for( i = 1; i <= iMax; i++ )
  {
    x     = delta*i;
    // B[i] *= sumG/sumB; 
    // writef<<x<<"\t"<<B[i]<<"\t"<<G[i]<<G4endl;
  }

  iStat  = 100000;

  for( i = 1; i <= iStat; i++ )
  {
    x = CLHEP::RandGamma::shoot(a,a/xm);
    x += delta;

    for( j = 0; j < iMax; j++ )
    {
      if( x <= delta*j ) 
      {
        sG[j] += 1.;
        break;
      }
    } 
  }

  // normalise on 1.

  sumG = 0.;

  for( i = 0; i < iMax; i++ )
  {   
    sumG += delta*sG[i];
  }

  sumB = 0.;
 
  for( i = 1; i <= iMax; i++ )
  {   
    x     = delta*i;

    sG[i] /= sumG;
    // sG[i] /= delta*iStat;

    sumB    += sG[i]*delta;
    sGint[i] = 1. - sumB;


    // writef<<x<<"\t"<<B[i]<<"\t"<<G[i]<<"\t"<<sG[i]<<G4endl;
    writef<<x<<"\t"<<Bint[i]<<"\t"<<Gint[i]<<"\t"<<sGint[i]<<G4endl;
  } 


  // elements, materials using NIST

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




  return 1;
} // end of main
