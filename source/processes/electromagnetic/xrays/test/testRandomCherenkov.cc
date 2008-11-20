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
// Unit test for Cherenkov models in media with random fluctuations
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


G4double epsilon = 1.2;   //  1.05;

G4double cof = (fine_structure_const/hbarc)*1.*eV; // *50.*0.9*0.3;

// cof *= 50.*0.9*0.3;  // T=0.9, Q=0.3, 1 eV of energy range

using namespace std;


///////////////////////////////////////////////////////////
//
// CR spectrum in transparent medium with fluctuations

G4double RandomSpectrum(G4double ratio, G4double x)
{
  G4double result;
  G4double beta2 = x/(1+x);
  G4complex a = G4complex(epsilon*std::sqrt(1.-ratio*ratio)*beta2, ratio*epsilon*beta2);
  G4complex b = 1. - 1./a;
  G4complex ln = std::log(1./(1.-a));
  G4complex bln = b*ln;

  result = imag(bln);
  result *= cof/pi;

  return result;
}

/////////////////////////////////////////////////////////////////////////
//
// Tamm-Frank spectrum of transparent uniform medium

G4double TammFrank(G4double x)
{
  G4double result;
  G4double beta2 = x/(1+x);
  G4double a = epsilon*beta2;

  if( a < 1.) result = DBL_MIN;
  else        result = cof*(1.-1./a);

  return result;
}





int main()
{

  G4int i, iMax;

  std::ofstream writef("crspectrum.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  G4double betagamma, bg2, r, momentum; // , picr, kcr;
  G4double tf, cr1, cr2, cr3, cr4;
  iMax = 100; // 150;
  
  writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {

    momentum = std::exp(i*0.1)*0.1*MeV;

    betagamma = momentum;    //   /139.57018*MeV;
    bg2 = betagamma*betagamma;

    r   = 0.7;
    cr1 = RandomSpectrum(r,bg2); 

    r   = 0.1;
    cr2 = RandomSpectrum(r,bg2); 

    // picr = RandomSpectrum(r,bg2); 

    //betagamma = momentum/493.677*MeV;
    // bg2 = betagamma*betagamma;
    // kcr = RandomSpectrum(r,bg2); 


    // r   = 0.01;
    // cr3 = RandomSpectrum(r,bg2); 

    r   = 0.001;
    cr4 = RandomSpectrum(r,bg2);

    tf = TammFrank(bg2); 
 
    // G4cout<<momentum<<"\t"<<picr<<"\t"<<kcr<<G4endl;
    // writef<<momentum<<"\t"<<picr<<"\t"<<kcr<<G4endl;

    G4cout<<momentum<<"\t"<<tf<<"\t"<<cr4<<"\t"<<cr2<<"\t"<<cr1<<G4endl;
    writef<<momentum<<"\t"<<tf<<"\t"<<cr4<<"\t"<<cr2<<"\t"<<cr1<<G4endl;


  }

  return 1;
} // end of main
