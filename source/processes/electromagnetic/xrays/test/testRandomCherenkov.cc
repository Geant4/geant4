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


G4double epsilon = 1.2;

using namespace std;

G4double randomSpectrum(G4double ratio, G4double x)
{
  G4double order = 6.*x;
  return 1. - std::exp(-order);
}



G4double TammFrank(G4double x)
{
  G4double order = 6.*x;
  return 1. - std::exp(-order);
}





int main()
{

  G4int i, iMax;

  std::ofstream writef("crspectrum.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  G4double betagamma, bg2, cr1, cr2, cr3, cr4;

  iMax = 100;
  
  // writef<<iMax<<G4endl;

  for( i = 0; i < iMax; i++ )
  {

    betagamma = std::exp(i*0.1)*0.1;
    G4cout<<"betagamma = "<<betagamma<<G4endl;
  }




  return 1;
} // end of main
