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
#include "G4VKM_NuclearDensity.hh"
#include "G4KM_NuclearFermiDensity.hh"
#include "G4KM_NuclearShellModelDensity.hh"
#include "G4VNuclearDensity.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4Fancy3DNucleus.hh"



int main(int argc, char ** argv)
{
  // Get nucleus
  G4int A, Z;
  if(argc != 3)
  {
    cerr << "Wrong arguments. Please, enter A and Z for the nucleus: ";
    cin >> A >> Z;
  }
  else
  {
    A = atoi(argv[1]);
    Z = atoi(argv[2]);
  }

  G4VNuclearDensity *g4Density;
  G4VKM_NuclearDensity * nuclearDensity;
  if(A < 17)
  {
    g4Density = new G4NuclearShellModelDensity(A, Z);
    nuclearDensity = new G4KM_NuclearShellModelDensity(A);
  }
  else
  {
    g4Density = new G4NuclearFermiDensity(A, Z);
    nuclearDensity = new G4KM_NuclearFermiDensity(A);
  }

  G4Fancy3DNucleus nucleus;
  nucleus.Init(A, Z);
  G4double radius = nucleus.GetOuterRadius();
  G4double step = radius/1000.;
  G4double r = 0;
  while(r < 2*radius)
    {
      G4ThreeVector pos(0, 0, r);
      cout << nuclearDensity->GetDensity(pos)*fermi*fermi*fermi << " "
	   << g4Density->GetDensity(pos)*fermi*fermi*fermi << " "
	   << nuclearDensity->GetDeriv(pos)*fermi*fermi*fermi*fermi << " "
	   << g4Density->GetDeriv(pos)*fermi*fermi*fermi*fermi << " "
	   << r/fermi << endl;
      r += step;
    }

  return 0;
}

