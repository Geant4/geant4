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
#include "globals.hh"
#include "G4FermiMomentum.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4Fancy3DNucleus.hh"


int main(int argc, char ** argv)
{
  G4int A, Z, nEvents;

  if(argc != 4)
  {
    cout << "Input A and Z: ";
    cin >> A >> Z;
    cout << "Input number of events: ";
    cin >> nEvents;
  }
  else
  {
    A = atoi(argv[1]);
    Z = atoi(argv[2]);
    nEvents = atoi(argv[3]);
  }

  G4VNuclearDensity * density;
  if(A < 17)
    density = new G4NuclearShellModelDensity(A, Z);
  else
    density = new G4NuclearFermiDensity(A, Z);

  G4FermiMomentum fermiMom;
  fermiMom.Init(A, Z);

  G4Fancy3DNucleus nucleus;
  nucleus.Init(A, Z);
  G4double radius = nucleus.GetOuterRadius();
  G4double step = radius/1000.;
  G4double r = 0;
  while(r < 2*radius)
  {
    G4ThreeVector pos(0, 0, r);
    G4double d = density->GetDensity(pos);
    G4double p = fermiMom.GetFermiMomentum(d);
    G4ThreeVector mom = fermiMom.GetMomentum(d);
    cout << p/MeV << " " << (1/MeV)*mom.x() << " " << (1/MeV)*mom.y() << " "
	 << (1/MeV)*mom.z() << " " << (1/MeV)*mom.mag() << " " << r/fermi
	 << endl;
    r += step;
  }
  return 0;
}
