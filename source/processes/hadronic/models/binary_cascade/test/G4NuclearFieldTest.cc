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
#include "G4VNuclearField.hh"
#include "G4AntiProtonField.hh"
#include "G4PionZeroField.hh"
#include "G4ProtonField.hh"
#include "G4KaonMinusField.hh"
#include "G4SigmaMinusField.hh"
#include "G4KaonPlusField.hh"
#include "G4SigmaPlusField.hh"
#include "G4KaonZeroField.hh"
#include "G4SigmaZeroField.hh"
#include "G4NeutronField.hh"
#include "G4PionMinusField.hh"
#include "G4VNuclearField.hh"
#include "G4PionPlusField.hh"
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

  G4V3DNucleus * nucleus = new G4Fancy3DNucleus;
  nucleus->Init(A, Z);
  G4VNuclearField * pro = new G4ProtonField(nucleus);
  G4VNuclearField * neu = new G4NeutronField(nucleus);
  G4VNuclearField * pbar = new G4AntiProtonField(nucleus);
  G4VNuclearField * pip = new G4PionPlusField(nucleus);
  G4VNuclearField * pi0 = new G4PionZeroField(nucleus);
  G4VNuclearField * pim = new G4PionMinusField(nucleus);
  G4VNuclearField * kp = new G4KaonPlusField(nucleus);
  G4VNuclearField * k0 = new G4KaonZeroField(nucleus);
  G4VNuclearField * km = new G4KaonMinusField(nucleus);
  G4VNuclearField * sp = new G4SigmaPlusField(nucleus);
  G4VNuclearField * s0 = new G4SigmaZeroField(nucleus);
  G4VNuclearField * sm = new G4SigmaMinusField(nucleus);

  G4double radius = nucleus->GetOuterRadius();
  G4double step = radius/1000.;
  G4double r = 0;
  
  G4cout << "C Nucradius proton Neutron Pbar Pi+ Pi0 Pi- K+ K0 K- Sigma+ Sigma0 Sigma-" << G4endl;
  G4cout << "C " << radius << " " 
  	   << pro->GetBarrier()/MeV << " "
	   << neu->GetBarrier()/MeV << " "
	   << pbar->GetBarrier()/MeV << " "
	   << pip->GetBarrier()/MeV << " "
	   << pi0->GetBarrier()/MeV << " "
	   << pim->GetBarrier()/MeV << " "
	   << kp->GetBarrier()/MeV << " "
	   << k0->GetBarrier()/MeV << " "
	   << km->GetBarrier()/MeV << " "
	   << sp->GetBarrier()/MeV << " "
	   << s0->GetBarrier()/MeV << " "
	   << sm->GetBarrier()/MeV << " "
	   << G4endl;

  G4cout << "C radius proton Neutron Pbar Pi+ Pi0 Pi- K+ K0 K- Sigma+ Sigma0 Sigma-" << G4endl;
  while(r < 1.1*radius)
    {
      G4ThreeVector pos(0, 0, r);
      G4cout << r/fermi << " "
           << pro->GetField(pos)/MeV << " "
	   << neu->GetField(pos)/MeV << " "
	   << pbar->GetField(pos)/MeV << " "
	   << pip->GetField(pos)/MeV << " "
	   << pi0->GetField(pos)/MeV << " "
	   << pim->GetField(pos)/MeV << " "
	   << kp->GetField(pos)/MeV << " "
	   << k0->GetField(pos)/MeV << " "
	   << km->GetField(pos)/MeV << " "
	   << sp->GetField(pos)/MeV << " "
	   << s0->GetField(pos)/MeV << " "
	   << sm->GetField(pos)/MeV << " "
	   << G4endl;
      r += step;
    }

  return 0;
}

