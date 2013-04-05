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
#include "G4Neutron.hh"
#include "G4Nucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4NucleiPropertiesTable.hh"

main()
{
  G4Nucleus theNucleus;
  theNucleus.SetParameters(12, 6);
  G4double theBaseA = 12;
  G4int theBaseZ = 6;
  G4double eps = 0.000001;
  G4double targetMass = ( G4NucleiPropertiesTable::GetAtomicMass(theBaseZ+eps, theBaseA+eps)-
                          theBaseZ* G4Electron::Electron()->GetPDGMass() ) /
                          G4Neutron::Neutron()->GetPDGMass();
  
  G4ReactionProduct theNeutron( G4Neutron::NeutronDefinition() );
  G4ReactionProduct theTarget;
  
  G4int N=1000;
  G4int i, j;
  G4double start = 10000.*eV;
  for(j=0; j<100; j++)
  {
    G4ThreeVector aMom (start, 0., 0.);
    start *=1.2;
    G4double mean = 0;
    G4double eKin = aMom*aMom/(2.*G4Neutron::NeutronDefinition()->GetPDGMass());
    for(i=0; i<N; i++)
    {
      theNeutron.SetMomentum( aMom );
      theNeutron.SetKineticEnergy( eKin );
      theTarget = theNucleus.GetThermalNucleus(targetMass);
      theNeutron.Lorentz(theNeutron, -1*theTarget);
      G4double eKinetic = theNeutron.GetKineticEnergy();
      G4double velocity = sqrt(2.*eKinetic/G4Neutron::NeutronDefinition()->GetPDGMass());
      G4cout << "Velocity = "<<velocity<<endl;
      mean += velocity;
    }
    mean /= G4double(N);
    G4cout << "E, Mean = "<<eKin<<" "<<mean<<endl;
  }
}
