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
#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4FermiMomentum.hh"

#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "globals.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4KineticTrackVectorSTL.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

G4V3DNucleus * the3DNucleus;
G4KineticTrackVectorSTL theTargetList;

void BuildTargetList();
G4ThreeVector GetSpherePoint(G4double r);

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

  G4KineticTrackVectorSTL::iterator iter;
  for(G4int i = 0; i < nEvents; ++i)
  {
// create the nucleus
    the3DNucleus = new G4Fancy3DNucleus;
    the3DNucleus->Init(A, Z);
// create theTargetList
    BuildTargetList();
    G4LorentzVector totalMom(0, 0, 0, 0);
    for(iter = theTargetList.begin(); iter != theTargetList.end(); ++iter)
    {
      G4ThreeVector pos = (1/fermi)*(*iter)->GetPosition();
      G4LorentzVector mom = (1/MeV)*(*iter)->Get4Momentum();
      totalMom+=mom;
      cout << pos.x() << " " << pos.y() << " " << pos.z() << " "
	   << pos.mag() << " " << mom.x() << " " << mom.y() << " "
	   << mom.z() << " " << mom.t() << " " << mom.vect().mag()
	   << " " << mom.mag() << endl;
    }
    cerr << totalMom.x() << " " << totalMom.y() << " "
	 << totalMom.z() << " " << totalMom.t() << " "
	 << totalMom.vect().mag() << " " << totalMom.mag() << " "
	 << (1/MeV) * G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A)
	 << endl;
  }
  return 0;
}

void BuildTargetList()
{
  G4KineticTrackVectorSTL::iterator iter;
  for(iter = theTargetList.begin(); iter != theTargetList.end(); ++iter)
    delete (*iter);

  theTargetList.clear();

  if(!the3DNucleus->StartLoop())
  {
    G4cerr << "G4HadronKineticModel::BuildTargetList(): StartLoop() error!"
	   << G4endl;
    return;
  }
  G4Nucleon * nucleon;
  G4ParticleDefinition * definition;
  G4ThreeVector pos;
  G4LorentzVector mom;
  while((nucleon = the3DNucleus->GetNextNucleon()) != NULL)
  {
    definition = nucleon->GetDefinition();
    pos = nucleon->GetPosition();
    mom = nucleon->GetMomentum();
    G4KineticTrack * kt = new G4KineticTrack(definition, 0., pos, mom);
    theTargetList.push_back(kt);
  }
}






