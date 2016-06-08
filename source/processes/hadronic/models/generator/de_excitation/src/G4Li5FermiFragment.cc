// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4Li5FermiFragment.hh"


G4Li5FermiFragment::G4Li5FermiFragment()
{
}

G4Li5FermiFragment::G4Li5FermiFragment(const G4Li5FermiFragment &right)
{
  G4Exception("G4Li5FermiFragment::copy_constructor meant to not be accessable");
}


G4Li5FermiFragment::~G4Li5FermiFragment()
{
}


const G4Li5FermiFragment & G4Li5FermiFragment::operator=(const G4Li5FermiFragment &right)
{
  G4Exception("G4Li5FermiFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4Li5FermiFragment::operator==(const G4Li5FermiFragment &right) const
{
  return false;
}

G4bool G4Li5FermiFragment::operator!=(const G4Li5FermiFragment &right) const
{
  return true;
}


G4FragmentVector * G4Li5FermiFragment::GetFragment(const G4LorentzVector & aMomentum)
  // Li5 ----> alpha + proton
{
  const G4int NumSubFrag = 2;
  G4double Masses[NumSubFrag];
  G4double Charges[NumSubFrag];
  G4double AtomNum[NumSubFrag];


  Masses[0] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4); // alpha
  Masses[1] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1); // proton
  
  AtomNum[0] = 4;
  AtomNum[1] = 1;

  Charges[0] = 2;
  Charges[1] = 1;

//   G4double AvalKineticE = G4NucleiPropertiesTable::GetMassExcess(Z,A) + ExcitEnergy - // Li5
//     G4NucleiPropertiesTable::GetMassExcess(1,1) - // proton
//     G4NucleiPropertiesTable::GetMassExcess(2,4); // alpha
  G4double AvalKineticE =  sqrt(aMomentum.e()*aMomentum.e() - 
				aMomentum.vect().mag2())  - // Li5
		   Masses[0] - // proton
		  Masses[1]; // alpha

  G4RWTPtrOrderedVector<G4LorentzVector> * SubFragsMomentum =
    FragmentsMomentum(AvalKineticE, NumSubFrag,Masses);



  G4FragmentVector * theResult = new G4FragmentVector;

  for (G4int i = 0; i < NumSubFrag; i++) {
    
    // Lorentz boost
    SubFragsMomentum->at(i)->boost(aMomentum.boostVector());

    theResult->insert(new G4Fragment(AtomNum[i],Charges[i],*(SubFragsMomentum->at(i))));
  }

  SubFragsMomentum->clearAndDestroy();
  delete SubFragsMomentum;

  return theResult;
}
