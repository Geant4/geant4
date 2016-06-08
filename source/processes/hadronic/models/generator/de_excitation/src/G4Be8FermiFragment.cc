// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4Be8FermiFragment.hh"


G4Be8FermiFragment::G4Be8FermiFragment()
{
}

G4Be8FermiFragment::G4Be8FermiFragment(const G4Be8FermiFragment &right)
{
  G4Exception("G4Be8FermiFragment::copy_constructor meant to not be accessable");
}


G4Be8FermiFragment::~G4Be8FermiFragment()
{
}


const G4Be8FermiFragment & G4Be8FermiFragment::operator=(const G4Be8FermiFragment &right)
{
  G4Exception("G4Be8FermiFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be8FermiFragment::operator==(const G4Be8FermiFragment &right) const
{
  return false;
}

G4bool G4Be8FermiFragment::operator!=(const G4Be8FermiFragment &right) const
{
  return true;
}



G4FragmentVector * G4Be8FermiFragment::GetFragment(const G4LorentzVector & aMomentum)
  // Be8 ----> alpha + alpha 
{
  const G4int NumSubFrag = 2;
  G4double Masses[NumSubFrag];
  G4double Charges[NumSubFrag];
  G4double AtomNum[NumSubFrag];


  Masses[0] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4); // alpha
  Masses[1] = Masses[0]; // alpha
  
  AtomNum[0] = 4;
  AtomNum[1] = 4;

  Charges[0] = 2;
  Charges[1] = 2;

//   G4double AvalKineticE = G4NucleiPropertiesTable::GetMassExcess(Z,A) + ExcitEnergy - // Be8
//     2.0*G4NucleiPropertiesTable::GetMassExcess(2,4); // alphas
  G4double AvalKineticE =  sqrt(aMomentum.e()*aMomentum.e() - 
				aMomentum.vect().mag2())  -// Be8
    2.0*Masses[0]; // alphas


  RWTPtrOrderedVector<G4LorentzVector> * SubFragsMomentum =
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
