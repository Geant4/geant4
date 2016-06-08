//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4He5FermiFragment.cc,v 1.6 2001/08/01 17:05:31 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4He5FermiFragment.hh"


G4He5FermiFragment::G4He5FermiFragment()
{
}

G4He5FermiFragment::G4He5FermiFragment(const G4He5FermiFragment &right)
{
    G4Exception("G4He5FermiFragment::copy_constructor meant to not be accessable");
}


G4He5FermiFragment::~G4He5FermiFragment()
{
}


const G4He5FermiFragment & G4He5FermiFragment::operator=(const G4He5FermiFragment &right)
{
    G4Exception("G4He5FermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4He5FermiFragment::operator==(const G4He5FermiFragment &right) const
{
    return false;
}

G4bool G4He5FermiFragment::operator!=(const G4He5FermiFragment &right) const
{
    return true;
}




G4FragmentVector * G4He5FermiFragment::GetFragment(const G4LorentzVector & aMomentum)
    // He5 ----> alpha + neutron
{
    const G4int NumSubFrag = 2;
    G4double Masses[NumSubFrag];
    G4double Charges[NumSubFrag];
    G4double AtomNum[NumSubFrag];


    Masses[0] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4); // alpha
    Masses[1] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(0,1); // neutron
  
    AtomNum[0] = 4;
    AtomNum[1] = 1;

    Charges[0] = 2;
    Charges[1] = 0;

    //   G4double AvalKineticE = G4NucleiPropertiesTable::GetMassExcess(Z,A) + ExcitEnergy - // He5
    //     G4NucleiPropertiesTable::GetMassExcess(0,1) - // neutron
    //     G4NucleiPropertiesTable::GetMassExcess(2,4); // alpha
    G4double AvalKineticE = sqrt(aMomentum.e()*aMomentum.e() - 
				 aMomentum.vect().mag2())  - // He5
	Masses[1] - // neutron
	Masses[0]; // alpha

    G4std::deque<G4LorentzVector*> * SubFragsMomentum =
	FragmentsMomentum(AvalKineticE, NumSubFrag,Masses);
  


    G4FragmentVector * theResult = new G4FragmentVector;

    for (G4int i = 0; i < NumSubFrag; i++) {

	// Lorentz boost
	SubFragsMomentum->operator[](i)->boost(aMomentum.boostVector());
    
	theResult->push_back(new G4Fragment(AtomNum[i],Charges[i],*(SubFragsMomentum->operator[](i))));
    }
  
    //  SubFragsMomentum->clearAndDestroy();
    while (!SubFragsMomentum->empty()) {
	delete SubFragsMomentum->back();
	SubFragsMomentum->pop_back();
    }
    delete SubFragsMomentum;

    return theResult;
}
