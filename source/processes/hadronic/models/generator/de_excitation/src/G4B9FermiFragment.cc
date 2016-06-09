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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4B9FermiFragment.cc,v 1.10 2003/06/16 17:06:14 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4B9FermiFragment.hh"


G4B9FermiFragment::G4B9FermiFragment()
{
}

G4B9FermiFragment::G4B9FermiFragment(const G4B9FermiFragment &) : G4UnstableFermiFragment()
{
    G4Exception("G4B9FermiFragment::copy_constructor meant to not be accessable");
}


G4B9FermiFragment::~G4B9FermiFragment()
{
}


const G4B9FermiFragment & G4B9FermiFragment::operator=(const G4B9FermiFragment &)
{
    G4Exception("G4B9FermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4B9FermiFragment::operator==(const G4B9FermiFragment &) const
{
    return false;
}

G4bool G4B9FermiFragment::operator!=(const G4B9FermiFragment &) const
{
    return true;
}



G4FragmentVector * G4B9FermiFragment::GetFragment(const G4LorentzVector & aMomentum)
    // B9 ----> alpha + alpha + proton
{
    const G4int NumSubFrag = 3;
    G4double Masses[NumSubFrag];
    G4double Charges[NumSubFrag];
    G4double AtomNum[NumSubFrag];


    Masses[0] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4); // alpha
    Masses[1] = Masses[0]; // alpha
    Masses[2] = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1); // proton
  
    AtomNum[0] = 4;
    AtomNum[1] = 4;
    AtomNum[2] = 1;

    Charges[0] = 2;
    Charges[1] = 2;
    Charges[2] = 1;

    //   G4double AvalKineticE = G4NucleiPropertiesTable::GetMassExcess(Z,A) + ExcitEnergy - // B9
    //     G4NucleiPropertiesTable::GetMassExcess(1,1) - // proton
    //     2.0*G4NucleiPropertiesTable::GetMassExcess(2,4);
    G4double AvalKineticE =  sqrt(aMomentum.e()*aMomentum.e() - 
				  aMomentum.vect().mag2())  - // B9
	Masses[2] - // proton
	2.0*Masses[0];


    std::deque<G4LorentzVector*> * SubFragsMomentum =
	FragmentsMomentum(AvalKineticE, NumSubFrag,Masses);

    G4FragmentVector * theResult = new G4FragmentVector;

    for (G4int i = 0; i < NumSubFrag; i++) {
    
	// Lorentz boost
	SubFragsMomentum->operator[](i)->boost(aMomentum.boostVector());
    
	theResult->push_back(new G4Fragment(static_cast<G4int>(AtomNum[i]),
					    static_cast<G4int>(Charges[i]),
					    *(SubFragsMomentum->operator[](i))));
    }

    //  SubFragsMomentum->clearAndDestroy();
    while (!SubFragsMomentum->empty()) {
	delete SubFragsMomentum->back();
	SubFragsMomentum->pop_back();
    }
    delete SubFragsMomentum;

    return theResult;
}
