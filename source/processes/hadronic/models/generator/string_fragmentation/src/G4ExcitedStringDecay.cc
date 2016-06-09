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
// Historic fragment from M.Komogorov; clean-up still necessary @@@
#include "G4ExcitedStringDecay.hh"

G4ExcitedStringDecay::G4ExcitedStringDecay()
{
	theStringDecay=NULL;
}

G4ExcitedStringDecay::G4ExcitedStringDecay(G4VLongitudinalStringDecay * aStringDecay)
: theStringDecay(aStringDecay)
{}

G4ExcitedStringDecay::G4ExcitedStringDecay(const G4ExcitedStringDecay &) : G4VStringFragmentation()
{
} 

G4ExcitedStringDecay::~G4ExcitedStringDecay()
{
}

const G4ExcitedStringDecay & G4ExcitedStringDecay::operator=(const G4ExcitedStringDecay &)
{
  G4Exception("G4ExcitedStringDecay::operator= meant to not be accessable");
  return *this;
}

int G4ExcitedStringDecay::operator==(const G4ExcitedStringDecay &) const
{
  return 0;
}

int G4ExcitedStringDecay::operator!=(const G4ExcitedStringDecay &) const
{
  return 1;
}

G4bool G4ExcitedStringDecay::
EnergyAndMomentumCorrector(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMom)   
    {
    const int    nAttemptScale = 500;
    const double ErrLimit = 1.E-5;
    if (Output->empty())
       return TRUE;
    G4LorentzVector SumMom;
    G4double        SumMass = 0;     
    G4double        TotalCollisionMass = TotalCollisionMom.m(); 
    // Calculate sum hadron 4-momenta and summing hadron mass
    unsigned int cHadron;
    for(cHadron = 0; cHadron < Output->size(); cHadron++)
        {
        SumMom  += Output->operator[](cHadron)->Get4Momentum();
        SumMass += Output->operator[](cHadron)->GetDefinition()->GetPDGMass();
        }
    if (SumMass > TotalCollisionMass) return FALSE;
    SumMass = SumMom.m2();
    if (SumMass < 0) return FALSE;
    SumMass = sqrt(SumMass);

     // Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
    G4ThreeVector Beta = -SumMom.boostVector();
    Output->Boost(Beta);

    // Scale total c.m.s. hadron energy (hadron system mass).
    // It should be equal interaction mass
    G4double Scale = 1;
    for(G4int cAttempt = 0; cAttempt < nAttemptScale; cAttempt++)
        {
        G4double Sum = 0;
        for(cHadron = 0; cHadron < Output->size(); cHadron++)
            {
            G4LorentzVector HadronMom = Output->operator[](cHadron)->Get4Momentum();
            HadronMom.setVect(Scale*HadronMom.vect());
            G4double E = sqrt(HadronMom.vect().mag2() + sqr(Output->operator[](cHadron)->GetDefinition()->GetPDGMass()));
            HadronMom.setE(E);
            Output->operator[](cHadron)->Set4Momentum(HadronMom);
            Sum += E;
            }   
        Scale = TotalCollisionMass/Sum;    
        if (Scale - 1 <= ErrLimit) goto LTRUE;
        }
    
    G4cout << "G4VPartonStringModel::EnergyAndMomentumCorrector"<<G4endl;
    G4cout << "   Increase number of attempts or increase ERRLIMIT"<<G4endl;
 
 LTRUE:   
    // Compute c.m.s. interaction velocity and KTV back boost   
    Beta = TotalCollisionMom.boostVector(); 
    Output->Boost(Beta);
    return TRUE;
    }



