// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4EvaporationProbability.hh"




G4EvaporationProbability::G4EvaporationProbability(const G4EvaporationProbability &right)
{
 G4Exception("G4EvaporationProbability::copy_constructor meant to not be accessable");
}




const G4EvaporationProbability & G4EvaporationProbability::
operator=(const G4EvaporationProbability &right)
{
  G4Exception("G4EvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4EvaporationProbability::operator==(const G4EvaporationProbability &right) const
{
  return false;
}

G4bool G4EvaporationProbability::operator!=(const G4EvaporationProbability &right) const
{
  return true;
}

G4double G4EvaporationProbability::EmissionProbability(const G4Fragment & fragment, const G4double anEnergy)
{
	G4double EmissionProbability = 0.0;
	G4double MaximalKineticEnergy = anEnergy;
	if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
		EmissionProbability = CalcProbability(fragment,MaximalKineticEnergy);
// 		// Next there is a loop over excited states for this channel summing probabilities
// 		G4double SavedGamma = Gamma;
// 		for (G4int i = 0; i < ExcitationEnergies->length(); i++) {
// 			if (ExcitationSpins->operator()(i) < 0.1) continue;
// 			Gamma = ExcitationSpins->operator()(i)*A; // A is the channel mass number
// 			// substract excitation energies
// 			MaximalKineticEnergy -= ExcitationEnergies->operator()(i);
// 			// update probability
// 			EmissionProbability += CalcProbability(fragment,MaximalKineticEnergy);
// 				EmissionProbability += tmp;
// 			}
// 		// restore Gamma and MaximalKineticEnergy
// 		MaximalKineticEnergy = SavedMaximalKineticEnergy;
// 		Gamma = SavedGamma;
// 		}
	}
	return EmissionProbability;
}

G4double G4EvaporationProbability::CalcProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy)
  // Calculate integrated probability (width) for rvaporation channel
{	
	G4double ResidualA = G4double(fragment.GetA() - theA);
	G4double ResidualZ = G4double(fragment.GetZ() - theZ);
	G4double U = fragment.GetExcitationEnergy();
	
	G4double NuclearMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusMass(theZ,theA);


	G4double delta0 = PairingCorrection(fragment.GetA(),fragment.GetZ());

	G4double SystemEntropy = 2.0*sqrt(theEvapLDPptr->LevelDensityParameter(fragment.GetA(),fragment.GetZ(),U)*(U-delta0));
								  
	// compute the integrated probability of evaporation channel
	G4double RN = 1.5*fermi;

	G4double Alpha = CalcAlphaParam(fragment);
	G4double Beta = CalcBetaParam(fragment);
	
	G4double Rmax = MaximalKineticEnergy;
	G4double a = theEvapLDPptr->LevelDensityParameter(ResidualA,ResidualZ,Rmax);
	G4double GlobalFactor = G4double(Gamma) * (Alpha/(a*a)) *
									(NuclearMass*RN*RN*pow(ResidualA,2./3.))/
									(2.*pi* hbar_Planck*hbar_Planck);
	G4double Term1 = (2.0*Beta*a-3.0)/2.0 + Rmax*a;
	G4double Term2 = (2.0*Beta*a-3.0)*sqrt(Rmax*a) + 2.0*a*Rmax;
	
	G4double ExpTerm1 = 0.0;
	if (SystemEntropy <= 600.0) ExpTerm1 = exp(-SystemEntropy);
	
	G4double ExpTerm2 = 2.*sqrt(a*Rmax) - SystemEntropy;
	if (ExpTerm2 > 700.0) ExpTerm2 = 700.0;
	ExpTerm2 = exp(ExpTerm2);
	
	G4double Width = GlobalFactor*(Term1*ExpTerm1 + Term2*ExpTerm2);
	
	return Width;
}



G4double G4EvaporationProbability::PairingCorrection(const G4int A, const G4int Z) const
{
	const G4double PairingConstant = 12.0*MeV;
	const G4int N = A - Z;
	G4double Pair = (1.0 - G4double(Z) + 2.0*(Z/2)) + (1.0 - G4double(N) + 2.0*(N/2));
	G4double PCorrection = Pair*PairingConstant/sqrt(G4double(A));
	return PCorrection;
}


