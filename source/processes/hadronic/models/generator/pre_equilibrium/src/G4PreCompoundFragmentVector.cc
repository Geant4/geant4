#include "G4PreCompoundFragmentVector.hh"

#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundProton.hh"
#include "G4PreCompoundDeuteron.hh"
#include "G4PreCompoundTriton.hh"
#include "G4PreCompoundHe3.hh"
#include "G4PreCompoundAlpha.hh"

G4PreCompoundFragmentVector::G4PreCompoundFragmentVector() :
	TotalEmissionProbability(0.0)
{
  // neutron 
  theChannels.insert(new G4PreCompoundNeutron());
  // proton
  theChannels.insert(new G4PreCompoundProton());
  // deuterium
  theChannels.insert(new G4PreCompoundDeuteron());
  // triton
  theChannels.insert(new G4PreCompoundTriton());
  // helium3
  theChannels.insert(new G4PreCompoundHe3());
  // alpha
  theChannels.insert(new G4PreCompoundAlpha());
}


G4PreCompoundFragmentVector::~G4PreCompoundFragmentVector()
{
	theChannels.clearAndDestroy();
}

const G4PreCompoundFragmentVector & G4PreCompoundFragmentVector::operator=(const G4PreCompoundFragmentVector &right)
{
    G4Exception("G4PreCompoundFragmentVector::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundFragmentVector::operator==(const G4PreCompoundFragmentVector &right) const
{
  return false;
}

G4bool G4PreCompoundFragmentVector::operator!=(const G4PreCompoundFragmentVector &right) const
{
  return true;
}



G4double G4PreCompoundFragmentVector::CalculateProbabilities(const G4Fragment & aFragment)
{
	TotalEmissionProbability = 0.0;
	for (G4int i = 0; i < theChannels.entries(); i++) {
		theChannels(i)->CalcExcitonLevelDensityRatios(aFragment.GetNumberOfExcitons(),
																	 aFragment.GetNumberOfParticles());	
		theChannels(i)->CalcCondensationProbability(aFragment.GetA());
		// Calculate emission probailities
		if (aFragment.GetNumberOfParticles() <= theChannels(i)->GetA()-0.01) {
			// if number of particles less than a fragment atomic number 
			// set probability to emit a fragment 0
			theChannels(i)->SetEmissionProbability(0.0);
		} else if (aFragment.GetNumberOfExcitons() <= theChannels(i)->GetA()+0.01 && 
															aFragment.GetNumberOfExcitons() != 1) {
			theChannels(i)->SetEmissionProbability(0.0);
		} else if (aFragment.GetNumberOfCharged() <= theChannels(i)->GetZ()-0.01) {
			// if number of charged particles (protons) is less than charge of fragment
			// set probability to emit a fragment 0
			theChannels(i)->SetEmissionProbability(0.0);
		} else if (theChannels(i)->GetMaximalKineticEnergy() <= 0.0) {
			// if the energy threshold for emitted fragment is less or equal 0
			// set probability to emit a fragment 0
			theChannels(i)->SetEmissionProbability(0.0);
		} else {
			// Compute total (integrated over kinetic energy) emission 
			// probability of a fragment and
			// Summing channel emission probabilities
			TotalEmissionProbability += theChannels(i)->CalcEmissionProbability(aFragment);
		}
	}
	return TotalEmissionProbability;
}


G4VPreCompoundFragment * G4PreCompoundFragmentVector::ChooseFragment(void)
{
	const G4int NumOfFrags = theChannels.entries();
	G4double * running = new G4double[NumOfFrags];
	running[0] = theChannels(0)->GetEmissionProbability();
        G4int i;
	for (i = 1; i < NumOfFrags; i++) {
		running[i]=running[i-1]+theChannels(i)->GetEmissionProbability();
	}
	
	// Choose an emission channel
	G4double aChannel = G4UniformRand()*TotalEmissionProbability;
	G4int ChosenChannel = -1;
	for (i = 0; i < NumOfFrags; i++) {
		if (aChannel <= running[i]) {
			ChosenChannel = i;
			break;
		}
	}
	delete [] running;
	if (ChosenChannel < 0) 
		G4Exception("G4PreCompoundFragmentVector::ChooseFragment: I can't determine a channel");
	
	return theChannels(ChosenChannel);
}


