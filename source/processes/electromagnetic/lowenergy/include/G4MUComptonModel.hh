
// *********************************************************************
// |                                                                   |
// |             G4MUComptonModel -- Geant4 Monash University          |
// |                   low energy Compton scattering model.            |
// |             J. M. C. Brown, Monash University, Australia          |
// |                    ## Unpolarised photons only ##                 |
// |                                                                   |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | The following is a Geant4 class to simulate the process of        |
// | bound electron Compton scattering. General code structure is      |
// | based on G4LowEnergyCompton.cc and G4LivermoreComptonModel.cc.    |
// | Algorithms for photon energy, and ejected Compton electron        |
// | direction taken from:                                             |
// |                                                                   |
// | J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin,    |
// | "The Monash University low energy Compton scattering model",      |
// | IEEE Transactions on Nuclear Science, in preparation.             |
// |                                                                   |
// | The author acknowledges the work of the Geant4 collaboration      |
// | in developing the following algorithms that have been employed    |
// | or adapeted for the present software:                             |    
// |                                                                   |
// |  # sampling of photon scattering angle,                           |
// |  # target element selection in composite materials,               |
// |  # target shell selection in element,                             |
// |  # and sampling of bound electron momentum from Compton profiles. |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | History:                                                          |
// | --------                                                          |
// |                                                                   |
// | Nov. 2011 JMCB       - First version                              |
// | Feb. 2012 JMCB       - Migration to Geant4 9.5                    |
// *********************************************************************

#ifndef G4MUComptonModel_h
#define G4MUComptonModel_h 1

#include "G4VEmModel.hh"
#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"

class G4ParticleChangeForGamma;
class G4VCrossSectionHandler;
class G4VAtomDeexcitation;
class G4VEMDataSet;

class G4MUComptonModel : public G4VEmModel
{

public:

  G4MUComptonModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "MUCompton");
  
  virtual ~G4MUComptonModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom( const G4ParticleDefinition*,
                                               G4double kinEnergy, 
                                               G4double Z, 
                                               G4double A=0, 
                                               G4double cut=0,
                                               G4double emax=DBL_MAX );

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:

  G4double lowEnergyLimit;  
  G4double highEnergyLimit; 
  G4bool isInitialised;
  G4int verboseLevel;
  
  G4VEMDataSet* scatterFunctionData;
  G4VCrossSectionHandler* crossSectionHandler;
  
  G4VAtomDeexcitation*    fAtomDeexcitation;

  G4ShellData shellData;
  G4DopplerProfile profileData;

  G4MUComptonModel & operator=(const  G4MUComptonModel &right);
  G4MUComptonModel(const  G4MUComptonModel&);

  
};

//****************************************************************************

#endif
