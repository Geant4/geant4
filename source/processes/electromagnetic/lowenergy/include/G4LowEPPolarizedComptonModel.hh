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
// ********************************************************************
// *********************************************************************
// |                                                                   |
// |      G4LowEPPolarizedComptonModel-- Geant4 Monash University      |
// |         polarised low energy Compton scattering model.            |
// |          J. M. C. Brown, Monash University, Australia             |
// |                                                                   |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | The following is a Geant4 class to simulate the process of        |
// | bound electron Compton scattering. General code structure is      |
// | based on G4LowEnergyCompton.cc and                                |
// | G4LivermorePolarizedComptonModel.cc.                              |
// | Algorithms for photon energy, and ejected Compton electron        |
// | direction taken from:                                             |
// |                                                                   |
// | J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin,    |
// | "A low energy bound atomic electron Compton scattering model      |
// |  for Geant4", NIMB, Vol. 338, 77-88, 2014.                        |
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
// | Jan. 2015 JMCB       - 1st Version based on G4LowEPComptonModel   |
// | Feb. 2016 JMCB       - Geant4 10.2 FPE fix for bug 1676           |
// | Nov. 2016 JMCB       - Polarisation tracking fix in collaboration |
// |                        of Dr. Merlin Reynaard Kole,               |
// |                        University of Geneva                       |
// |                                                                   |
// *********************************************************************



#ifndef G4LowEPPolarizedComptonModel_h
#define G4LowEPPolarizedComptonModel_h 1

#include "G4VEmModel.hh"
#include "G4LPhysicsFreeVector.hh"
#include <limits>
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4ForceCondition.hh"


class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;
class G4ShellData;
class G4DopplerProfile;

class G4LowEPPolarizedComptonModel : public G4VEmModel
{

public:

  G4LowEPPolarizedComptonModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "LowEPComptonModel");
  
  virtual ~G4LowEPPolarizedComptonModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*,
                               G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

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

private:

  void ReadData(size_t Z, const char* path = 0);

  G4double ComputeScatteringFunction(G4double x, G4int Z);

  G4LowEPPolarizedComptonModel & operator=(const  G4LowEPPolarizedComptonModel &right);
  G4LowEPPolarizedComptonModel(const  G4LowEPPolarizedComptonModel&);

  G4ThreeVector GetRandomPolarization(G4ThreeVector& direction0); // Random Polarization
  G4ThreeVector GetPerpendicularPolarization(const G4ThreeVector& direction0, const G4ThreeVector& polarization0) const;

  G4ThreeVector SetNewPolarization(G4double LowEPPCepsilon, G4double sinT2,
                                   G4double phi, G4double cosTheta);
  G4double SetPhi(G4double, G4double);

  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1,
                         G4ThreeVector& polarization0, G4ThreeVector& polarization1);

  void SystemOfRefChangeElect(G4ThreeVector& pdirection, G4ThreeVector& edirection,
                         G4ThreeVector& ppolarization);


  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a); 
  G4bool isInitialised;
  G4int verboseLevel;
 
  //G4double lowestEnergy;
 
  G4ParticleChangeForGamma* fParticleChange;
  G4VAtomDeexcitation*      fAtomDeexcitation;

  static G4ShellData*       shellData;
  static G4DopplerProfile*  profileData;

  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100];

  static const G4double ScatFuncFitParam[101][9];
  
};

//****************************************************************************

#endif
