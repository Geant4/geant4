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
// $Id: G4PenelopeBremsstrahlungModel.cc,v 1.8 2010-11-25 09:44:05 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
// --------
// 05 Dec 2008   L Pandola    Migration from process to model
// 25 Mar 2008   L Pandola    Fixed .unit() call
// 16 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - added MinEnergyCut method
//                  - do not change track status
// 14 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they are checked later on
//
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeBremsstrahlungContinuous.hh"
#include "G4PenelopeBremsstrahlungAngular.hh"
#include "G4eBremsstrahlungSpectrum.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4DataVector.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4LogLogInterpolation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungModel::G4PenelopeBremsstrahlungModel(const G4ParticleDefinition*,
							     const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),energySpectrum(0),
   angularData(0),stoppingPowerData(0),crossSectionHandler(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  
  verboseLevel= 0;
   
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
   
  //These vectors do not change when materials or cut change.
  //Therefore I can read it at the constructor
  angularData = new std::map<G4int,G4PenelopeBremsstrahlungAngular*>;
 
  //These data do not depend on materials and cuts.
  G4DataVector eBins;
 
  eBins.push_back(1.0e-12);
  eBins.push_back(0.05);
  eBins.push_back(0.075);
  eBins.push_back(0.1);
  eBins.push_back(0.125);
  eBins.push_back(0.15);
  eBins.push_back(0.2);
  eBins.push_back(0.25);
  eBins.push_back(0.3);
  eBins.push_back(0.35);
  eBins.push_back(0.40);
  eBins.push_back(0.45);
  eBins.push_back(0.50);
  eBins.push_back(0.55);
  eBins.push_back(0.60);
  eBins.push_back(0.65);
  eBins.push_back(0.70);
  eBins.push_back(0.75);
  eBins.push_back(0.80);
  eBins.push_back(0.85);
  eBins.push_back(0.90);
  eBins.push_back(0.925);
  eBins.push_back(0.95);
  eBins.push_back(0.97);
  eBins.push_back(0.99);
  eBins.push_back(0.995);
  eBins.push_back(0.999);
  eBins.push_back(0.9995);
  eBins.push_back(0.9999);
  eBins.push_back(0.99995);
  eBins.push_back(0.99999);
  eBins.push_back(1.0);
 
  const G4String dataName("/penelope/br-sp-pen.dat");
  energySpectrum = new G4eBremsstrahlungSpectrum(eBins,dataName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungModel::~G4PenelopeBremsstrahlungModel()
{
  if (crossSectionHandler)
    delete crossSectionHandler;
  
  if (energySpectrum) 
    delete energySpectrum;
  
  if (angularData)
    {
      std::map <G4int,G4PenelopeBremsstrahlungAngular*>::iterator i;
      for (i=angularData->begin();i != angularData->end();i++)
	if (i->second) delete i->second;
      delete angularData;
    }
  
  if (stoppingPowerData)
    {
      std::map <std::pair<G4int,G4double>,G4PenelopeBremsstrahlungContinuous*>::iterator j;
      for (j=stoppingPowerData->begin();j != stoppingPowerData->end();j++)
	if (j->second) delete j->second;
      delete stoppingPowerData;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::Initialise(const G4ParticleDefinition* particle,
					       const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeBremsstrahlungModel::Initialise()" << G4endl;
  
  // Delete everything, but angular data (do not depend on cuts)
  if (crossSectionHandler)
    {
      crossSectionHandler->Clear();
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }

  if (stoppingPowerData)
    {
      std::map <std::pair<G4int,G4double>,G4PenelopeBremsstrahlungContinuous*>::iterator j;
      for (j=stoppingPowerData->begin();j != stoppingPowerData->end();j++)
	if (j->second) 
	  { 
	    delete j->second;
	    j->second = 0;
	  }
      delete stoppingPowerData;
      stoppingPowerData = 0;
    }
  
  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Clear();
  //
  if (particle==G4Electron::Electron())
      crossSectionHandler->LoadData("brem/br-cs-");
  else
    crossSectionHandler->LoadData("penelope/br-cs-pos-"); //cross section for positrons
    
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata = 
    crossSectionHandler->BuildMeanFreePathForMaterials();
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;   

  if (verboseLevel > 2)
    G4cout << "Loaded cross section files for PenelopeBremsstrahlungModel" << G4endl;
  
  if (verboseLevel > 0) {
    G4cout << "Penelope Bremsstrahlung model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / keV << " keV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  } 

  //This has to be invoked AFTER the crossSectionHandler has been created, 
  //because it makes use of ComputeCrossSectionPerAtom()
  InitialiseElementSelectors(particle,cuts);

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
						     const G4MaterialCutsCouple*)
{
  return 250.*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
								   G4double kinEnergy,
								   G4double Z,
								   G4double,
								   G4double cutEnergy,
								   G4double)
{
  // Penelope model to calculate cross section for hard bremsstrahlung emission 
  // (gamma energy > cutEnergy).
  //
  // The total bremsstrahlung cross section is read from database, following data 
  // reported in the EEDL library
  //  D.E.Cullen et al., Report UCRL-50400 (Lawrence Livermore National Laboratory) (1989)
  // The probability to have photon emission above a given threshold is calculated 
  // analytically using the differential cross section model dSigma/dW = F(x)/x, where 
  // W is the outgoing photon energy and x = W/E is the ratio of the photon energy to the 
  // incident energy. The function F(x) is tabulated (for all elements) using 32 points in x 
  // ranging from 1e-12 to 1. Data are derived from 
  //  S.M.Seltzer and M.J.Berger, At.Data Nucl.Data Tables 35,345 (1986)
  // Differential cross sections for electrons and positrons dSigma/dW are assumed to scale 
  // with a function S(Z,E) which does not depend on W; therefore, only overall cross section 
  // changes but not the shape of the photon energy spectrum.
  //

  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopeBremsstrahlungModel" << G4endl;
  
  G4int iZ = (G4int) Z;

  // VI - not needed in run time
  // if (!crossSectionHandler)
  //  {
  //    G4cout << "G4PenelopeBremsstrahlungModel::ComputeCrossSectionPerAtom" << G4endl;
  //    G4cout << "The cross section handler is not correctly initialized" << G4endl;
  //    G4Exception();
  //  }
  G4double totalCs = crossSectionHandler->FindValue(iZ,kinEnergy);
  G4double cs = totalCs * energySpectrum->Probability(iZ,cutEnergy,kinEnergy,kinEnergy);

  if (verboseLevel > 2)
    {
      G4cout << "Bremsstrahlung cross section at " << kinEnergy/MeV << " MeV for Z=" << Z <<
	" and energy > " << cutEnergy/keV << " keV --> " << cs/barn << " barn" << G4endl;
      G4cout << "Total bremsstrahlung cross section at " << kinEnergy/MeV << " MeV for Z=" << 
      Z << " --> " << totalCs/barn << " barn" << G4endl;
    }
  return cs;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double 
G4PenelopeBremsstrahlungModel::ComputeDEDXPerVolume(const G4Material* theMaterial,
						    const G4ParticleDefinition* theParticle,
						    G4double kineticEnergy,
						    G4double cutEnergy)
{
  // Penelope model to calculate the stopping power (in [Energy]/[Length]) for soft 
  // bremsstrahlung emission (gamma energy < cutEnergy).
  //
  // The actual calculation is performed by the helper class 
  // G4PenelopeBremsstrahlungContinuous and its method CalculateStopping(). Notice:
  // CalculateStopping() gives the stopping cross section, namely the first momentum of 
  // dSigma/dW, restricted to W < cut (W = gamma energy) This is dimensionally:
  //  [Energy]*[Surface]
  // The calculation is performed by interpolation (in E = incident energy and 
  // x=W/E) from the tabulated data derived from
  //  M.J.Berger and S.M.Seltzer, Report NBSIR 82-2550 (National Bureau of Standards) (1982);
  // for electrons.
  // For positrons, dSigma/dW are assumed to scale with a function S(Z,E) with respect to electrons. 
  // An analytical approximation for the scaling function S(Z,E) is given in
  //  L.Kim et al., Phys.Rev.A 33,3002 (1986)
  //
  if (!stoppingPowerData)
    stoppingPowerData = new std::map<std::pair<G4int,G4double>,
      G4PenelopeBremsstrahlungContinuous*>;
  
  const G4ElementVector* theElementVector = theMaterial->GetElementVector();
  const G4double* theAtomicNumDensityVector = theMaterial->GetAtomicNumDensityVector();

  G4double sPower = 0.0;

  //Loop on the elements of the material
  for (size_t iel=0;iel<theMaterial->GetNumberOfElements();iel++)
    {
      G4int iZ = (G4int) ((*theElementVector)[iel]->GetZ());
      G4PenelopeBremsstrahlungContinuous* theContinuousCalculator = 
	GetStoppingPowerData(iZ,cutEnergy,theParticle);
      sPower += theContinuousCalculator->CalculateStopping(kineticEnergy)*
	theAtomicNumDensityVector[iel];
    }
  
   if (verboseLevel > 2)
    {
      G4cout << "Bremsstrahlung stopping power at " << kineticEnergy/MeV 
	     << " MeV for material " << theMaterial->GetName() 
	     << " and energy < " << cutEnergy/keV << " keV --> " 
	     << sPower/(keV/mm) << " keV/mm" << G4endl;
    }

  return sPower;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						      const G4MaterialCutsCouple* couple,
						      const G4DynamicParticle* aDynamicParticle,
						      G4double cutG,G4double)
{
  // Penelope model to sample the final state for hard bremsstrahlung emission 
  // (gamma energy < cutEnergy).
  // The energy distributionof the emitted photons is sampled according to the F(x) 
  // function tabulated in the database from 
  //  S.M.Seltzer and M.J.Berger, At.Data Nucl.Data Tables 35,345 (1986)
  // The database contains the function F(x) (32 points) for 57 energies of the
  // incident electron between 1 keV and 100 GeV. For other primary energies,
  // logarithmic interpolation is used to obtain the values of the function F(x).
  // The double differential cross section dSigma/(dW dOmega), with W=photon energy,
  // is described as 
  //  dSigma/(dW dOmega) = dSigma/dW * p(Z,E,x,cosTheta)
  // where the shape function p depends on atomic number, incident energy and x=W/E.
  // Numerical values of the shape function, calculated by partial-waves methods, have been 
  // reported in 
  //  L.Kissel et al., At.Data Nucl.Data.Tab. 28,381 (1983); 
  // for Z=2,8,13,47,79 and 92; E=1,5,10,50,100 and 500 keV; x=0,0.6,0.8 and 0.95. The 
  // function p(Z,E,x,cosTheta) is approximated by a Lorentz-dipole function as reported in 
  //  Penelope - A Code System for Monte Carlo Simulation of Electron and Photon Transport, 
  //  Workshop Proceedings Issy-les-Moulineaux, France, 5-7 November 2001, AEN-NEA;
  // The analytical function contains two adjustable parameters that are obtained by fitting 
  // the complete solution from 
  //  L.Kissel et al., At.Data Nucl.Data.Tab. 28,381 (1983); 
  // This allows the evaluation of p(Z,E,x,cosTheta) for any choice of Z, E and x. The actual 
  // sampling of cos(theta) is performed in the helper class
  //  G4PenelopeBremsstrahlungAngular, method ExtractCosTheta()
  // Energy and direction of the primary particle are updated according to energy-momentum 
  // conservation. For positrons, it is sampled the same final state as for electrons.
  //
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4PenelopeBremsstrahlungModel" << G4endl;

  G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
  const G4ParticleDefinition* theParticle = aDynamicParticle->GetDefinition();

  if (kineticEnergy <= fIntrinsicLowEnergyLimit)
   {
     fParticleChange->SetProposedKineticEnergy(0.);
     fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy);
     return ;
   }
 
  G4ParticleMomentum particleDirection0 = aDynamicParticle->GetMomentumDirection();
  //This is the momentum
  G4ThreeVector initialMomentum =  aDynamicParticle->GetMomentum();
 
  //One can use Vladimir's selector! 
  if (verboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;
  // atom can be selected effitiantly if element selectors are initialised
  const G4Element* anElement = SelectRandomAtom(couple,theParticle,kineticEnergy);
  G4int iZ = (G4int) anElement->GetZ();
  if (verboseLevel > 2)
    G4cout << "Selected " << anElement->GetName() << G4endl;
  //

  //Sample gamma's energy according to the spectrum
  G4double gammaEnergy = energySpectrum->SampleEnergy(iZ,cutG,kineticEnergy,kineticEnergy);
  
  //Now sample cosTheta for the Gamma
  G4double cosThetaPrimary = GetAngularDataForZ(iZ)->ExtractCosTheta(kineticEnergy,gammaEnergy);
  
  G4double residualPrimaryEnergy = kineticEnergy-gammaEnergy;
  if (residualPrimaryEnergy < 0)
    {
      //Ok we have a problem, all energy goes with the gamma
      gammaEnergy += residualPrimaryEnergy;
      residualPrimaryEnergy = 0.0;
    }

  //Get primary kinematics
  G4double sinTheta = std::sqrt(1. - cosThetaPrimary*cosThetaPrimary);
  G4double phi  = twopi * G4UniformRand(); 
  G4ThreeVector gammaDirection1(sinTheta* std::cos(phi),
				sinTheta* std::sin(phi),
				cosThetaPrimary);
  
  gammaDirection1.rotateUz(particleDirection0);
  
  //Produce final state according to momentum conservation
  G4ThreeVector particleDirection1 = initialMomentum - gammaEnergy*gammaDirection1;
  particleDirection1 = particleDirection1.unit(); //normalize    

  //Update the primary particle
  if (residualPrimaryEnergy > 0.)
    {
      fParticleChange->ProposeMomentumDirection(particleDirection1);
      fParticleChange->SetProposedKineticEnergy(residualPrimaryEnergy);
    }
  else
    {
      fParticleChange->SetProposedKineticEnergy(0.);
    }

  //Now produce the photon
  G4DynamicParticle* theGamma = new G4DynamicParticle(G4Gamma::Gamma(),
						      gammaDirection1,
						      gammaEnergy);
  fvect->push_back(theGamma);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeBremsstrahlung" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << residualPrimaryEnergy/keV << " keV" << G4endl;
      G4cout << "Bremsstrahlung photon " << gammaEnergy/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (residualPrimaryEnergy+gammaEnergy)/keV 
	     << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(residualPrimaryEnergy+gammaEnergy-kineticEnergy);
      if (energyDiff > 0.05*keV)
        G4cout << "Warning from G4PenelopeBremsstrahlung: problem with energy conservation: " <<
          (residualPrimaryEnergy+gammaEnergy)/keV <<
          " keV (final) vs. " <<
          kineticEnergy/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungAngular* G4PenelopeBremsstrahlungModel::GetAngularDataForZ(G4int iZ)
{  
  if (!angularData)
    angularData = new std::map<G4int,G4PenelopeBremsstrahlungAngular*>;

  if (angularData->count(iZ)) //the material already exists
    return angularData->find(iZ)->second;

  //Otherwise create a new object, store it and return it
  G4PenelopeBremsstrahlungAngular* theAngular = new G4PenelopeBremsstrahlungAngular(iZ);
  angularData->insert(std::make_pair(iZ,theAngular));

  if (angularData->count(iZ)) //the material should exist now
    return angularData->find(iZ)->second;
  else
    {
      G4Exception("Problem in G4PenelopeBremsstrahlungModel::GetAngularDataForZ()");
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungContinuous* 
G4PenelopeBremsstrahlungModel::GetStoppingPowerData(G4int iZ,G4double energyCut,
						    const G4ParticleDefinition* 
						    theParticle)
{  
  if (!stoppingPowerData)
    stoppingPowerData = new std::map<std::pair<G4int,G4double>,G4PenelopeBremsstrahlungContinuous*>;

  std::pair<G4int,G4double> theKey = std::make_pair(iZ,energyCut);

  if (stoppingPowerData->count(theKey)) //the material already exists
    return stoppingPowerData->find(theKey)->second;

  //Otherwise create a new object, store it and return it
  G4String theParticleName = theParticle->GetParticleName();
  G4PenelopeBremsstrahlungContinuous* theContinuous = new 
    G4PenelopeBremsstrahlungContinuous(iZ,energyCut,LowEnergyLimit(),HighEnergyLimit(),theParticleName);
  stoppingPowerData->insert(std::make_pair(theKey,theContinuous));

  if (stoppingPowerData->count(theKey)) //the material should exist now
    return stoppingPowerData->find(theKey)->second;
  else
    {
      G4Exception("Problem in G4PenelopeBremsstrahlungModel::GetStoppingPowerData()");
      return 0;
    }
}
