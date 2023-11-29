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
//
// G4MicroElecInelasticModel.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//
//          - Inelastic cross-sections of low energy electrons in silicon
//	    for the simulation of heavy ion tracks with theGeant4-DNA toolkit,
//	    NSS Conf. Record 2010, pp. 80-85.
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    NIM B, vol. 288, pp. 66 - 73, 2012.
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for protons and
//	    heavy ions in Si, NIM B, vol. 287, pp. 124 - 129, 2012.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MicroElecInelasticModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4MicroElecSiStructure.hh"
#include "G4LossTableManager.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4MicroElecCrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4NistManager.hh"
#include "G4LogLogInterpolation.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecInelasticModel::G4MicroElecInelasticModel(const G4ParticleDefinition*,
                                                     const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
  nistSi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if( verboseLevel>0 )
  {
    G4cout << "MicroElec inelastic model is constructed " << G4endl;
  }
  
  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fAtomDeexcitation = nullptr;
  fParticleChangeForGamma = nullptr;
  
  // default generator
  SetAngularDistribution(new G4DeltaAngle());
  
  // Selection of computation method
  fasterCode = true; //false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MicroElecInelasticModel::~G4MicroElecInelasticModel()
{
  // Cross section
  for (auto & pos : tableData)
  {
    G4MicroElecCrossSectionDataSet* table = pos.second;
    delete table;
  }
  
  // Final state
  eVecm.clear();
  pVecm.clear();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecInelasticModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{
  
  if (verboseLevel > 3)
    G4cout << "Calling G4MicroElecInelasticModel::Initialise()" << G4endl;
  
  // Energy limits
  
  G4String fileElectron("microelec/sigma_inelastic_e_Si");
  G4String fileProton("microelec/sigma_inelastic_p_Si");
  
  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4String electron = electronDef->GetParticleName();
  G4String proton = protonDef->GetParticleName();;
  
  G4double scaleFactor = 1e-18 * cm *cm;

  const char* path = G4FindDataDir("G4LEDATA");

  // *** ELECTRON
  tableFile[electron] = fileElectron;  
  lowEnergyLimit[electron] = 16.7 * eV;
  highEnergyLimit[electron] = 100.0 * MeV;
  
  // Cross section  
  G4MicroElecCrossSectionDataSet* tableE = new G4MicroElecCrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableE->LoadData(fileElectron);
  
  tableData[electron] = tableE;
  
  // Final state
  
  std::ostringstream eFullFileName;
  
  if (fasterCode) eFullFileName << path << "/microelec/sigmadiff_cumulated_inelastic_e_Si.dat";
  else eFullFileName << path << "/microelec/sigmadiff_inelastic_e_Si.dat";
  
  std::ifstream eDiffCrossSection(eFullFileName.str().c_str());
  
  if (!eDiffCrossSection)
  {
     if (fasterCode) G4Exception("G4MicroElecInelasticModel::Initialise","em0003",
     FatalException,"Missing data file:/microelec/sigmadiff_cumulated_inelastic_e_Si.dat");

     else G4Exception("G4MicroElecInelasticModel::Initialise","em0003",
     FatalException,"Missing data file:/microelec/sigmadiff_inelastic_e_Si.dat");   
  }

  // Clear the arrays for re-initialization case (MT mode)
  // Octobre 22nd, 2014 - Melanie Raine  
  eTdummyVec.clear();
  pTdummyVec.clear();
  
  eVecm.clear();
  pVecm.clear();
  
  for (int j=0; j<6; j++)
  {
    eProbaShellMap[j].clear();
    pProbaShellMap[j].clear();

    eDiffCrossSectionData[j].clear();
    pDiffCrossSectionData[j].clear();

    eNrjTransfData[j].clear();
    pNrjTransfData[j].clear();
  }
  
  //
  eTdummyVec.push_back(0.);
  while(!eDiffCrossSection.eof())
  {
    G4double tDummy;
    G4double eDummy;
    eDiffCrossSection>>tDummy>>eDummy;
    if (tDummy != eTdummyVec.back()) eTdummyVec.push_back(tDummy);
    
    G4double tmp;
    for (int j=0; j<6; j++)
    {
      eDiffCrossSection>> tmp;
      
      eDiffCrossSectionData[j][tDummy][eDummy] = tmp;
      
      if (fasterCode)
      {
        eNrjTransfData[j][tDummy][eDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
        eProbaShellMap[j][tDummy].push_back(eDiffCrossSectionData[j][tDummy][eDummy]);
      }
      else
      {
      // SI - only if eof is not reached !
	if (!eDiffCrossSection.eof()) eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;
	eVecm[tDummy].push_back(eDummy);
      }
      
    }
  }
  //
  
  // *** PROTON   
  tableFile[proton] = fileProton;  
  lowEnergyLimit[proton] = 50. * keV;
  highEnergyLimit[proton] = 10. * GeV;
  
  // Cross section
  G4MicroElecCrossSectionDataSet* tableP = new G4MicroElecCrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableP->LoadData(fileProton);
  tableData[proton] = tableP;
  
  // Final state  
  std::ostringstream pFullFileName;
  
  if (fasterCode) pFullFileName << path << "/microelec/sigmadiff_cumulated_inelastic_p_Si.dat";
  else pFullFileName << path << "/microelec/sigmadiff_inelastic_p_Si.dat";

  std::ifstream pDiffCrossSection(pFullFileName.str().c_str());
  
  if (!pDiffCrossSection)
  {
    if (fasterCode) G4Exception("G4MicroElecInelasticModel::Initialise","em0003",
      FatalException,"Missing data file:/microelec/sigmadiff_cumulated_inelastic_p_Si.dat");

    else G4Exception("G4MicroElecInelasticModel::Initialise","em0003",
      FatalException,"Missing data file:/microelec/sigmadiff_inelastic_p_Si.dat");
  }
  
  pTdummyVec.push_back(0.);
  while(!pDiffCrossSection.eof())
  {
    G4double tDummy;
    G4double eDummy;
    pDiffCrossSection>>tDummy>>eDummy;
    if (tDummy != pTdummyVec.back()) pTdummyVec.push_back(tDummy);
    for (int j=0; j<6; j++)
    {
      pDiffCrossSection>>pDiffCrossSectionData[j][tDummy][eDummy];
      
      if (fasterCode)
      {
        pNrjTransfData[j][tDummy][pDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
        pProbaShellMap[j][tDummy].push_back(pDiffCrossSectionData[j][tDummy][eDummy]);
      }
      else
      {
	// SI - only if eof is not reached !
	if (!pDiffCrossSection.eof()) pDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;
	pVecm[tDummy].push_back(eDummy);
      }
    }
  }
  
  if (particle==electronDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[electron]);
    SetHighEnergyLimit(highEnergyLimit[electron]);
  }
  
  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }
  
  if( verboseLevel>0 )
  {
    G4cout << "MicroElec Inelastic model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / keV << " keV - "
    << HighEnergyLimit() / MeV << " MeV for "
    << particle->GetParticleName()
	   << " with mass (amu) " << particle->GetPDGMass()/proton_mass_c2
	   << " and charge " << particle->GetPDGCharge()
    << G4endl << G4endl ;
  }
  
  //
  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();
  
  if (isInitialised) { return; }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MicroElecInelasticModel::CrossSectionPerVolume(const G4Material* material,
                                                          const G4ParticleDefinition* particleDefinition,
                                                          G4double ekin,
                                                          G4double,
                                                          G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4MicroElecInelasticModel" << G4endl;
  
  G4double density = material->GetTotNbOfAtomsPerVolume();
  
  // Calculate total cross section for model 
  G4double lowLim = 0;
  G4double highLim = 0;
  G4double sigma=0;
  
  const G4String& particleName = particleDefinition->GetParticleName();
  G4String nameLocal = particleName ;
  
  G4double Zeff2 = 1.0;
  G4double Mion_c2 = particleDefinition->GetPDGMass();
  
  if (Mion_c2 > proton_mass_c2)
  {
    G4ionEffectiveCharge EffCharge ;
    G4double Zeff = EffCharge.EffectiveCharge(particleDefinition, material,ekin);
    Zeff2 = Zeff*Zeff;
    
    if (verboseLevel > 3)
      G4cout << "Before scaling : " << G4endl
      << "Particle : " << nameLocal << ", mass : " << Mion_c2/proton_mass_c2 << "*mp, charge " << Zeff
      << ", Ekin (eV) = " << ekin/eV << G4endl ;
    
    ekin *= proton_mass_c2/Mion_c2 ;
    nameLocal = "proton" ;
    
    if (verboseLevel > 3)
      G4cout << "After scaling : " << G4endl
      << "Particle : " << nameLocal  << ", Ekin (eV) = " << ekin/eV << G4endl ;
  }
  
  if (material == nistSi || material->GetBaseMaterial() == nistSi)
  {    
    auto pos1 = lowEnergyLimit.find(nameLocal);
    if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }
    
    auto pos2 = highEnergyLimit.find(nameLocal);
    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }
    
    if (ekin >= lowLim && ekin < highLim)
    {
      auto pos = tableData.find(nameLocal);      
      if (pos != tableData.end())
      {
        G4MicroElecCrossSectionDataSet* table = pos->second;
        if (table != nullptr)
        {
          sigma = table->FindValue(ekin);
        }
      }
      else
      {
        G4Exception("G4MicroElecInelasticModel::CrossSectionPerVolume","em0002",
		    FatalException,"Model not applicable to particle type.");
      }
    }
    else
    {
      if (nameLocal!="e-")
      {
        // G4cout << "Particle : " << nameLocal << ", Ekin (eV) = " << ekin/eV << G4endl;
        // G4cout << "### Warning: particle energy out of bounds! ###" << G4endl;
      }
    }
    
    if (verboseLevel > 3)
    {
      G4cout << "---> Kinetic energy (eV)=" << ekin/eV << G4endl;
      G4cout << " - Cross section per Si atom (cm^2)=" << sigma*Zeff2/cm2 << G4endl;
      G4cout << " - Cross section per Si atom (cm^-1)=" << sigma*density*Zeff2/(1./cm) << G4endl;
    }
    
  } // if (SiMaterial)
  return sigma*density*Zeff2;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecInelasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                  const G4MaterialCutsCouple* couple,
                                                  const G4DynamicParticle* particle,
                                                  G4double,
                                                  G4double)
{
  
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4MicroElecInelasticModel" << G4endl;
  
  G4double lowLim = 0;
  G4double highLim = 0;
  
  G4double ekin = particle->GetKineticEnergy();
  G4double k = ekin ;
  
  G4ParticleDefinition* PartDef = particle->GetDefinition();
  const G4String& particleName = PartDef->GetParticleName();
  G4String nameLocal2 = particleName ;
  G4double particleMass = particle->GetDefinition()->GetPDGMass();
  
  if (particleMass > proton_mass_c2)
  {
    k *= proton_mass_c2/particleMass ;
    PartDef = G4Proton::ProtonDefinition();
    nameLocal2 = "proton" ;
  }
  
  auto pos1 = lowEnergyLimit.find(nameLocal2); 
  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }
  
  auto pos2 = highEnergyLimit.find(nameLocal2);
  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }
  
  if (k >= lowLim && k < highLim)
  {
    G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
    G4double totalEnergy = ekin + particleMass;
    G4double pSquare = ekin * (totalEnergy + particleMass);
    G4double totalMomentum = std::sqrt(pSquare);
    
    G4int Shell = 0;
    
    /* if (!fasterCode)*/ Shell = RandomSelect(k,nameLocal2);
    
    // SI: The following protection is necessary to avoid infinite loops :
    //  sigmadiff_ionisation_e_born.dat has non zero partial xs at 18 eV for shell 3 (ionizationShell ==2)
    //  sigmadiff_cumulated_ionisation_e_born.dat has zero cumulated partial xs at 18 eV for shell 3 (ionizationShell ==2)
    //  this is due to the fact that the max allowed transfered energy is (18+10.79)/2=17.025 eV and only transfered energies
    //  strictly above this value have non zero partial xs in sigmadiff_ionisation_e_born.dat (starting at trans = 17.12 eV)    
    
    G4double bindingEnergy = SiStructure.Energy(Shell);
    
    if (verboseLevel > 3)
    {
      G4cout << "---> Kinetic energy (eV)=" << k/eV << G4endl ;
      G4cout << "Shell: " << Shell << ", energy: " << bindingEnergy/eV << G4endl;
    }
    
    // sample deexcitation
    std::size_t secNumberInit = 0;  // need to know at a certain point the energy of secondaries
    std::size_t secNumberFinal = 0; // So I'll make the difference and then sum the energies
    
    //SI: additional protection if tcs interpolation method is modified
    if (k<bindingEnergy) return;
    
    G4int Z = 14;
    
    if(fAtomDeexcitation && Shell > 2) {
      
      G4AtomicShellEnumerator as = fKShell;
      
      if (Shell == 4)
      {
        as = G4AtomicShellEnumerator(1);
      }
      else if (Shell == 3)
      {
        as = G4AtomicShellEnumerator(3);
      }
      
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      secNumberInit = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
      secNumberFinal = fvect->size();
    }
    
    G4double secondaryKinetic=-1000*eV;

    if (!fasterCode)
    {
      secondaryKinetic = RandomizeEjectedElectronEnergy(PartDef,k,Shell);
    }
    else
    {
      secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulatedDcs(PartDef,k,Shell);
    }
    
    if (verboseLevel > 3)
    {
      G4cout << "Ionisation process" << G4endl;
      G4cout << "Shell: " << Shell << " Kin. energy (eV)=" << k/eV
      << " Sec. energy (eV)=" << secondaryKinetic/eV << G4endl;
    }
    
    G4ThreeVector deltaDirection =
    GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic,
                                                      Z, Shell,
                                                      couple->GetMaterial());
    
    if (particle->GetDefinition() == G4Electron::ElectronDefinition())
    {
      G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));    
      G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
      G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
      G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;
      
      G4ThreeVector direction;
      direction.set(finalPx,finalPy,finalPz);
      
      fParticleChangeForGamma->ProposeMomentumDirection(direction.unit()) ;
    }
    else 
      fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection) ;
    
    // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
    G4double deexSecEnergy = 0;
    for (std::size_t j=secNumberInit; j < secNumberFinal; ++j) {
      deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();}
    
    fParticleChangeForGamma->SetProposedKineticEnergy(ekin-bindingEnergy-secondaryKinetic);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy-deexSecEnergy);
    
    if (secondaryKinetic>0)
    {
      G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
      fvect->push_back(dp);
    }     
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition,
                                                                   G4double k, G4int shell)
{
  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    G4double maximumEnergyTransfer=0.;
    if ((k+SiStructure.Energy(shell))/2. > k) maximumEnergyTransfer=k;
    else maximumEnergyTransfer = (k+SiStructure.Energy(shell))/2.;
    
    G4double crossSectionMaximum = 0.;
    
    G4double minEnergy = SiStructure.Energy(shell);
    G4double maxEnergy = maximumEnergyTransfer;
    G4int nEnergySteps = 100;
    
    G4double value(minEnergy);
    G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
    G4int step(nEnergySteps);
    while (step>0)
    {
      step--;
      G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
      if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
      value*=stpEnergy;
    }
    
    G4double secondaryElectronKineticEnergy=0.;
    do
    {
      secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer-SiStructure.Energy(shell));
    } while(G4UniformRand()*crossSectionMaximum >
            DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+SiStructure.Energy(shell))/eV,shell));
    
    return secondaryElectronKineticEnergy;    
  }
  
  if (particleDefinition == G4Proton::ProtonDefinition())
  {
    G4double maximumEnergyTransfer = 4.* (electron_mass_c2 / proton_mass_c2) * k;
    G4double crossSectionMaximum = 0.;
    
    G4double minEnergy = SiStructure.Energy(shell);
    G4double maxEnergy = maximumEnergyTransfer;
    G4int nEnergySteps = 100;
    
    G4double value(minEnergy);
    G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
    G4int step(nEnergySteps);
    while (step>0)
    {
      step--;
      G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
      if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
      value*=stpEnergy;
    }
    
    G4double secondaryElectronKineticEnergy = 0.;
    do
    {
      secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer-SiStructure.Energy(shell));
      
    } while(G4UniformRand()*crossSectionMaximum >
            DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+SiStructure.Energy(shell))/eV,shell));
    return secondaryElectronKineticEnergy;
  }
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The following section is not used anymore but is kept for memory
// GetAngularDistribution()->SampleDirectionForShell is used instead

/*void G4MicroElecInelasticModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition,
 G4double k,
 G4double secKinetic,
 G4double & cosTheta,
 G4double & phi )
 {
 if (particleDefinition == G4Electron::ElectronDefinition())
 {
 phi = twopi * G4UniformRand();
 G4double sin2O = (1.-secKinetic/k) / (1.+secKinetic/(2.*electron_mass_c2));
 cosTheta = std::sqrt(1.-sin2O);
 }
 
 if (particleDefinition == G4Proton::ProtonDefinition())
 {
 G4double maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
 phi = twopi * G4UniformRand();
 cosTheta = std::sqrt(secKinetic / maxSecKinetic);
 }
 
 else
 {
 G4double maxSecKinetic = 4.* (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
 phi = twopi * G4UniformRand();
 cosTheta = std::sqrt(secKinetic / maxSecKinetic);
 }
 }
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition,
                                                           G4double k,
                                                           G4double energyTransfer,
                                                           G4int LevelIndex)
{
  G4double sigma = 0.;
  
  if (energyTransfer >= SiStructure.Energy(LevelIndex))
  {
    G4double valueT1 = 0;
    G4double valueT2 = 0;
    G4double valueE21 = 0;
    G4double valueE22 = 0;
    G4double valueE12 = 0;
    G4double valueE11 = 0;
    
    G4double xs11 = 0;
    G4double xs12 = 0;
    G4double xs21 = 0;
    G4double xs22 = 0;
    
    if (particleDefinition == G4Electron::ElectronDefinition())
    {
      // k should be in eV and energy transfer eV also      
      auto t2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);
      auto t1 = t2-1;
      // The following condition avoids situations where energyTransfer >last vector element
      if (energyTransfer <= eVecm[(*t1)].back() && energyTransfer <= eVecm[(*t2)].back() )
      {
        auto e12 = std::upper_bound(eVecm[(*t1)].begin(),eVecm[(*t1)].end(), energyTransfer);
        auto e11 = e12-1;
        auto e22 = std::upper_bound(eVecm[(*t2)].begin(),eVecm[(*t2)].end(), energyTransfer);
        auto e21 = e22-1;
        
        valueT1  =*t1;
        valueT2  =*t2;
        valueE21 =*e21;
        valueE22 =*e22;
        valueE12 =*e12;
        valueE11 =*e11;
        
        xs11 = eDiffCrossSectionData[LevelIndex][valueT1][valueE11];
        xs12 = eDiffCrossSectionData[LevelIndex][valueT1][valueE12];
        xs21 = eDiffCrossSectionData[LevelIndex][valueT2][valueE21];
        xs22 = eDiffCrossSectionData[LevelIndex][valueT2][valueE22];
      }      
    }
    
    if (particleDefinition == G4Proton::ProtonDefinition())
    {
      // k should be in eV and energy transfer eV also
      auto t2 = std::upper_bound(pTdummyVec.begin(),pTdummyVec.end(), k);
      auto t1 = t2-1;
      if (energyTransfer <= pVecm[(*t1)].back() && energyTransfer <= pVecm[(*t2)].back() )
      {
        auto e12 = std::upper_bound(pVecm[(*t1)].begin(),pVecm[(*t1)].end(), energyTransfer);
        auto e11 = e12-1;
        
        auto e22 = std::upper_bound(pVecm[(*t2)].begin(),pVecm[(*t2)].end(), energyTransfer);
        auto e21 = e22-1;
        
        valueT1  =*t1;
        valueT2  =*t2;
        valueE21 =*e21;
        valueE22 =*e22;
        valueE12 =*e12;
        valueE11 =*e11;
        
        xs11 = pDiffCrossSectionData[LevelIndex][valueT1][valueE11];
        xs12 = pDiffCrossSectionData[LevelIndex][valueT1][valueE12];
        xs21 = pDiffCrossSectionData[LevelIndex][valueT2][valueE21];
        xs22 = pDiffCrossSectionData[LevelIndex][valueT2][valueE22];
      }
    }
    
    sigma = QuadInterpolator(     valueE11, valueE12,
				  valueE21, valueE22,
				  xs11, xs12,
				  xs21, xs22,
				  valueT1, valueT2,
				  k, energyTransfer);
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::Interpolate(G4double e1,
                                                G4double e2,
                                                G4double e,
                                                G4double xs1,
                                                G4double xs2)
{
  G4double value = 0.;

  // Log-log interpolation by default
  if (e1 != 0 && e2 != 0 && (std::log10(e2) - std::log10(e1)) != 0
      && !fasterCode)
  {  
    G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
    G4double b = std::log10(xs2) - a*std::log10(e2);
    G4double sigma = a*std::log10(e) + b;
    value = (std::pow(10.,sigma));
    
  }
  
  // Switch to log-lin interpolation for faster code
  if ((e2 - e1) != 0 && xs1 != 0 && xs2 != 0 && fasterCode)
  {
    G4double d1 = std::log10(xs1);
    G4double d2 = std::log10(xs2);
    value = std::pow(10., (d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
  }
  
  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero
  if ((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0)) // && fasterCode)
  {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::QuadInterpolator(G4double e11, G4double e12,
                                                     G4double e21, G4double e22,
                                                     G4double xs11, G4double xs12,
                                                     G4double xs21, G4double xs22,
                                                     G4double t1, G4double t2,
                                                     G4double t, G4double e)
{
  G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
  G4double value = Interpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4MicroElecInelasticModel::RandomSelect(G4double k, const G4String& particle )
{
  G4int level = 0;
 
  auto pos = tableData.find(particle); 
  if (pos != tableData.cend())
  {
    G4MicroElecCrossSectionDataSet* table = pos->second;
    
    if (table != nullptr)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
      const G4int n = (G4int)table->NumberOfComponents();
      G4int i(n);
      G4double value = 0.;
      
      while (i>0)
      {
        --i;
        valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
        value += valuesBuffer[i];
      }
      
      value *= G4UniformRand();
      
      i = n;
      
      while (i > 0)
      {
        --i;
        
        if (valuesBuffer[i] > value)
        {
          delete[] valuesBuffer;
          return i;
        }
        value -= valuesBuffer[i];
      }
      
      if (valuesBuffer) delete[] valuesBuffer;
      
    }
  }
  else
  {
    G4Exception("G4MicroElecInelasticModel::RandomSelect","em0002",FatalException,"Model not applicable to particle type.");
  }
  
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::RandomizeEjectedElectronEnergyFromCumulatedDcs(G4ParticleDefinition* particleDefinition,
                                                                                   G4double k,
                                                                                   G4int shell)
{

  G4double secondaryElectronKineticEnergy = 0.;
  G4double random = G4UniformRand();
  secondaryElectronKineticEnergy = TransferedEnergy(particleDefinition,
                                                    k / eV,
                                                    shell,
                                                    random) * eV
      - SiStructure.Energy(shell);
  
  if (secondaryElectronKineticEnergy < 0.)
    return 0.;
  //
  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecInelasticModel::TransferedEnergy(G4ParticleDefinition* particleDefinition,
                                                     G4double k,
                                                     G4int ionizationLevelIndex,
                                                     G4double random)
{
  G4double nrj = 0.;

  G4double valueK1 = 0;
  G4double valueK2 = 0;
  G4double valuePROB21 = 0;
  G4double valuePROB22 = 0;
  G4double valuePROB12 = 0;
  G4double valuePROB11 = 0;

  G4double nrjTransf11 = 0;
  G4double nrjTransf12 = 0;
  G4double nrjTransf21 = 0;
  G4double nrjTransf22 = 0;
  
  G4double maximumEnergyTransfer1 = 0;  
  G4double maximumEnergyTransfer2 = 0;  
  G4double maximumEnergyTransferP = 4.* (electron_mass_c2 / proton_mass_c2) * k;
  G4double bindingEnergy = SiStructure.Energy(ionizationLevelIndex)*1e6;

  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    // k should be in eV
    auto k2 = std::upper_bound(eTdummyVec.begin(),
			       eTdummyVec.end(),
			       k);
    auto k1 = k2 - 1;

    /*
     G4cout << "----> k=" << k
     << " " << *k1
     << " " << *k2
     << " " << random
     << " " << ionizationLevelIndex
     << " " << eProbaShellMap[ionizationLevelIndex][(*k1)].back()
     << " " << eProbaShellMap[ionizationLevelIndex][(*k2)].back()
     << G4endl;
     */

    // SI : the following condition avoids situations where random >last vector element
    if (random <= eProbaShellMap[ionizationLevelIndex][(*k1)].back()
        && random <= eProbaShellMap[ionizationLevelIndex][(*k2)].back())
    {
      auto prob12 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k1)].end(),
                           random);
      auto prob11 = prob12 - 1;
      auto prob22 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);
      auto prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;
      valuePROB12 = *prob12;
      valuePROB11 = *prob11;
      
      // The following condition avoid getting transfered energy < binding energy and forces cumxs = 1 for maximum energy transfer.
      if(valuePROB11 == 0) nrjTransf11 = bindingEnergy; 
      else nrjTransf11 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
      if(valuePROB12 == 1) 
      {	
	if ((valueK1+bindingEnergy)/2. > valueK1) maximumEnergyTransfer1=valueK1;
    	else maximumEnergyTransfer1 = (valueK1+bindingEnergy)/2.;
	
	nrjTransf12 = maximumEnergyTransfer1;
      }
      else nrjTransf12 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];

      if(valuePROB21 == 0) nrjTransf21 = bindingEnergy;
      else nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      if(valuePROB22 == 1) 
      {	
	if ((valueK2+bindingEnergy)/2. > valueK2) maximumEnergyTransfer2=valueK2;
    	else maximumEnergyTransfer2 = (valueK2+bindingEnergy)/2.;
	
	nrjTransf22 = maximumEnergyTransfer2;
      }
      else
	nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];
    }
    // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)
    if (random > eProbaShellMap[ionizationLevelIndex][(*k1)].back())
    {
      auto prob22 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      auto prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;

      nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      G4double interpolatedvalue2 = Interpolate(valuePROB21,
                                                valuePROB22,
                                                random,
                                                nrjTransf21,
                                                nrjTransf22);

      // zeros are explicitly set
      G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);     
      return value;
    }
  }
  //
  else if (particleDefinition == G4Proton::ProtonDefinition())
  {
    // k should be in eV
    auto k2 = std::upper_bound(pTdummyVec.begin(),
			       pTdummyVec.end(),
			       k);

    auto k1 = k2 - 1;

    /*
     G4cout << "----> k=" << k
     << " " << *k1
     << " " << *k2
     << " " << random
     << " " << ionizationLevelIndex
     << " " << pProbaShellMap[ionizationLevelIndex][(*k1)].back()
     << " " << pProbaShellMap[ionizationLevelIndex][(*k2)].back()
     << G4endl;
     */

    // SI : the following condition avoids situations where random > last vector element,
    //      for eg. when the last element is zero
    if (random <= pProbaShellMap[ionizationLevelIndex][(*k1)].back()
        && random <= pProbaShellMap[ionizationLevelIndex][(*k2)].back())
    {
      auto prob12 =
          std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
                           pProbaShellMap[ionizationLevelIndex][(*k1)].end(),
                           random);

      auto prob11 = prob12 - 1;

      auto prob22 =
	std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
			 pProbaShellMap[ionizationLevelIndex][(*k2)].end(),
			 random);
      
      auto prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;
      valuePROB12 = *prob12;
      valuePROB11 = *prob11;

      // The following condition avoid getting transfered energy < binding energy and forces cumxs = 1 for maximum energy transfer.
      if(valuePROB11 == 0) nrjTransf11 = bindingEnergy; 
      else nrjTransf11 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
      if(valuePROB12 == 1) nrjTransf12 = maximumEnergyTransferP;
      else nrjTransf12 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];
      if(valuePROB21 == 0) nrjTransf21 = bindingEnergy;
      else nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      if(valuePROB22 == 1) nrjTransf22 = maximumEnergyTransferP;
      else nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

    }

    // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)
    if (random > pProbaShellMap[ionizationLevelIndex][(*k1)].back())
    {
      auto prob22 =
          std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           pProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      auto prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;

      nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      G4double interpolatedvalue2 = Interpolate(valuePROB21,
                                                valuePROB22,
                                                random,
                                                nrjTransf21,
                                                nrjTransf22);

      // zeros are explicitly set
      G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
     
      return value;
    }
  }
  // End electron and proton cases

  G4double nrjTransfProduct = nrjTransf11 * nrjTransf12 * nrjTransf21
      * nrjTransf22;

  if (nrjTransfProduct != 0.)
  {
    nrj = QuadInterpolator(valuePROB11,
                           valuePROB12,
                           valuePROB21,
                           valuePROB22,
                           nrjTransf11,
                           nrjTransf12,
                           nrjTransf21,
                           nrjTransf22,
                           valueK1,
                           valueK2,
                           k,
                           random);
  }

  return nrj;
}


