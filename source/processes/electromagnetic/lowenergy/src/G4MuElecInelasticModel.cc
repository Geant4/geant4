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
// G4MuElecInelasticModel.cc, 2011/08/29 A.Valentin, M. Raine
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

#include "G4MuElecInelasticModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4ionEffectiveCharge.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuElecInelasticModel::G4MuElecInelasticModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),fAtomDeexcitation(0),isInitialised(false)
{
  
   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The name of the class G4MuElecInelasticModel is changed to G4MicroElecInelasticModel. " << G4endl;
   G4cout << "   The obsolete class will be REMOVED with the next release of Geant4. " << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
   
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
    G4cout << "MuElec inelastic model is constructed " << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fParticleChangeForGamma = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuElecInelasticModel::~G4MuElecInelasticModel()
{  
  // Cross section
  
  std::map< G4String,G4MuElecCrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4MuElecCrossSectionDataSet* table = pos->second;
    delete table;
  }
  
  // Final state
  
  eVecm.clear();
  pVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuElecInelasticModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4MuElecInelasticModel::Initialise()" << G4endl;

  // Energy limits

  G4String fileElectron("microelec/sigma_inelastic_e_Si");
  G4String fileProton("microelec/sigma_inelastic_p_Si");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

  G4String electron;
  G4String proton;
  
  G4double scaleFactor = 1e-18 * cm *cm;

  char *path = getenv("G4LEDATA");

  // *** ELECTRON
    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    lowEnergyLimit[electron] = 16.7 * eV; 
    highEnergyLimit[electron] = 100.0 * MeV;

    // Cross section
    
    G4MuElecCrossSectionDataSet* tableE = new G4MuElecCrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
      
    tableData[electron] = tableE;
    
    // Final state
    
    std::ostringstream eFullFileName;
    eFullFileName << path << "/microelec/sigmadiff_inelastic_e_Si.dat";
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

    if (!eDiffCrossSection)
    { 
            G4Exception("G4MuElecInelasticModel::Initialise","em0003",FatalException,"Missing data file:/microelec/sigmadiff_inelastic_e_Si.dat");
    }
      
    eTdummyVec.push_back(0.);
    while(!eDiffCrossSection.eof())
    {
      double tDummy;
      double eDummy;
      eDiffCrossSection>>tDummy>>eDummy;
      if (tDummy != eTdummyVec.back()) eTdummyVec.push_back(tDummy);
      for (int j=0; j<6; j++)
      {
        eDiffCrossSection>>eDiffCrossSectionData[j][tDummy][eDummy];

        // SI - only if eof is not reached !
        if (!eDiffCrossSection.eof()) eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

        eVecm[tDummy].push_back(eDummy);

      }
    }
    //

  // *** PROTON

    proton = protonDef->GetParticleName();

    tableFile[proton] = fileProton;

    lowEnergyLimit[proton] = 50. * keV;
    highEnergyLimit[proton] = 10. * GeV;

    // Cross section
    
    G4MuElecCrossSectionDataSet* tableP = new G4MuElecCrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableP->LoadData(fileProton);
      
    tableData[proton] = tableP;
    
    // Final state

    std::ostringstream pFullFileName;
    pFullFileName << path << "/microelec/sigmadiff_inelastic_p_Si.dat";
    std::ifstream pDiffCrossSection(pFullFileName.str().c_str());
    
    if (!pDiffCrossSection)
    { 
            G4Exception("G4MuElecInelasticModel::Initialise","em0003",FatalException,"Missing data file:/microelec/sigmadiff_inelastic_p_Si.dat");
    }
      
    pTdummyVec.push_back(0.);
    while(!pDiffCrossSection.eof())
    {
      double tDummy;
      double eDummy;
      pDiffCrossSection>>tDummy>>eDummy;
      if (tDummy != pTdummyVec.back()) pTdummyVec.push_back(tDummy);
      for (int j=0; j<6; j++)
      {
        pDiffCrossSection>>pDiffCrossSectionData[j][tDummy][eDummy];

        // SI - only if eof is not reached !
        if (!pDiffCrossSection.eof()) pDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

        pVecm[tDummy].push_back(eDummy); 
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
    G4cout << "MuElec Inelastic model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
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

G4double G4MuElecInelasticModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* particleDefinition,
					   G4double ekin,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4MuElecInelasticModel" << G4endl;

  G4double density = material->GetTotNbOfAtomsPerVolume();

 /* if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != G4Electron::ElectronDefinition()
      &&
      particleDefinition != G4GenericIon::GenericIonDefinition()
     )
   	    
    return 0;*/
  
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

    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(nameLocal);
    if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }
  
    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(nameLocal);
    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

    if (ekin >= lowLim && ekin < highLim)
    {
      std::map< G4String,G4MuElecCrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(nameLocal);
	
      if (pos != tableData.end())
      {
        G4MuElecCrossSectionDataSet* table = pos->second;
        if (table != 0)
        {
	  sigma = table->FindValue(ekin);
        }
      }
      else
      {
	G4Exception("G4MuElecInelasticModel::CrossSectionPerVolume","em0002",FatalException,"Model not applicable to particle type.");
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

void G4MuElecInelasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* /*couple*/,
					      const G4DynamicParticle* particle,
					      G4double,
					      G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4MuElecInelasticModel" << G4endl;

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

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(nameLocal2);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(nameLocal2);

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

    G4int Shell = RandomSelect(k,nameLocal2);
    G4double bindingEnergy = SiStructure.Energy(Shell);
    if (verboseLevel > 3)
    {
	G4cout << "---> Kinetic energy (eV)=" << k/eV << G4endl ;
	G4cout << "Shell: " << Shell << ", energy: " << bindingEnergy/eV << G4endl;
    }

   // sample deexcitation

    G4int secNumberInit = 0;  // need to know at a certain point the energy of secondaries   
    G4int secNumberFinal = 0; // So I'll make the difference and then sum the energies

    if(fAtomDeexcitation && Shell > 2) {
      G4int Z = 14;
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

    G4double secondaryKinetic = RandomizeEjectedElectronEnergy(PartDef,k,Shell);

    if (verboseLevel > 3)
    {
	G4cout << "Ionisation process" << G4endl;
	G4cout << "Shell: " << Shell << " Kin. energy (eV)=" << k/eV 
	<< " Sec. energy (eV)=" << secondaryKinetic/eV << G4endl;
    }

    G4double cosTheta = 0.;
    G4double phi = 0.; 
    RandomizeEjectedElectronDirection(PartDef, k, secondaryKinetic, cosTheta, phi);

    G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
    G4double dirX = sinTheta*std::cos(phi);
    G4double dirY = sinTheta*std::sin(phi);
    G4double dirZ = cosTheta;
    G4ThreeVector deltaDirection(dirX,dirY,dirZ);
    deltaDirection.rotateUz(primaryDirection);

    //if (particle->GetDefinition() == G4Electron::ElectronDefinition())
    //{
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
    //}
    //else fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection) ;

    // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
    G4double deexSecEnergy = 0;
    for (G4int j=secNumberInit; j < secNumberFinal; j++) {
      deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();} 

    fParticleChangeForGamma->SetProposedKineticEnergy(ekin-bindingEnergy-secondaryKinetic);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy-deexSecEnergy);

    G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
    fvect->push_back(dp);

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuElecInelasticModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
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

    } while(G4UniformRand()*crossSectionMaximum >= 
	      DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+SiStructure.Energy(shell))/eV,shell));
    return secondaryElectronKineticEnergy;
  }
 
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuElecInelasticModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double G4MuElecInelasticModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition, 
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

      std::vector<double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);
      std::vector<double>::iterator t1 = t2-1;
      // SI : the following condition avoids situations where energyTransfer >last vector element
      if (energyTransfer <= eVecm[(*t1)].back() && energyTransfer <= eVecm[(*t2)].back() )
      {
        std::vector<double>::iterator e12 = std::upper_bound(eVecm[(*t1)].begin(),eVecm[(*t1)].end(), energyTransfer);
        std::vector<double>::iterator e11 = e12-1;

        std::vector<double>::iterator e22 = std::upper_bound(eVecm[(*t2)].begin(),eVecm[(*t2)].end(), energyTransfer);
        std::vector<double>::iterator e21 = e22-1;

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
      std::vector<double>::iterator t2 = std::upper_bound(pTdummyVec.begin(),pTdummyVec.end(), k);
      std::vector<double>::iterator t1 = t2-1;
      if (energyTransfer <= pVecm[(*t1)].back() && energyTransfer <= pVecm[(*t2)].back() )
      {
        std::vector<double>::iterator e12 = std::upper_bound(pVecm[(*t1)].begin(),pVecm[(*t1)].end(), energyTransfer);
		std::vector<double>::iterator e11 = e12-1;

        std::vector<double>::iterator e22 = std::upper_bound(pVecm[(*t2)].begin(),pVecm[(*t2)].end(), energyTransfer);
        std::vector<double>::iterator e21 = e22-1;
 
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

   G4double xsProduct = xs11 * xs12 * xs21 * xs22;
   if (xsProduct != 0.)
   {
     sigma = QuadInterpolator(     valueE11, valueE12, 
				   valueE21, valueE22, 
				   xs11, xs12, 
				   xs21, xs22, 
				   valueT1, valueT2, 
				   k, energyTransfer);
   }

 }
  
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuElecInelasticModel::LogLogInterpolate(G4double e1, 
						       G4double e2, 
						       G4double e, 
						       G4double xs1, 
						       G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuElecInelasticModel::QuadInterpolator(G4double e11, G4double e12, 
						      G4double e21, G4double e22, 
						      G4double xs11, G4double xs12, 
						      G4double xs21, G4double xs22, 
						      G4double t1, G4double t2, 
						      G4double t, G4double e)
{
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4MuElecInelasticModel::RandomSelect(G4double k, const G4String& particle )
{   
  G4int level = 0;

  std::map< G4String,G4MuElecCrossSectionDataSet*,std::less<G4String> >::iterator pos;
  pos = tableData.find(particle);

  if (pos != tableData.end())
  {
    G4MuElecCrossSectionDataSet* table = pos->second;

    if (table != 0)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
      const size_t n(table->NumberOfComponents());
      size_t i(n);
      G4double value = 0.;
	    
      while (i>0)
      { 
	i--;
	valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
	value += valuesBuffer[i];
      }
	    
      value *= G4UniformRand();
    
      i = n;
	    
      while (i > 0)
      {
	i--;

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
    G4Exception("G4MuElecInelasticModel::RandomSelect","em0002",FatalException,"Model not applicable to particle type.");
  }
      
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


