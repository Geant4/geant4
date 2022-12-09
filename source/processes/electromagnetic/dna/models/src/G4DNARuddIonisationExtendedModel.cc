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
// Modified by Z. Francis, S. Incerti to handle HZE 
// && inverse rudd function sampling 26-10-2010

#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"

static G4Pow * gpow = G4Pow::GetInstance();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::G4DNARuddIonisationExtendedModel(const G4ParticleDefinition*,
                                                                   const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
    slaterEffectiveCharge[0]=0.;
    slaterEffectiveCharge[1]=0.;
    slaterEffectiveCharge[2]=0.;
    sCoefficient[0]=0.;
    sCoefficient[1]=0.;
    sCoefficient[2]=0.;

    lowEnergyLimitForA[1] = 0 * eV;
    lowEnergyLimitForA[2] = 0 * eV;
    lowEnergyLimitForA[3] = 0 * eV;
    lowEnergyLimitOfModelForA[1] = 100 * eV;
    lowEnergyLimitOfModelForA[4] = 1 * keV;
    lowEnergyLimitOfModelForA[5] = 0.5 * MeV; // For A = 3 or above, limit is MeV/uma
    killBelowEnergyForA[1] = lowEnergyLimitOfModelForA[1];
    killBelowEnergyForA[4] = lowEnergyLimitOfModelForA[4];
    killBelowEnergyForA[5] = lowEnergyLimitOfModelForA[5];

    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 )
    {
        G4cout << "Rudd ionisation model is constructed " << G4endl;
    }

    // Define default angular generator
    SetAngularDistribution(new G4DNARuddAngle());

    // Mark this model as "applicable" for atomic deexcitation
    SetDeexcitationFlag(true);

    // Selection of stationary mode

    statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::~G4DNARuddIonisationExtendedModel()
{  
  if(isIon) {
    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
        G4DNACrossSectionDataSet* table = pos->second;
        delete table;
    }
  } else {
    delete mainTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::Initialise(const G4ParticleDefinition* particle,
                                                  const G4DataVector& /*cuts*/)
{
    if (isInitialised) { return; }
    if (verboseLevel > 3)
        G4cout << "Calling G4DNARuddIonisationExtendedModel::Initialise()" << G4endl;

    // Energy limits
    G4String fileProton("dna/sigma_ionisation_p_rudd");
    G4String fileHydrogen("dna/sigma_ionisation_h_rudd");
    G4String fileAlphaPlusPlus("dna/sigma_ionisation_alphaplusplus_rudd");
    G4String fileAlphaPlus("dna/sigma_ionisation_alphaplus_rudd");
    G4String fileHelium("dna/sigma_ionisation_he_rudd");
    G4String fileCarbon("dna/sigma_ionisation_c_rudd");
    G4String fileNitrogen("dna/sigma_ionisation_n_rudd");
    G4String fileOxygen("dna/sigma_ionisation_o_rudd");
    G4String fileSilicon("dna/sigma_ionisation_si_rudd");
    G4String fileIron("dna/sigma_ionisation_fe_rudd");

    G4String pname = particle->GetParticleName();

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
    protonDef = G4Proton::ProtonDefinition();
    hydrogenDef = instance->GetIon("hydrogen");
    alphaPlusPlusDef = G4Alpha::Alpha();
    alphaPlusDef = instance->GetIon("alpha+");
    heliumDef = instance->GetIon("helium");

    carbonDef = instance->GetIon("carbon");
    nitrogenDef = instance->GetIon("nitrogen");
    oxygenDef = instance->GetIon("oxygen");
    siliconDef = instance->GetIon("silicon");
    ironDef = instance->GetIon("iron");

    G4String carbon;
    G4String nitrogen;
    G4String oxygen;
    G4String silicon;
    G4String iron;

    G4double scaleFactor = 1 * m*m;
    massC12 = carbonDef->GetPDGMass();

    // LIMITS AND DATA

    // **********************************************************************************************

    if(pname == "proton") {
      localMinEnergy = lowEnergyLimitForA[1];

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileProton);

    // **********************************************************************************************
    
    } else if(pname == "hydrogen") {

      localMinEnergy = lowEnergyLimitForA[1];

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileHydrogen);

    // **********************************************************************************************
    
    } else if(pname == "alpha") {

      localMinEnergy = lowEnergyLimitForA[4];

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileAlphaPlusPlus);

    // **********************************************************************************************
    
    } else if(pname == "alpha+") {

      localMinEnergy = lowEnergyLimitForA[4];

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileAlphaPlus);

    // **********************************************************************************************

    } else if(pname == "helium") {

      localMinEnergy = lowEnergyLimitForA[4];

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileHelium);

    // **********************************************************************************************

    } else if(pname == "GenericIon") {
    
      isIon = true;
      carbon = carbonDef->GetParticleName();
      localMinEnergy = lowEnergyLimitForA[5]*massC12/proton_mass_c2;

      // Cross section
      mainTable = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
      mainTable->LoadData(fileCarbon);

      tableData[carbon] = mainTable;

    // **********************************************************************************************
    
      oxygen = oxygenDef->GetParticleName();
      tableFile[oxygen] = fileOxygen;

      // Cross section
      G4DNACrossSectionDataSet* tableOxygen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                           eV, scaleFactor );
      tableOxygen->LoadData(fileOxygen);
      tableData[oxygen] = tableOxygen;

    // **********************************************************************************************
    
      nitrogen = nitrogenDef->GetParticleName();
      tableFile[nitrogen] = fileNitrogen;

      // Cross section
      G4DNACrossSectionDataSet* tableNitrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                             eV, scaleFactor );
      tableNitrogen->LoadData(fileNitrogen);
      tableData[nitrogen] = tableNitrogen;

    // **********************************************************************************************

      silicon = siliconDef->GetParticleName();
      tableFile[silicon] = fileSilicon;
    
      // Cross section
      G4DNACrossSectionDataSet* tableSilicon = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                            eV, scaleFactor );
      tableSilicon->LoadData(fileSilicon);
      tableData[silicon] = tableSilicon;
     
    // **********************************************************************************************
    
      iron = ironDef->GetParticleName();
      tableFile[iron] = fileIron;

      // Cross section

      G4DNACrossSectionDataSet* tableIron = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV, scaleFactor );
      tableIron->LoadData(fileIron);
      tableData[iron] = tableIron;
    }
    // **********************************************************************************************

    if( verboseLevel>0 )
    {
      G4cout << "Rudd ionisation model is initialized " << G4endl
	     << "Energy range for model: "
	     << LowEnergyLimit() / keV << " keV - "
	     << HighEnergyLimit() / keV << " keV for "
	     << pname
	     << " internal low energy limit E(keV)=" << localMinEnergy / keV
	     << G4endl;
    }

    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

    //
    fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

    fParticleChangeForGamma = GetParticleChangeForGamma();
    isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARuddIonisationExtendedModel::CrossSectionPerVolume(const G4Material* material,
                                                                 const G4ParticleDefinition* partDef,
                                                                 G4double k,
                                                                 G4double,
                                                                 G4double)
{
    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNARuddIonisationExtendedModel" << G4endl;

    currParticle = GetDNAIonParticleDefinition(partDef);
    currentScaledEnergy = k;
    G4double e = k;
    G4double q2 = 1.0;
    currentTable = mainTable;

    if (isIon){
      if (currParticle == nullptr) {//not DNA particle
        currentScaledEnergy *= massC12/partDef->GetPDGMass();
        G4double q = partDef->GetPDGCharge()/(eplus*6);
        q2 *= q*q;
        e = currentScaledEnergy;
        currParticle = carbonDef;
      }
      G4String pname = currParticle->GetParticleName();
      auto goodTable = tableData.find(pname);
      currentTable = goodTable->second;
    }
    // below low the ion should be stopped
    if (currentScaledEnergy < localMinEnergy) { return DBL_MAX; }

    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];
    G4double sigma = 0.0;
    if (nullptr != currentTable) {
      sigma = currentTable->FindValue(e)*q2;
    } else {
      G4cout << "G4DNARuddIonisationExtendedModel - no data table for " 
	     << partDef->GetParticleName() << G4endl;
      G4Exception("G4DNARuddIonisationExtendedModel::CrossSectionPerVolume(...)","em0002",
		  FatalException,"Data table is not available for the model.");
    }
    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "G4DNARuddIonisationExtendedModel - XS INFO START" << G4endl;
      G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << partDef->GetParticleName() << G4endl;
      G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
      G4cout << "G4DNARuddIonisationExtendedModel - XS INFO END" << G4endl;
    }
    return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                         const G4MaterialCutsCouple* couple,
                                                         const G4DynamicParticle* particle,
                                                         G4double,
                                                         G4double)
{
    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNARuddIonisationExtendedModel" << G4endl;

    // stop ion with energy below low energy limit
    G4double k = particle->GetKineticEnergy();
    if (currentScaledEnergy < localMinEnergy) {
      fParticleChangeForGamma->SetProposedKineticEnergy(0.);
      fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
      return;
    }

    // sampling of final state
    G4ParticleDefinition* definition = particle->GetDefinition();
    const G4ThreeVector& primaryDirection = particle->GetMomentumDirection();

    G4int ionizationShell = RandomSelect(currentScaledEnergy);

    // sample deexcitation
    // here we assume that H_{2}O electronic levels are the same as Oxygen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    G4double bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

    //SI: additional protection if tcs interpolation method is modified
    if (k < bindingEnergy) return;
    //

    G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);

    // is ionisation possible?
    G4double scatteredEnergy = k - bindingEnergy - secondaryKinetic;
    if(scatteredEnergy < 0.0) { return; }

        G4int Z = 8;
	
	G4ThreeVector deltaDirection = 
	  GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic, 
							    Z, ionizationShell,
							    couple->GetMaterial());

        G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
        fvect->push_back(dp);

        fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
		
        size_t secNumberInit = 0;// need to know at a certain point the energy of secondaries
        size_t secNumberFinal = 0;// So I'll make the diference and then sum the energies

        // SI: only atomic deexcitation from K shell is considered
        if(fAtomDeexcitation != nullptr && ionizationShell == 4)
        {
          const G4AtomicShell* shell 
            = fAtomDeexcitation->GetAtomicShell(Z, G4AtomicShellEnumerator(0));
          secNumberInit = fvect->size();
          fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
          secNumberFinal = fvect->size();

          if(secNumberFinal > secNumberInit) 
          {
	    for (size_t i=secNumberInit; i<secNumberFinal; ++i) 
            {
              //Check if there is enough residual energy 
              if (bindingEnergy >= ((*fvect)[i])->GetKineticEnergy())
              {
                //Ok, this is a valid secondary: keep it
	        bindingEnergy -= ((*fvect)[i])->GetKineticEnergy();
              }
              else
              {
 	        //Invalid secondary: not enough energy to create it!
 	        //Keep its energy in the local deposit
                delete (*fvect)[i]; 
                (*fvect)[i]=0;
              }
	    } 
          }
        }

        //bindingEnergy has been decreased 
        //by the amount of energy taken away by deexc. products
        if (!statCode)
        {
          fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy);
        }
        else
        {
          fParticleChangeForGamma->SetProposedKineticEnergy(k);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy);
        }

        // TEST //////////////////////////
        // if (secondaryKinetic<0) abort();
        // if (scatteredEnergy<0) abort();
        // if (k-scatteredEnergy-secondaryKinetic-deexSecEnergy<0) abort();
        // if (k-scatteredEnergy<0) abort();     
        /////////////////////////////////

        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
                                                               ionizationShell,
                                                               theIncomingTrack);
    // SI - not useful since low energy of model is 0 eV
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
                                                                          G4double k,
                                                                          G4int shell)
{
    //-- Fast sampling method -----
    G4double proposed_energy;
    G4double random1;
    G4double value_sampling;
    G4double max1;

    do
    {
        proposed_energy = ProposedSampledEnergy(particleDefinition, k, shell); // Proposed energy by inverse function sampling

        max1=0.;

        for(G4double en=0.; en<20.; en+=1.) if(RejectionFunction(particleDefinition, k, en, shell) > max1)
            max1=RejectionFunction(particleDefinition, k, en, shell);

        random1 = G4UniformRand()*max1;

        value_sampling = RejectionFunction(particleDefinition, k, proposed_energy, shell);

    } while(random1 > value_sampling);

    return(proposed_energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::RejectionFunction(G4ParticleDefinition* particleDefinition, 
                                                             G4double k,
                                                             G4double proposed_ws,
                                                             G4int ionizationLevelIndex)
{
    const G4int j=ionizationLevelIndex;
    G4double Bj_energy, alphaConst;
    G4double Ry = 13.6*eV;
    const G4double Gj[5] = {0.99, 1.11, 1.11, 0.52, 1.};

    // const G4double Bj[5] = {12.61*eV, 14.73*eV, 18.55*eV, 32.20*eV, 539.7*eV}; //Ding Paper

    // Following values provided by M. Dingfelder (priv. comm)
    const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

    if (j == 4)
    {
        alphaConst = 0.66;
        //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
        Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
        //---
    }
    else
    {
        alphaConst = 0.64;
        Bj_energy = Bj[ionizationLevelIndex];
    }

    G4double energyTransfer = proposed_ws + Bj_energy;
    proposed_ws/=Bj_energy;

    G4double tau = 0.;
    G4double A_ion = 0.;
    tau = (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
    A_ion = particleDefinition->GetAtomicMass();

    G4double v2;
    G4double beta2;

    if((tau/MeV)<5.447761194e-2)
    {
        v2 = tau / Bj_energy;
        beta2 = 2.*tau / electron_mass_c2;
    }
    // Relativistic
    else
    {
        v2 = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
        beta2 =1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion));
    }

    G4double v = std::sqrt(v2);
    G4double wc = 4.*v2 - 2.*v - (Ry/(4.*Bj_energy));
    G4double rejection_term = 1.+G4Exp(alphaConst*(proposed_ws - wc) / v);
    rejection_term = (1./rejection_term)*CorrectionFactor(particleDefinition,k,ionizationLevelIndex) * Gj[j];
    //* (S/Bj_energy) ; Not needed anymore

    G4bool isHelium = false;

    if (    particleDefinition == protonDef
            || particleDefinition == hydrogenDef
            )
    {
        return(rejection_term);
    }

    else if(particleDefinition->GetAtomicMass() > 4) // anything above Helium
    {
        G4double Z = particleDefinition->GetAtomicNumber();

        G4double x = 100.*std::sqrt(beta2)/gpow->powA(Z, 2./3.);
        G4double Zeffion = Z*(1.-G4Exp(-1.316*x+0.112*x*x-0.0650*x*x*x));
        rejection_term*=Zeffion*Zeffion;
    }

    else if (particleDefinition == alphaPlusPlusDef )
    {
        isHelium = true;
        slaterEffectiveCharge[0]=0.;
        slaterEffectiveCharge[1]=0.;
        slaterEffectiveCharge[2]=0.;
        sCoefficient[0]=0.;
        sCoefficient[1]=0.;
        sCoefficient[2]=0.;
    }

    else if (particleDefinition == alphaPlusDef )
    {
        isHelium = true;
        slaterEffectiveCharge[0]=2.0;
        // The following values are provided by M. Dingfelder (priv. comm)
        slaterEffectiveCharge[1]=2.0;
        slaterEffectiveCharge[2]=2.0;
        //
        sCoefficient[0]=0.7;
        sCoefficient[1]=0.15;
        sCoefficient[2]=0.15;
    }

    else if (particleDefinition == heliumDef )
    {
        isHelium = true;
        slaterEffectiveCharge[0]=1.7;
        slaterEffectiveCharge[1]=1.15;
        slaterEffectiveCharge[2]=1.15;
        sCoefficient[0]=0.5;
        sCoefficient[1]=0.25;
        sCoefficient[2]=0.25;
    }

    if (isHelium)
    {

        G4double zEff = particleDefinition->GetPDGCharge() / eplus + particleDefinition->GetLeptonNumber();

        zEff -= ( sCoefficient[0] * S_1s(k, energyTransfer, slaterEffectiveCharge[0], 1.) +
                  sCoefficient[1] * S_2s(k, energyTransfer, slaterEffectiveCharge[1], 2.) +
                  sCoefficient[2] * S_2p(k, energyTransfer, slaterEffectiveCharge[2], 2.) );

        rejection_term*= zEff * zEff;
    }

    return (rejection_term);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4DNARuddIonisationExtendedModel::ProposedSampledEnergy(G4ParticleDefinition* particle, 
                                                                 G4double k,
                                                                 G4int ionizationLevelIndex)
{

    const G4int j=ionizationLevelIndex;

    G4double A1, B1, C1, D1, E1, A2, B2, C2, D2;
    //G4double alphaConst ;
    G4double Bj_energy;

    // const G4double Bj[5] = {12.61*eV, 14.73*eV, 18.55*eV, 32.20*eV, 539.7*eV}; //Ding Paper
    // Following values provided by M. Dingfelder (priv. comm)

    const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

    if (j == 4)
    {
        //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
        A1 = 1.25;
        B1 = 0.5;
        C1 = 1.00;
        D1 = 1.00;
        E1 = 3.00;
        A2 = 1.10;
        B2 = 1.30;
        C2 = 1.00;
        D2 = 0.00;
        //alphaConst = 0.66;
        //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
        Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
        //---
    }
    else
    {
        //Data For Liquid Water from Dingfelder (Protons in Water)
        A1 = 1.02;
        B1 = 82.0;
        C1 = 0.45;
        D1 = -0.80;
        E1 = 0.38;
        A2 = 1.07;
        //B2 = 14.6; From Ding Paper
        // Value provided by M. Dingfelder (priv. comm)
        B2 = 11.6;
        //
        C2 = 0.60;
        D2 = 0.04;
        //alphaConst = 0.64;

        Bj_energy = Bj[ionizationLevelIndex];
    }

    G4double tau = 0.;
    G4double A_ion = 0.;
    tau = (electron_mass_c2 / particle->GetPDGMass()) * k;

    A_ion = particle->GetAtomicMass();

    G4double v2;
    G4double beta2;
    if((tau/MeV)<5.447761194e-2)
    {
        v2 = tau / Bj_energy;
        beta2 = 2.*tau / electron_mass_c2;
    }
    // Relativistic
    else
    {
        v2 = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
        beta2 =1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion));
    }

    G4double v = std::sqrt(v2);
    //G4double wc = 4.*v2 - 2.*v - (Ry/(4.*Bj_energy));
    G4double L1 = (C1* gpow->powA(v,(D1))) / (1.+ E1*gpow->powA(v, (D1+4.)));
    G4double L2 = C2*gpow->powA(v,(D2));
    G4double H1 = (A1*G4Log(1.+v2)) / (v2+(B1/v2));
    G4double H2 = (A2/v2) + (B2/(v2*v2));
    G4double F1 = L1+H1;
    G4double F2 = (L2*H2)/(L2+H2);

    // ZF. generalized & relativistic version
    G4double maximumEnergy;

    //---- maximum kinetic energy , non relativistic ------
    if( (k/MeV)/(particle->GetPDGMass()/MeV)  <= 0.1 )
    {
        maximumEnergy = 4.* (electron_mass_c2 / particle->GetPDGMass()) * k;
    }
    //---- relativistic -----------------------------------
    else
    {
        G4double gamma = 1./sqrt(1.-beta2);
        maximumEnergy = 2.*electron_mass_c2*(gamma*gamma-1.)/
                (1.+2.*gamma*(electron_mass_c2/particle->GetPDGMass())+pow(electron_mass_c2/particle->GetPDGMass(), 2.) );
    }

    //either it is transfered energy or secondary electron energy ...
    //maximumEnergy-=Bj_energy;

    //-----------------------------------------------------
    G4double wmax = maximumEnergy/Bj_energy;
    G4double c = wmax*(F2*wmax+F1*(2.+wmax))/(2.*(1.+wmax)*(1.+wmax));
    c=1./c; //!!!!!!!!!!! manual calculus leads to  c=1/c
    G4double randVal = G4UniformRand();
    G4double proposed_ws = F1*F1*c*c + 2.*F2*c*randVal - 2.*F1*c*randVal;
    proposed_ws = -F1*c+2.*randVal+std::sqrt(proposed_ws);
    //  proposed_ws = -F1*c+2.*randVal-std::sqrt(proposed_ws);
    proposed_ws/= ( F1*c + F2*c - 2.*randVal );
    proposed_ws*=Bj_energy;

    return(proposed_ws);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::S_1s(G4double t, 
                                                G4double energyTransferred,
                                                G4double slaterEffectiveChg,
                                                G4double shellNumber)
{
    // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
    // Dingfelder, in Chattanooga 2005 proceedings, formula (7)

    G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
    G4double value = 1. - G4Exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::S_2s(G4double t,
                                                G4double energyTransferred,
                                                G4double slaterEffectiveChg,
                                                G4double shellNumber)
{
    // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
    // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

    G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
    G4double value =  1. - G4Exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

    return value;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::S_2p(G4double t, 
                                                G4double energyTransferred,
                                                G4double slaterEffectiveChg,
                                                G4double shellNumber)
{
    // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
    // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

    G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
    G4double value =  1. - G4Exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::R(G4double t,
                                             G4double energyTransferred,
                                             G4double slaterEffectiveChg,
                                             G4double shellNumber)
{
    // tElectron = m_electron / m_alpha * t
    // Dingfelder, in Chattanooga 2005 proceedings, p 4

    G4double tElectron = 0.511/3728. * t;
    // The following values are provided by M. Dingfelder (priv. comm)
    G4double H = 2.*13.60569172 * eV;
    G4double value = std::sqrt ( 2. * tElectron / H ) / ( energyTransferred / H ) *  (slaterEffectiveChg/shellNumber);

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::CorrectionFactor(G4ParticleDefinition* particleDefinition, G4double k, G4int shell) 
{
    // ZF Shortened

    if (particleDefinition == hydrogenDef && shell < 4)
    {
        G4double value = ((G4Log(k/eV)/gpow->logZ(10))-4.2)/0.5;
        // The following values are provided by M. Dingfelder (priv. comm)
        return((0.6/(1+G4Exp(value))) + 0.9);
    }
    else
    {
        return(1.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationExtendedModel::RandomSelect(G4double k)
{   
    G4int level = 0;
 
    // Retrieve data table corresponding to the current particle type

    G4DNACrossSectionDataSet* table = currentTable;

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
    return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::PartialCrossSection(const G4Track& track )
{
    G4double sigma = 0.;

    const G4DynamicParticle* particle = track.GetDynamicParticle();
    G4double k = particle->GetKineticEnergy();

    G4double lowLim = 0;
    G4double highLim = 0;

    const G4String& particleName = particle->GetDefinition()->GetParticleName();

    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(particleName);

    if (pos1 != lowEnergyLimit.end())
    {
        lowLim = pos1->second;
    }

    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);

    if (pos2 != highEnergyLimit.end())
    {
        highLim = pos2->second;
    }

    if (k >= lowLim && k <= highLim)
    {
        std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
        pos = tableData.find(particleName);

        if (pos != tableData.end())
        {
            G4DNACrossSectionDataSet* table = pos->second;
            if (table != 0)
            {
                sigma = table->FindValue(k);
            }
        }
        else
        {
            G4Exception("G4DNARuddIonisationExtendedModel::PartialCrossSection","em0002",
                        FatalException,"Model not applicable to particle type.");
        }
    }

    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::Sum(G4double /* energy */, const G4String& /* particle */)
{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4DNARuddIonisationExtendedModel::GetDNAIonParticleDefinition(const G4ParticleDefinition* particleDefinition)
{
  //for proton, hydrogen, alphas
  if(particleDefinition->GetAtomicMass() <= 4)
  {
    return const_cast<G4ParticleDefinition*>(particleDefinition);
  }
  else{
    auto instance = G4DNAGenericIonsManager::Instance();

    auto PDGEncoding = particleDefinition->GetPDGEncoding();
    if(PDGEncoding == 1000140280){
      return instance->GetIon("silicon");
    }else if(PDGEncoding == 1000260560){
      return instance->GetIon("iron");
    }else if(PDGEncoding == 1000080160){
      return instance->GetIon("oxygen");
    }else if(PDGEncoding == 1000070140){
      return instance->GetIon("nitrogen");
    }else if(PDGEncoding == 1000060120){
      return instance->GetIon("carbon");
    }
    //if there is no DNA particle, get nullptr
    return nullptr;
  }
}
