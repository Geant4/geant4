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
// $Id: G4DNARuddIonisationExtendedModel.cc 104430 2017-05-31 07:43:44Z gcosmo $
// GEANT4 tag $Name:  $
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

//SEB
#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::G4DNARuddIonisationExtendedModel(const G4ParticleDefinition*,
                                                                   const G4String& nam)
    :G4VEmModel(nam),isInitialised(false)
{
    //  nistwater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity = 0;

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
    fAtomDeexcitation = 0;
    fParticleChangeForGamma = 0;

    // Selection of stationary mode

    statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::~G4DNARuddIonisationExtendedModel()
{  
    // Cross section

    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
        G4DNACrossSectionDataSet* table = pos->second;
        delete table;
    }

    // The following removal is forbidden G4VEnergyLossModel takes care of deletion
    // however coverity will signal this as an error
    // if (fAtomDeexcitation) {delete  fAtomDeexcitation;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::Initialise(const G4ParticleDefinition* particle,
                                                  const G4DataVector& /*cuts*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling G4DNARuddIonisationExtendedModel::Initialise()" << G4endl;

    // Energy limits

    G4String fileProton("dna/sigma_ionisation_p_rudd");
    G4String fileHydrogen("dna/sigma_ionisation_h_rudd");
    G4String fileAlphaPlusPlus("dna/sigma_ionisation_alphaplusplus_rudd");
    G4String fileAlphaPlus("dna/sigma_ionisation_alphaplus_rudd");
    G4String fileHelium("dna/sigma_ionisation_he_rudd");
    G4String fileLithium("dna/sigma_ionisation_li_rudd");
    G4String fileBeryllium("dna/sigma_ionisation_be_rudd");
    G4String fileBoron("dna/sigma_ionisation_b_rudd");
    G4String fileCarbon("dna/sigma_ionisation_c_rudd");
    G4String fileNitrogen("dna/sigma_ionisation_n_rudd");
    G4String fileOxygen("dna/sigma_ionisation_o_rudd");
    G4String fileSilicon("dna/sigma_ionisation_si_rudd");
    G4String fileIron("dna/sigma_ionisation_fe_rudd");

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
    G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
    G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
    G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
    G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
    G4ParticleDefinition* heliumDef = instance->GetIon("helium");

    //SEB
    //G4ParticleDefinition* carbonDef = instance->GetIon("carbon");
    //G4ParticleDefinition* nitrogenDef = instance->GetIon("nitrogen");
    //G4ParticleDefinition* oxygenDef = instance->GetIon("oxygen");
    //G4ParticleDefinition* siliconDef = instance->GetIon("silicon");
    //G4ParticleDefinition* ironDef = instance->GetIon("iron");
    G4ParticleDefinition* lithiumDef =  G4IonTable::GetIonTable()->GetIon(3,7);
    G4ParticleDefinition* berylliumDef =  G4IonTable::GetIonTable()->GetIon(4,9);
    G4ParticleDefinition* boronDef =  G4IonTable::GetIonTable()->GetIon(5,11);
    G4ParticleDefinition* carbonDef =  G4IonTable::GetIonTable()->GetIon(6,12);
    G4ParticleDefinition* nitrogenDef =  G4IonTable::GetIonTable()->GetIon(7,14);
    G4ParticleDefinition* oxygenDef =  G4IonTable::GetIonTable()->GetIon(8,16);
    G4ParticleDefinition* siliconDef = G4IonTable::GetIonTable()->GetIon(14,28);
    G4ParticleDefinition* ironDef =  G4IonTable::GetIonTable()->GetIon(26,56);
    //

    G4String proton;
    G4String hydrogen;
    G4String alphaPlusPlus;
    G4String alphaPlus;
    G4String helium;
    G4String lithium;
    G4String beryllium;
    G4String boron;
    G4String carbon;
    G4String nitrogen;
    G4String oxygen;
    G4String silicon;
    G4String iron;

    G4double scaleFactor = 1 * m*m;

    // LIMITS AND DATA

    // **********************************************************************************************

    proton = protonDef->GetParticleName();
    tableFile[proton] = fileProton;
    lowEnergyLimit[proton] = lowEnergyLimitForA[1];
    highEnergyLimit[proton] = 500. * keV;

    // Cross section

    G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableProton->LoadData(fileProton);
    tableData[proton] = tableProton;

    // **********************************************************************************************
    
    hydrogen = hydrogenDef->GetParticleName();
    tableFile[hydrogen] = fileHydrogen;

    lowEnergyLimit[hydrogen] = lowEnergyLimitForA[1];
    highEnergyLimit[hydrogen] = 100. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                           eV,
                                                                           scaleFactor );
    tableHydrogen->LoadData(fileHydrogen);

    tableData[hydrogen] = tableHydrogen;

    // **********************************************************************************************
    
    alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
    tableFile[alphaPlusPlus] = fileAlphaPlusPlus;

    lowEnergyLimit[alphaPlusPlus] = lowEnergyLimitForA[4];
    highEnergyLimit[alphaPlusPlus] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                                eV,
                                                                                scaleFactor );
    tableAlphaPlusPlus->LoadData(fileAlphaPlusPlus);

    tableData[alphaPlusPlus] = tableAlphaPlusPlus;

    // **********************************************************************************************
    
    alphaPlus = alphaPlusDef->GetParticleName();
    tableFile[alphaPlus] = fileAlphaPlus;

    lowEnergyLimit[alphaPlus] = lowEnergyLimitForA[4];
    highEnergyLimit[alphaPlus] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                            eV,
                                                                            scaleFactor );
    tableAlphaPlus->LoadData(fileAlphaPlus);
    tableData[alphaPlus] = tableAlphaPlus;

    // **********************************************************************************************
    
    helium = heliumDef->GetParticleName();
    tableFile[helium] = fileHelium;

    lowEnergyLimit[helium] = lowEnergyLimitForA[4];
    highEnergyLimit[helium] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableHelium->LoadData(fileHelium);
    tableData[helium] = tableHelium;

    // **********************************************************************************************
    
    lithium = lithiumDef->GetParticleName();
    tableFile[lithium] = fileLithium;

    //SEB
    //lowEnergyLimit[carbon] = lowEnergyLimitForA[5] * particle->GetAtomicMass();
    //highEnergyLimit[carbon] = 1e6* particle->GetAtomicMass() * MeV;
    lowEnergyLimit[lithium] = 0.5*7*MeV;
    highEnergyLimit[lithium] = 1e6*7*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableLithium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableLithium->LoadData(fileLithium);
    tableData[lithium] = tableLithium;

    // **********************************************************************************************
    
    beryllium = berylliumDef->GetParticleName();
    tableFile[beryllium] = fileBeryllium;

    //SEB
    //lowEnergyLimit[carbon] = lowEnergyLimitForA[5] * particle->GetAtomicMass();
    //highEnergyLimit[carbon] = 1e6* particle->GetAtomicMass() * MeV;
    lowEnergyLimit[beryllium] = 0.5*9*MeV;
    highEnergyLimit[beryllium] = 1e6*9*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableBeryllium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableBeryllium->LoadData(fileBeryllium);
    tableData[beryllium] = tableBeryllium;

    // **********************************************************************************************
    
    boron = boronDef->GetParticleName();
    tableFile[boron] = fileBoron;

    //SEB
    //lowEnergyLimit[carbon] = lowEnergyLimitForA[5] * particle->GetAtomicMass();
    //highEnergyLimit[carbon] = 1e6* particle->GetAtomicMass() * MeV;
    lowEnergyLimit[boron] = 0.5*11*MeV;
    highEnergyLimit[boron] = 1e6*11*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableBoron = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableBoron->LoadData(fileBoron);
    tableData[boron] = tableBoron;

    // **********************************************************************************************
    
    carbon = carbonDef->GetParticleName();
    tableFile[carbon] = fileCarbon;

    //SEB
    //lowEnergyLimit[carbon] = lowEnergyLimitForA[5] * particle->GetAtomicMass();
    //highEnergyLimit[carbon] = 1e6* particle->GetAtomicMass() * MeV;
    lowEnergyLimit[carbon] = 0.5*12*MeV;
    highEnergyLimit[carbon] = 1e6*12*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableCarbon = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableCarbon->LoadData(fileCarbon);
    tableData[carbon] = tableCarbon;

    // **********************************************************************************************
    
    oxygen = oxygenDef->GetParticleName();
    tableFile[oxygen] = fileOxygen;

    //SEB
    //lowEnergyLimit[oxygen] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    //highEnergyLimit[oxygen] = 1e6* particle->GetAtomicMass()* MeV;
    lowEnergyLimit[oxygen] = 0.5*16*MeV;
    highEnergyLimit[oxygen] = 1e6*16*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableOxygen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                         eV,
                                                                         scaleFactor );
    tableOxygen->LoadData(fileOxygen);
    tableData[oxygen] = tableOxygen;

    // **********************************************************************************************
    
    nitrogen = nitrogenDef->GetParticleName();
    tableFile[nitrogen] = fileNitrogen;

    //SEB
    //lowEnergyLimit[nitrogen] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    //highEnergyLimit[nitrogen] = 1e6* particle->GetAtomicMass()* MeV;
    lowEnergyLimit[nitrogen] = 0.5*14*MeV;
    highEnergyLimit[nitrogen] = 1e6*14*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableNitrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                           eV,
                                                                           scaleFactor );
    tableNitrogen->LoadData(fileNitrogen);
    tableData[nitrogen] = tableNitrogen;

    // **********************************************************************************************

    silicon = siliconDef->GetParticleName();
    tableFile[silicon] = fileSilicon;
    
    //lowEnergyLimit[silicon] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    //highEnergyLimit[silicon] = 1e6* particle->GetAtomicMass()* MeV;
    lowEnergyLimit[silicon] = 0.5*28*MeV;
    highEnergyLimit[silicon] = 1e6*28*MeV;
    //
    
    // Cross section
    
    G4DNACrossSectionDataSet* tableSilicon = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                          eV,
                                                                          scaleFactor );
    tableSilicon->LoadData(fileSilicon);
    tableData[silicon] = tableSilicon;
     
    // **********************************************************************************************
    
    iron = ironDef->GetParticleName();
    tableFile[iron] = fileIron;

    //SEB
    //lowEnergyLimit[iron] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    //highEnergyLimit[iron] = 1e6* particle->GetAtomicMass()* MeV;
    lowEnergyLimit[iron] = 0.5*56*MeV;
    highEnergyLimit[iron] = 1e6*56*MeV;
    //

    // Cross section

    G4DNACrossSectionDataSet* tableIron = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
                                                                       eV,
                                                                       scaleFactor );
    tableIron->LoadData(fileIron);
    tableData[iron] = tableIron;

    // **********************************************************************************************

    //SEB: not anymore
    // ZF Following lines can be replaced by:
    //SetLowEnergyLimit(lowEnergyLimit[particle->GetParticleName()]);
    //SetHighEnergyLimit(highEnergyLimit[particle->GetParticleName()]);
    // at least for HZE
    
  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if (particle==hydrogenDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[hydrogen]);
    SetHighEnergyLimit(highEnergyLimit[hydrogen]);
  }

  if (particle==heliumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[helium]);
    SetHighEnergyLimit(highEnergyLimit[helium]);
  }

  if (particle==alphaPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  if (particle==alphaPlusPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlusPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlusPlus]);
  }

  if (particle==lithiumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[lithium]);
    SetHighEnergyLimit(highEnergyLimit[lithium]);
  }

  if (particle==berylliumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[beryllium]);
    SetHighEnergyLimit(highEnergyLimit[beryllium]);
  }

  if (particle==boronDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[boron]);
    SetHighEnergyLimit(highEnergyLimit[boron]);
  }

  if (particle==carbonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[carbon]);
    SetHighEnergyLimit(highEnergyLimit[carbon]);
  }

  if (particle==nitrogenDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[nitrogen]);
    SetHighEnergyLimit(highEnergyLimit[nitrogen]);
  }

  if (particle==oxygenDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[oxygen]);
    SetHighEnergyLimit(highEnergyLimit[oxygen]);
  }

  if (particle==siliconDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[silicon]);
    SetHighEnergyLimit(highEnergyLimit[silicon]);
  }

  if (particle==ironDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[iron]);
    SetHighEnergyLimit(highEnergyLimit[iron]);
  }

    //----------------------------------------------------------------------

    if( verboseLevel>0 )
    {
        G4cout << "Rudd ionisation model is initialized " << G4endl
               << "Energy range: "
               << LowEnergyLimit() / eV << " eV - "
               << HighEnergyLimit() / keV << " keV for "
               << particle->GetParticleName()
               << G4endl;
    }

    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

    //

    fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

    if (isInitialised) { return; }
    fParticleChangeForGamma = GetParticleChangeForGamma();
    isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARuddIonisationExtendedModel::CrossSectionPerVolume(const G4Material* material,
                                                                 const G4ParticleDefinition* particleDefinition,
                                                                 G4double k,
                                                                 G4double,
                                                                 G4double)
{
    //SEB: particleDefinition->GetParticleName() is for eg. Fe56
    //     particleDefinition->GetPDGMass() is correct
    //     particleDefinition->GetAtomicNumber() is correct

    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNARuddIonisationExtendedModel" << G4endl;

    // Calculate total cross section for model

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();

    if (
            particleDefinition != G4Proton::ProtonDefinition()
            &&
            particleDefinition != instance->GetIon("hydrogen")
            &&
            particleDefinition != instance->GetIon("alpha++")
            &&
            particleDefinition != instance->GetIon("alpha+")
            &&
            particleDefinition != instance->GetIon("helium")
            &&
            //SEB
            //particleDefinition != instance->GetIon("carbon")
            //&&
            //particleDefinition != instance->GetIon("nitrogen")
            //&&
            //particleDefinition != instance->GetIon("oxygen")
            //&&
            //particleDefinition != instance->GetIon("iron")
            particleDefinition != G4IonTable::GetIonTable()->GetIon(3,7)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(4,9)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(5,11)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(6,12)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(7,14)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(8,16)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(14,28)
            &&
            particleDefinition != G4IonTable::GetIonTable()->GetIon(26,56)
            //
            )

        return 0;

    G4double lowLim = 0;

    if (     particleDefinition == G4Proton::ProtonDefinition()
             ||  particleDefinition == instance->GetIon("hydrogen")
             )

        lowLim = lowEnergyLimitOfModelForA[1];

    else if (     particleDefinition == instance->GetIon("alpha++")
                  ||       particleDefinition == instance->GetIon("alpha+")
                  ||       particleDefinition == instance->GetIon("helium")
                  )

        lowLim = lowEnergyLimitOfModelForA[4];

    else lowLim = lowEnergyLimitOfModelForA[5];

    G4double highLim = 0;
    G4double sigma=0;


    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0)
//    if (material == nistwater || material->GetBaseMaterial() == nistwater)
    {
        const G4String& particleName = particleDefinition->GetParticleName();

        std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
        pos2 = highEnergyLimit.find(particleName);

        if (pos2 != highEnergyLimit.end())
        {
            highLim = pos2->second;
        }

        if (k <= highLim)
        {

            //SI : XS must not be zero otherwise sampling of secondaries method ignored

            if (k < lowLim) k = lowLim;

            //

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
                G4Exception("G4DNARuddIonisationExtendedModel::CrossSectionPerVolume","em0002",
                            FatalException,"Model not applicable to particle type.");
            }

        } // if (k >= lowLim && k < highLim)

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "G4DNARuddIonisationExtendedModel - XS INFO START" << G4endl;
            G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
            G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
            //G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
            G4cout << "G4DNARuddIonisationExtendedModel - XS INFO END" << G4endl;

        }

    } // if (waterMaterial)

    return sigma*waterDensity;
//    return sigma*material->GetAtomicNumDensityVector()[1];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                         const G4MaterialCutsCouple* couple,
                                                         const G4DynamicParticle* particle,
                                                         G4double,
                                                         G4double)
{
    //SEB: particle->GetDefinition()->GetParticleName() is for eg. Fe56
    //     particle->GetDefinition()->GetPDGMass() is correct
    //     particle->GetDefinition()->GetAtomicNumber() is correct
    //     particle->GetDefinition()->GetAtomicMass() is correct

    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNARuddIonisationExtendedModel" << G4endl;

    G4double lowLim = 0;
    G4double highLim = 0;

    // ZF: the following line summarizes the commented part

    if(particle->GetDefinition()->GetAtomicMass() <= 4) lowLim = killBelowEnergyForA[particle->GetDefinition()->GetAtomicMass()];
    
    else lowLim = killBelowEnergyForA[5]*particle->GetDefinition()->GetAtomicMass();

    /* 

    if(particle->GetDefinition()->GetAtomicMass() >= 5) lowLim = killBelowEnergyForA[5]*particle->GetDefinition()->GetAtomicMass();
    
    if (     particle->GetDefinition() == G4Proton::ProtonDefinition()
       ||  particle->GetDefinition() == instance->GetIon("hydrogen")
       )

       lowLim = killBelowEnergyForA[1];
       
    if (     particle->GetDefinition() == instance->GetIon("alpha++")
       ||  particle->GetDefinition() == instance->GetIon("alpha+")
       ||  particle->GetDefinition() == instance->GetIon("helium")
       )

       lowLim = killBelowEnergyForA[4];
    
    */
    //

    G4double k = particle->GetKineticEnergy();

    const G4String& particleName = particle->GetDefinition()->GetParticleName();

    // SI - the following is useless since lowLim is already defined
    /*
    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(particleName);

    if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }
    */

    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);

    if (pos2 != highEnergyLimit.end())highLim = pos2->second;

    if (k >= lowLim && k <= highLim) 
    // SI: no strict limits, like in the non extended version of the model
    {
        G4ParticleDefinition* definition = particle->GetDefinition();
        G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
        /*
        G4double particleMass = definition->GetPDGMass();
        G4double totalEnergy = k + particleMass;
        G4double pSquare = k*(totalEnergy+particleMass);
        G4double totalMomentum = std::sqrt(pSquare);
        */

        G4int ionizationShell = RandomSelect(k,particleName);

        // sample deexcitation
        // here we assume that H_{2}O electronic levels are the same as Oxygen.
        // this can be considered true with a rough 10% error in energy on K-shell,

        G4int secNumberInit = 0;  // need to know at a certain point the energy of secondaries
        G4int secNumberFinal = 0; // So I'll make the diference and then sum the energies
        G4double bindingEnergy = 0;
        bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

        //SI: additional protection if tcs interpolation method is modified
        if (k<bindingEnergy) return;
        //

        G4int Z = 8;

        if(fAtomDeexcitation) {
            G4AtomicShellEnumerator as = fKShell;

            if (ionizationShell <5 && ionizationShell >1)
            {
                as = G4AtomicShellEnumerator(4-ionizationShell);
            }
            else if (ionizationShell <2)
            {
                as = G4AtomicShellEnumerator(3);
            }

            //	DEBUG
            //	if (ionizationShell == 4) {
            //
            //	  G4cout << "Z: " << Z << " as: " << as
            //               << " ionizationShell: " << ionizationShell << " bindingEnergy: "<< bindingEnergy/eV << G4endl;
            //	  G4cout << "Press <Enter> key to continue..." << G4endl;
            //	  G4cin.ignore();
            //	}


            const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
            secNumberInit = fvect->size();
            fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
            secNumberFinal = fvect->size();
        }

        G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);

	G4ThreeVector deltaDirection = 
	  GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic, 
							    Z, ionizationShell,
							    couple->GetMaterial());

        // SI: the following lines are not needed anymore
        /*
        G4double cosTheta = 0.;
        G4double phi = 0.;
        RandomizeEjectedElectronDirection(definition, k,secondaryKinetic, cosTheta, phi, ionizationShell);

        G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
        G4double dirX = sinTheta*std::cos(phi);
        G4double dirY = sinTheta*std::sin(phi);
        G4double dirZ = cosTheta;
        G4ThreeVector deltaDirection(dirX,dirY,dirZ);
        deltaDirection.rotateUz(primaryDirection);
        */
  
        // Ignored for ions on electrons
        /*
        G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));

        G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
        G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
        G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
        G4double finalMomentum = std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
        finalPx /= finalMomentum;
        finalPy /= finalMomentum;
        finalPz /= finalMomentum;

        G4ThreeVector direction;
        direction.set(finalPx,finalPy,finalPz);

        fParticleChangeForGamma->ProposeMomentumDirection(direction.unit()) ;
        */

        fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
        G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;
        G4double deexSecEnergy = 0;
        for (G4int j=secNumberInit; j < secNumberFinal; j++) {

            deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();

        }

        if (!statCode)
        {
          fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy-secondaryKinetic-deexSecEnergy);
        }
        else
        {
          fParticleChangeForGamma->SetProposedKineticEnergy(k);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy);
        }

        G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
        fvect->push_back(dp);

        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
                                                               ionizationShell,
                                                               theIncomingTrack);
    }

    // SI - not useful since low energy of model is 0 eV

    if (k < lowLim)
    {
        fParticleChangeForGamma->SetProposedKineticEnergy(0.);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }

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

// The following section is not used anymore but is kept for memory
// GetAngularDistribution()->SampleDirectionForShell is used instead

/*
void G4DNARuddIonisationExtendedModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
                                                                         G4double k,
                                                                         G4double secKinetic,
                                                                         G4double & cosTheta,
                                                                         G4double & phi,
                                                                         G4int shell )
{
    G4double maxSecKinetic = 0.;
    G4double maximumEnergyTransfer = 0.;

    // ZF. generalized & relativistic version

    if( (k/MeV)/(particleDefinition->GetPDGMass()/MeV)  <= 0.1 )
    {
        maximumEnergyTransfer= 4.* (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
        maximumEnergyTransfer+=waterStructure.IonisationEnergy(shell);
    }
    else
    {
        G4double approx_nuc_number = particleDefinition->GetPDGMass() / proton_mass_c2;
        G4double en_per_nucleon = k/approx_nuc_number;
        G4double beta2 = 1. - 1./pow( (1.+(en_per_nucleon/electron_mass_c2)*(electron_mass_c2/proton_mass_c2)), 2.);
        G4double gamma = 1./sqrt(1.-beta2);
        maximumEnergyTransfer = 2.*electron_mass_c2*(gamma*gamma-1.)/(1.+2.*gamma*(electron_mass_c2/particleDefinition->GetPDGMass())+pow(electron_mass_c2/particleDefinition->GetPDGMass(), 2.) );
        maximumEnergyTransfer+=waterStructure.IonisationEnergy(shell);
    }

    maxSecKinetic = maximumEnergyTransfer-waterStructure.IonisationEnergy(shell);

    phi = twopi * G4UniformRand();

    if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
    else cosTheta = (2.*G4UniformRand())-1.;

}
*/

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
    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
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

    if (    particleDefinition == G4Proton::ProtonDefinition()
            || particleDefinition == instance->GetIon("hydrogen")
            )
    {
        return(rejection_term);
    }

    else if(particleDefinition->GetAtomicMass() > 4) // anything above Helium
    {
        G4double Z = particleDefinition->GetAtomicNumber();

        G4double x = 100.*std::sqrt(beta2)/std::pow(Z,(2./3.));
        G4double Zeffion = Z*(1.-G4Exp(-1.316*x+0.112*x*x-0.0650*x*x*x));
        rejection_term*=Zeffion*Zeffion;
    }

    else if (particleDefinition == instance->GetIon("alpha++") )
    {
        isHelium = true;
        slaterEffectiveCharge[0]=0.;
        slaterEffectiveCharge[1]=0.;
        slaterEffectiveCharge[2]=0.;
        sCoefficient[0]=0.;
        sCoefficient[1]=0.;
        sCoefficient[2]=0.;
    }

    else if (particleDefinition == instance->GetIon("alpha+") )
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

    else if (particleDefinition == instance->GetIon("helium") )
    {
        isHelium = true;
        slaterEffectiveCharge[0]=1.7;
        slaterEffectiveCharge[1]=1.15;
        slaterEffectiveCharge[2]=1.15;
        sCoefficient[0]=0.5;
        sCoefficient[1]=0.25;
        sCoefficient[2]=0.25;
    }

   //    if (    particleDefinition == instance->GetIon("helium")
   //            || particleDefinition == instance->GetIon("alpha+")
   //            || particleDefinition == instance->GetIon("alpha++")
   //            )

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
    G4double L1 = (C1* std::pow(v,(D1))) / (1.+ E1*std::pow(v, (D1+4.)));
    G4double L2 = C2*std::pow(v,(D2));
    G4double H1 = (A1*std::log(1.+v2)) / (v2+(B1/v2));
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
    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();

    if (particleDefinition == instance->GetIon("hydrogen") && shell < 4)
    {
        G4double value = (std::log10(k/eV)-4.2)/0.5;
        // The following values are provided by M. Dingfelder (priv. comm)
        return((0.6/(1+G4Exp(value))) + 0.9);
    }
    else
    {
        return(1.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationExtendedModel::RandomSelect(G4double k, const G4String& particle )
{   

    G4int level = 0;

    // Retrieve data table corresponding to the current particle type

    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    pos = tableData.find(particle);

    if (pos != tableData.end())
    {
        G4DNACrossSectionDataSet* table = pos->second;

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
        G4Exception("G4DNARuddIonisationExtendedModel::RandomSelect","em0002",
                    FatalException,"Model not applicable to particle type.");
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

