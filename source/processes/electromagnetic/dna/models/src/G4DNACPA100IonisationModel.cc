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
// CPA100 ionisation model class for electrons
//
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#include "G4DNACPA100IonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100IonisationModel::G4DNACPA100IonisationModel(const G4ParticleDefinition*,
                                                       const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 )
    {
        G4cout << "CPA100 ionisation model is constructed " << G4endl;
    }

    // Mark this model as "applicable" for atomic deexcitation
    SetDeexcitationFlag(true);
    fAtomDeexcitation = 0;
    fParticleChangeForGamma = 0;
    fpMolWaterDensity = 0;
    
    // Selection of computation method    
    
    // useDcs = true if usage of dcs for sampling of secondaries
    // useDcs = false if usage of composition sampling (DEFAULT)
    
    useDcs = true;
    
    // if useDcs is true, one has the following choice
    // fasterCode = true for usage of cumulated dcs (DEFAULT)
    // fasterCode = false for usage of non-cumulated dcs
    
    fasterCode = true;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100IonisationModel::~G4DNACPA100IonisationModel()
{  
    // Cross section

    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
        G4DNACrossSectionDataSet* table = pos->second;
        delete table;
    }

    // Final state

    eVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100IonisationModel::Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& /*cuts*/)
{

    if (verboseLevel > 3)
        G4cout << "Calling G4DNACPA100IonisationModel::Initialise()" << G4endl;

    // Energy limits

    // The following file is proved by M. Terrissol et al. (sigion3)

    G4String fileElectron("dna/sigma_ionisation_e_cpa100_form_rel");

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();

    G4String electron;

    G4double scaleFactor = 1.e-20 * m*m;

    char *path = getenv("G4LEDATA");

    // *** ELECTRON

    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    lowEnergyLimit[electron] = 11. * eV; // Default - for composition sampling only
    //lowEnergyLimit[electron] = 10.985 * eV; // For composition sampling only
    if (useDcs) lowEnergyLimit[electron] = 11 * eV; // For dcs usage, they start at 11 eV
    
    highEnergyLimit[electron] = 255955. * eV;
    //highEnergyLimit[electron] = 1. * MeV; // Ionisation model goes up to 1 MeV but not other CPA100 models

    // Cross section
    
    G4DNACrossSectionDataSet* tableE = 
     new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    
    //G4DNACrossSectionDataSet* tableE =
    // new G4DNACrossSectionDataSet(new G4DNACPA100LogLogInterpolation, eV,scaleFactor );
    
    tableE->LoadData(fileElectron);

    tableData[electron] = tableE;
    
    // Final state
    
    // ******************************
    
    if (useDcs)
    {
    
    std::ostringstream eFullFileName;
   
    if (fasterCode)  eFullFileName << path << "/dna/sigmadiff_cumulated_ionisation_e_cpa100_rel.dat"; 

    if (!fasterCode) eFullFileName << path << "/dna/sigmadiff_ionisation_e_cpa100_rel.dat";
    
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

    if (!eDiffCrossSection)
    {
        if (fasterCode)  G4Exception("G4DNACPA100IonisationModel::Initialise","em0003",
                         FatalException,"Missing data file:/dna/sigmadiff_cumulated_ionisation_e_cpa100_rel.dat");
    
        if (!fasterCode) G4Exception("G4DNACPA100IonisationModel::Initialise","em0003",
                         FatalException,"Missing data file:/dna/sigmadiff_ionisation_e_cpa100_rel.dat");
    }

    // Clear the arrays for re-initialization case (MT mode)
    // March 25th, 2014 - Vaclav Stepan, Sebastien Incerti

    eTdummyVec.clear();
    eVecm.clear();
    eProbaShellMap->clear();
    eDiffCrossSectionData->clear();
    eNrjTransfData->clear();

    //

    eTdummyVec.push_back(0.);
    while(!eDiffCrossSection.eof())
    {
        double tDummy;
        double eDummy;
        eDiffCrossSection>>tDummy>>eDummy;
        if (tDummy != eTdummyVec.back()) eTdummyVec.push_back(tDummy);
        for (int j=0; j<5; j++)
        {
            eDiffCrossSection>>eDiffCrossSectionData[j][tDummy][eDummy];

            if (fasterCode) 
            {
              eNrjTransfData[j][tDummy][eDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
              eProbaShellMap[j][tDummy].push_back(eDiffCrossSectionData[j][tDummy][eDummy]);
            }

     // SI - only if eof is not reached
     if (!eDiffCrossSection.eof() && !fasterCode) eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;
            
     if (!fasterCode) eVecm[tDummy].push_back(eDummy);

       }
    }

    //

    } // end of if (useDcs)
   
    // ******************************
 
    //
    
    if (particle==electronDef)
    {
        SetLowEnergyLimit(lowEnergyLimit[electron]);
        SetHighEnergyLimit(highEnergyLimit[electron]);
    }

    if( verboseLevel>0 )
    {
        G4cout << "CPA100 ionisation model is initialized " << G4endl
               << "Energy range: "
               << LowEnergyLimit() / eV << " eV - "
               << HighEnergyLimit() / keV << " keV for "
               << particle->GetParticleName()
               << G4endl;
    }

    // Initialize water density pointer
    fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

    //
    fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

    if (isInitialised) { return; }
    fParticleChangeForGamma = GetParticleChangeForGamma();
    isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100IonisationModel::CrossSectionPerVolume( const G4Material* material,
                                                            const G4ParticleDefinition* particleDefinition,
                                                            G4double ekin,
                                                            G4double,
                                                            G4double)
{

    if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNACPA100IonisationModel" << G4endl;

    if (particleDefinition != G4Electron::ElectronDefinition()) return 0;

    // Calculate total cross section for model

    G4double lowLim = 0;
    G4double highLim = 0;
    G4double sigma=0;

    G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0)
        
    {
        const G4String& particleName = particleDefinition->GetParticleName();

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

        if (ekin > lowLim && ekin < highLim)
        {
            std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
            pos = tableData.find(particleName);

            if (pos != tableData.end())
            {
                G4DNACrossSectionDataSet* table = pos->second;
                if (table != 0)
                {
                    sigma = table->FindValue(ekin);
                }
            }
            else
            {
                G4Exception("G4DNACPA100IonisationModel::CrossSectionPerVolume","em0002",
                            FatalException,"Model not applicable to particle type.");
            }
        }

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "G4DNACPA100IonisationModel - XS INFO START" << G4endl;
            G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
            G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
            G4cout << "G4DNACPA100IonisationModel - XS INFO END" << G4endl;
        }

    } // if (waterMaterial)

    return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100IonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                 const G4MaterialCutsCouple* ,//must be set!
                                                 const G4DynamicParticle* particle,
                                                 G4double,
                                                 G4double)
{
    if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNACPA100IonisationModel" << G4endl;

    G4double lowLim = 0;
    G4double highLim = 0;

    G4double k = particle->GetKineticEnergy();

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

    if (k >= lowLim && k < highLim)
    {
        G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
        G4double particleMass = particle->GetDefinition()->GetPDGMass();
        G4double totalEnergy = k + particleMass;
        G4double pSquare = k * (totalEnergy + particleMass);
        G4double totalMomentum = std::sqrt(pSquare);

        G4int ionizationShell = -1;
 
        ionizationShell = RandomSelect(k,particleName); 

        //SI: PROTECTION FOR G4LOGLOGINTERPOLATION ON UPPER VALUE 
        if (k<waterStructure.IonisationEnergy(ionizationShell)) { return; } 
      
        // AM: sample deexcitation
        // here we assume that H_{2}O electronic levels are the same of Oxigen.
        // this can be considered true with a rough 10% error in energy on K-shell,

        G4int secNumberInit = 0;  // need to know at a certain point the enrgy of secondaries
        G4int secNumberFinal = 0; // So I'll make the diference and then sum the energies

        G4double bindingEnergy = 0;
        bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

        if(fAtomDeexcitation) {

            G4int Z = 8;
            G4AtomicShellEnumerator as = fKShell;

            if (ionizationShell <5 && ionizationShell >1)
            {
                as = G4AtomicShellEnumerator(4-ionizationShell);
            }
            else if (ionizationShell <2)
            {
                as = G4AtomicShellEnumerator(3);
            }

            // FOR DEBUG ONLY
            // if (ionizationShell == 4) {
            //
            //   G4cout << "Z: " << Z << " as: " << as
            //               << " ionizationShell: " << ionizationShell << " bindingEnergy: "<< bindingEnergy/eV << G4endl;
            //        G4cout << "Press <Enter> key to continue..." << G4endl;
            //   G4cin.ignore();
            // }

            const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
            secNumberInit = fvect->size();
            fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
            secNumberFinal = fvect->size();
        }

        G4double secondaryKinetic=-1000*eV;

        if (useDcs && !fasterCode)
          secondaryKinetic = RandomizeEjectedElectronEnergy(particle->GetDefinition(),k,ionizationShell);

        if (useDcs && fasterCode) 
          secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulatedDcs(particle->GetDefinition(),k,ionizationShell);

        if (!useDcs)
          secondaryKinetic = RandomizeEjectedElectronEnergyFromCompositionSampling(particle->GetDefinition(),k,ionizationShell);

        // Quick test
        /*
        FILE* myFile;
        myFile=fopen("nrj.txt","a");
        fprintf(myFile,"%e\n", secondaryKinetic/eV );
        fclose(myFile);
        */

        G4double cosTheta = 0.;
        G4double phi = 0.;
        RandomizeEjectedElectronDirection(particle->GetDefinition(), k,secondaryKinetic, cosTheta, phi);

        G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
        G4double dirX = sinTheta*std::cos(phi);
        G4double dirY = sinTheta*std::sin(phi);
        G4double dirZ = cosTheta;
        G4ThreeVector deltaDirection(dirX,dirY,dirZ);
        deltaDirection.rotateUz(primaryDirection);

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

        else fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection) ;

        // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
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
        
        // SI - 29/03/2014
        if (secondaryKinetic>0) 
        {
          G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
          fvect->push_back(dp);
        } 
        //

        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
                                                               ionizationShell,
                                                               theIncomingTrack);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
                                                                  G4double k, G4int shell)
{
    // G4cout << "*** SLOW computation for " << " " << particleDefinition->GetParticleName() << G4endl;

    if (particleDefinition == G4Electron::ElectronDefinition())
    {
        G4double maximumEnergyTransfer=0.;
        if ((k+waterStructure.IonisationEnergy(shell))/2. > k) maximumEnergyTransfer=k;
        else maximumEnergyTransfer = (k+waterStructure.IonisationEnergy(shell))/2.;

        // SI : original method
        /*
         G4double crossSectionMaximum = 0.;
         for(G4double value=waterStructure.IonisationEnergy(shell); value<=maximumEnergyTransfer; value+=0.1*eV)
         {
           G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
           if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
         }
        */

        // SI : alternative method

        G4double crossSectionMaximum = 0.;

        G4double minEnergy = waterStructure.IonisationEnergy(shell);
        G4double maxEnergy = maximumEnergyTransfer;
                
        // nEnergySteps can be optimized - 100 by default
        G4int nEnergySteps = 50;

        // *** METHOD 1 
        // FOR SLOW COMPUTATION ONLY
        /*   
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
        */

        // *** METHOD 2 : Faster method for CPA100 only since DCS is monotonously decreasing
        // FOR SLOW COMPUTATION ONLY
        
        G4double value(minEnergy);
        G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
        G4int step(nEnergySteps);
        G4double differentialCrossSection = 0.;
        while (step>0)
        {
            step--;
            differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
            if(differentialCrossSection >0) 
            {
              crossSectionMaximum=differentialCrossSection;
              break;
            }
            value*=stpEnergy;
        }
        
        //

        G4double secondaryElectronKineticEnergy=0.;
        do
        {
            secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer-waterStructure.IonisationEnergy(shell));
        } while(G4UniformRand()*crossSectionMaximum >
                DifferentialCrossSection(particleDefinition, k/eV,
                 (secondaryElectronKineticEnergy+waterStructure.IonisationEnergy(shell))/eV,shell));

        return secondaryElectronKineticEnergy;

    }

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNACPA100IonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition*, 
                                                                 G4double k,
                                                                 G4double secKinetic,
                                                                 G4double & cosTheta,
                                                                 G4double & phi )
{

    phi = twopi * G4UniformRand();
    G4double sin2O = (1.-secKinetic/k) / (1.+secKinetic/(2.*electron_mass_c2));
    cosTheta = std::sqrt(1.-sin2O);
        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition, 
                                                          G4double k,
                                                          G4double energyTransfer,
                                                          G4int ionizationLevelIndex)
{
    G4double sigma = 0.;

    if (energyTransfer >= waterStructure.IonisationEnergy(ionizationLevelIndex))
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

                xs11 = eDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE11];
                xs12 = eDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE12];
                xs21 = eDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE21];
                xs22 = eDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE22];
  
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

G4double G4DNACPA100IonisationModel::Interpolate(      G4double e1, 
                                                     G4double e2,
                                                     G4double e,
                                                     G4double xs1,
                                                     G4double xs2)
{

    G4double value = 0.;

    // Log-log interpolation by default
    
    if (e1!=0 && e2!=0  && (std::log10(e2)-std::log10(e1)) !=0 && !fasterCode && useDcs)
    {  
      G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
      G4double b = std::log10(xs2) - a*std::log10(e2);
      G4double sigma = a*std::log10(e) + b;
      value = (std::pow(10.,sigma));
    }
        
    // Switch to lin-lin interpolation
    /*  
    if ((e2-e1)!=0)
    {
      G4double d1 = xs1;
      G4double d2 = xs2;
      value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    }
    */
    
    // Switch to log-lin interpolation for faster code
    
    if ((e2-e1)!=0 && xs1 !=0 && xs2 !=0 && fasterCode && useDcs )
    {
      G4double d1 = std::log10(xs1);
      G4double d2 = std::log10(xs2);
      value = std::pow(10.,(d1 + (d2 - d1)*(e - e1)/ (e2 - e1)) );
    }

    // Switch to lin-lin interpolation for faster code
    // in case one of xs1 or xs2 (=cum proba) value is zero

    if ((e2-e1)!=0 && (xs1 ==0 || xs2 ==0) && fasterCode && useDcs )
    {
      G4double d1 = xs1;
      G4double d2 = xs2;
      value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    }

    /*
    G4cout 
    << e1 << " "
    << e2 << " "
    << e  << " "
    << xs1 << " "
    << xs2 << " "
    << value
    << G4endl;
    */
    
    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::QuadInterpolator(G4double e11, G4double e12, 
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

G4int G4DNACPA100IonisationModel::RandomSelect(G4double k, const G4String& particle )
{   
    G4int level = 0;

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

 //Verification
 /*
 G4double tmp=200*keV;
 G4cout <<  table->GetComponent(0)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(1)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(2)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(3)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(4)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  
 table->GetComponent(0)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(1)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(2)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(3)->FindValue(tmp)/(1e-20*m*m) 
 << G4endl;
 abort();
 */
 //
 //Dump
 //
 /*
 G4double minEnergy = 10.985  * eV;
 G4double maxEnergy = 255955. * eV;
 G4int nEnergySteps = 1000;
 G4double energy(minEnergy);
 G4double stpEnergy(std::pow(maxEnergy/energy, 1./static_cast<G4double>(nEnergySteps-1)));
 G4int step(nEnergySteps);
 system ("rm -rf ionisation-cpa100.out");
 FILE* myFile=fopen("ionisation-cpa100.out","a");
 while (step>0)
 {
   step--;
   fprintf (myFile,"%16.9le %16.9le %16.9le %16.9le %16.9le %16.9le %16.9le \n",
   energy/eV,
   table->GetComponent(0)->FindValue(energy)/(1e-20*m*m), 
   table->GetComponent(1)->FindValue(energy)/(1e-20*m*m), 
   table->GetComponent(2)->FindValue(energy)/(1e-20*m*m), 
   table->GetComponent(3)->FindValue(energy)/(1e-20*m*m), 
   table->GetComponent(4)->FindValue(energy)/(1e-20*m*m), 
   table->GetComponent(0)->FindValue(energy)/(1e-20*m*m)+ 
   table->GetComponent(1)->FindValue(energy)/(1e-20*m*m)+ 
   table->GetComponent(2)->FindValue(energy)/(1e-20*m*m)+ 
   table->GetComponent(3)->FindValue(energy)/(1e-20*m*m)+ 
   table->GetComponent(4)->FindValue(energy)/(1e-20*m*m) 
 );
 energy*=stpEnergy;
 }
 fclose (myFile);
 abort();
 */
 //
 // end of dump
 //
 // Test of diff XS
 // G4double nrj1 = .26827E+04; // in eV
 // G4double nrj2 =  .57991E+03; // in eV 
 // Shells run from 0 to 4
 // G4cout << DifferentialCrossSection(G4Electron::ElectronDefinition(), nrj1, nrj2, 0)/(1e-20*m*m) << G4endl;
 // G4cout << DifferentialCrossSection(G4Electron::ElectronDefinition(), nrj1, nrj2, 1)/(1e-20*m*m) << G4endl;
 // G4cout << DifferentialCrossSection(G4Electron::ElectronDefinition(), nrj1, nrj2, 2)/(1e-20*m*m) << G4endl;
 // G4cout << DifferentialCrossSection(G4Electron::ElectronDefinition(), nrj1, nrj2, 3)/(1e-20*m*m) << G4endl;
 // G4cout << DifferentialCrossSection(G4Electron::ElectronDefinition(), nrj1, nrj2, 4)/(1e-20*m*m) << G4endl;
 // abort();
 //

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
        G4Exception("G4DNACPA100IonisationModel::RandomSelect","em0002",
                    FatalException,"Model not applicable to particle type.");
    }

    return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergyFromCumulatedDcs
(G4ParticleDefinition* particleDefinition, G4double k, G4int shell)
{
   //G4cout << "*** FAST computation for " << " " << particleDefinition->GetParticleName() << G4endl;

   G4double secondaryElectronKineticEnergy = 0.;
 
   secondaryElectronKineticEnergy= 
   RandomTransferedEnergy(particleDefinition, k/eV, shell)*eV-waterStructure.IonisationEnergy(shell);
 
   //G4cout << RandomTransferedEnergy(particleDefinition, k/eV, shell) << G4endl;
   // SI - 29/03/2014
   if (secondaryElectronKineticEnergy<0.) return 0.;
   //

   return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::RandomTransferedEnergy
(G4ParticleDefinition* particleDefinition,G4double k, G4int ionizationLevelIndex)
{

        G4double random = G4UniformRand();
 
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

        if (particleDefinition == G4Electron::ElectronDefinition())
        {

            // k should be in eV

            std::vector<double>::iterator k2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);

            std::vector<double>::iterator k1 = k2-1;
          
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
            
  if (    random <= eProbaShellMap[ionizationLevelIndex][(*k1)].back() 
          && random <= eProbaShellMap[ionizationLevelIndex][(*k2)].back() )
            
  {
  
  std::vector<double>::iterator prob12 = 
   std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
                     eProbaShellMap[ionizationLevelIndex][(*k1)].end(), random);
                
  std::vector<double>::iterator prob11 = prob12-1;

                
  std::vector<double>::iterator prob22 = 
   std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                     eProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);
                
  std::vector<double>::iterator prob21 = prob22-1;

  valueK1  =*k1;
  valueK2  =*k2;
  valuePROB21 =*prob21;
  valuePROB22 =*prob22;
  valuePROB12 =*prob12;
  valuePROB11 =*prob11;

         
  /*
  G4cout << "        " << random << " " << valuePROB11 << " " 
   << valuePROB12 << " " << valuePROB21 << " " << valuePROB22 << G4endl;
  */

  nrjTransf11 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
  nrjTransf12 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];
  nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
  nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];
         
  /*
  G4cout << "        " << ionizationLevelIndex << " " 
   << random << " " <<valueK1 << " " << valueK2 << G4endl;
         
  G4cout << "        " << random << " " << nrjTransf11 << " " 
   << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;
  */
  
  }
     

  // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)
 
  if ( random > eProbaShellMap[ionizationLevelIndex][(*k1)].back() )

  {
  
  std::vector<double>::iterator prob22 = 
    
    std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
         eProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);
                
  std::vector<double>::iterator prob21 = prob22-1;

  valueK1  =*k1;
  valueK2  =*k2;
  valuePROB21 =*prob21;
  valuePROB22 =*prob22;
         
  //G4cout << "        " << random << " " << valuePROB21 << " " << valuePROB22 << G4endl;
                
  nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
  nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];
         
  G4double interpolatedvalue2 = Interpolate(valuePROB21, valuePROB22, random, nrjTransf21, nrjTransf22);
  
  // zero is explicitely set
  
  G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
  
  /*
  G4cout << "        " << ionizationLevelIndex << " " 
   << random << " " <<valueK1 << " " << valueK2 << G4endl;
         
  G4cout << "        " << random << " " << nrjTransf11 << " " 
   << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;
 
         G4cout << "ici" << " " << value << G4endl;
  */
  
  return value;
  }

 }
      
 //
 
 // End electron case
 
 G4double nrjTransfProduct = nrjTransf11 * nrjTransf12 * nrjTransf21 * nrjTransf22;
        
 //G4cout << "nrjTransfProduct=" << nrjTransfProduct << G4endl;

 if (nrjTransfProduct != 0.)
        {
            nrj = QuadInterpolator(       valuePROB11, valuePROB12,
                                          valuePROB21, valuePROB22,
                                          nrjTransf11, nrjTransf12,
                                          nrjTransf21, nrjTransf22,
                                          valueK1, valueK2,
                                          k, random);
        }
  
 //G4cout << nrj << endl;
  
        return nrj ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNACPA100IonisationModel::RandomizeEjectedElectronEnergyFromCompositionSampling
(G4ParticleDefinition*, G4double tt, G4int shell)
{
 //G4cout << "*** Rejection method for " << " " << particleDefinition->GetParticleName() << G4endl;

 // ***** METHOD 1 ***** (sequential)
 /*

 // ww is KINETIC ENERGY OF SECONDARY ELECTRON
 G4double un=1.;
 G4double deux=2.;
 
 G4double bb = waterStructure.IonisationEnergy(shell);
 G4double uu = waterStructure.UEnergy(shell);
 
 if (tt<=bb) return 0.;
 
 G4double t = tt/bb;
 G4double u = uu/bb;
 G4double tp1 = t + un;
 G4double tu1 = t + u + un;
 G4double tm1 = t - un;
 G4double tp12 = tp1 * tp1;
 G4double dlt = std::log(t);
 
 G4double a1 = t * tm1 / tu1 / tp12;
 G4double a2 = tm1 / tu1 / t / tp1 / deux;
 G4double a3 = dlt * (tp12 - deux * deux ) / tu1 / tp12;
 G4double ato = a1 + a2 + a3;
        
 // 15
 
 G4double r1 =G4UniformRand(); 
 G4double r2 =G4UniformRand();
 G4double r3 =G4UniformRand();

 while (r1<=a1/ato)  
 {
         G4double fx1=r2*tm1/tp1;
         G4double wx1=un/(un-fx1)-un;
         G4double gx1=(t-wx1)/t;
         if(r3 <= gx1) return wx1*bb;
         
         r1 =G4UniformRand(); 
         r2 =G4UniformRand();
         r3 =G4UniformRand();        
 
 }
       
 // 20
 
 while (r1<=(a1+a2)/ato)  
 {
         G4double fx2=tp1+r2*tm1;
         G4double wx2=t-t*tp1/fx2;
         G4double gx2=deux*(un-(t-wx2)/tp1);
         if(r3 <= gx2) return wx2*bb;
      
          // REPEAT 15
          r1 =G4UniformRand(); 
          r2 =G4UniformRand();
          r3 =G4UniformRand();
   
   while (r1<=a1/ato)  
   {
           G4double fx1=r2*tm1/tp1;
           G4double wx1=un/(un-fx1)-un;
           G4double gx1=(t-wx1)/t;
           if(r3 <= gx1) return wx1*bb;    
           r1 =G4UniformRand(); 
           r2 =G4UniformRand();
           r3 =G4UniformRand();       
   }
   // END 15
   
 }
 
 // 30
           
 G4double wx3=std::sqrt(un/(un-r2*(tp12-deux*deux)/tp12))-un;
 G4double gg3=(wx3+un)/(t-wx3);
 G4double gx3=(un+gg3*gg3*gg3)/deux;
       
 while (r3>gx3) 
 {
       
  // 15
 
  r1 =G4UniformRand(); 
  r2 =G4UniformRand();
  r3 =G4UniformRand();

  while (r1<=a1/ato)  
  {
          G4double fx1=r2*tm1/tp1;
          G4double wx1=un/(un-fx1)-un;
          G4double gx1=(t-wx1)/t;
          if(r3 <= gx1) return wx1*bb;
         
          r1 =G4UniformRand(); 
          r2 =G4UniformRand();
          r3 =G4UniformRand();        
 
  }
       
         // 20
 
  while (r1<=(a1+a2)/ato)  
  {
          G4double fx2=tp1+r2*tm1;
          G4double wx2=t-t*tp1/fx2;
          G4double gx2=deux*(un-(t-wx2)/tp1);
          if(r3 <= gx2)return wx2*bb;
      
          // REPEAT 15
          r1 =G4UniformRand(); 
          r2 =G4UniformRand();
          r3 =G4UniformRand();
   
   while (r1<=a1/ato)  
   {
           G4double fx1=r2*tm1/tp1;
           G4double wx1=un/(un-fx1)-un;
           G4double gx1=(t-wx1)/t;
           if(r3 <= gx1) return wx1*bb;
    
           r1 =G4UniformRand(); 
           r2 =G4UniformRand();
           r3 =G4UniformRand();       
   }
   //

  }

  wx3=std::sqrt(un/(un-r2*(tp12-deux*deux)/tp12))-un;
  gg3=(wx3+un)/(t-wx3);
  gx3=(un+gg3*gg3*gg3)/deux;

  }
       
  //
       
  return wx3*bb; 
  */

 // ***** METHOD 2 by M. C. Bordage ***** (optimized)

 G4double un=1.;
 G4double deux=2.;
 
 G4double bb = waterStructure.IonisationEnergy(shell);
 G4double uu = waterStructure.UEnergy(shell);
 
 if (tt<=bb) return 0.;
 
 G4double t = tt/bb;
 G4double u = uu/bb;
 G4double tp1 = t + un;
 G4double tu1 = t + u + un;
 G4double tm1 = t - un;
 G4double tp12 = tp1 * tp1;
 G4double dlt = std::log(t);
 
 G4double a1 = t * tm1 / tu1 / tp12;
 G4double a2 = tm1 / tu1 / t / tp1 / deux;
 G4double a3 = dlt * (tp12 - deux * deux ) / tu1 / tp12;
 G4double ato = a1 + a2 + a3;

 G4double A1 = a1/ato;
 G4double A2 = (a1+a2)/ato;
 G4int F = 0;
 G4double fx=0; 
 G4double gx=0;
 G4double gg=0;
 G4double wx=0;
 
 G4double r1=0;
 G4double r2=0;
 G4double r3=0;
 
 //
 
 do 
 {
   r1 =G4UniformRand(); 
   r2 =G4UniformRand();
   r3 =G4UniformRand();

   if (r1>A2)
    F=3;
   else if ((r1>A1) && (r1< A2))
    F=2;
   else
    F=1;

   switch (F)
   {   
    case 1:
    {
     fx=r2*tm1/tp1;
     wx=un/(un-fx)-un;
     gx=(t-wx)/t; 
     break;
    }      
  
    case 2:
    {
     fx=tp1+r2*tm1;
     wx=t-t*tp1/fx;
     gx=deux*(un-(t-wx)/tp1); 
     break;
    } 
  
    case 3:
    {
     fx=un-r2*(tp12-deux*deux)/tp12;
     wx=sqrt(un/fx)-un;
     gg=(wx+un)/(t-wx);
     gx=(un+gg*gg*gg)/deux;
     break;
    } 
   } // switch
   
  } while (r3>gx);
       
  return wx*bb; 

}

