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
// $Id$
//

#include "G4DNABornIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel::G4DNABornIonisationModel(const G4ParticleDefinition*,
                                                   const G4String& nam)
    :G4VEmModel(nam),isInitialised(false)
{
    //  nistwater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 )
    {
        G4cout << "Born ionisation model is constructed " << G4endl;
    }

    //Mark this model as "applicable" for atomic deexcitation
    SetDeexcitationFlag(true);
    fAtomDeexcitation = 0;
    fParticleChangeForGamma = 0;
    fpMolWaterDensity = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel::~G4DNABornIonisationModel()
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
    pVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& /*cuts*/)
{

    if (verboseLevel > 3)
        G4cout << "Calling G4DNABornIonisationModel::Initialise()" << G4endl;

    // Energy limits

    G4String fileElectron("dna/sigma_ionisation_e_born");
    G4String fileProton("dna/sigma_ionisation_p_born");

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
    G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

    G4String electron;
    G4String proton;

    G4double scaleFactor = (1.e-22 / 3.343) * m*m;

    char *path = getenv("G4LEDATA");

    // *** ELECTRON

    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    lowEnergyLimit[electron] = 11. * eV;
    highEnergyLimit[electron] = 1. * MeV;

    // Cross section
    
    G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);

    tableData[electron] = tableE;
    
    // Final state
    
    std::ostringstream eFullFileName;
    eFullFileName << path << "/dna/sigmadiff_ionisation_e_born.dat";
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

    if (!eDiffCrossSection)
    {
        G4Exception("G4DNABornIonisationModel::Initialise","em0003",
                    FatalException,"Missing data file:/dna/sigmadiff_ionisation_e_born.dat");
    }

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

            // SI - only if eof is not reached !
            if (!eDiffCrossSection.eof()) eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

            eVecm[tDummy].push_back(eDummy);

        }
    }

    // *** PROTON

    proton = protonDef->GetParticleName();

    tableFile[proton] = fileProton;

    lowEnergyLimit[proton] = 500. * keV;
    highEnergyLimit[proton] = 100. * MeV;

    // Cross section
    
    G4DNACrossSectionDataSet* tableP = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableP->LoadData(fileProton);

    tableData[proton] = tableP;
    
    // Final state

    std::ostringstream pFullFileName;
    pFullFileName << path << "/dna/sigmadiff_ionisation_p_born.dat";
    std::ifstream pDiffCrossSection(pFullFileName.str().c_str());
    
    if (!pDiffCrossSection)
    {
        G4Exception("G4DNABornIonisationModel::Initialise","em0003",
                    FatalException,"Missing data file:/dna/sigmadiff_ionisation_p_born.dat");
    }

    pTdummyVec.push_back(0.);
    while(!pDiffCrossSection.eof())
    {
        double tDummy;
        double eDummy;
        pDiffCrossSection>>tDummy>>eDummy;
        if (tDummy != pTdummyVec.back()) pTdummyVec.push_back(tDummy);
        for (int j=0; j<5; j++)
        {
            pDiffCrossSection>>pDiffCrossSectionData[j][tDummy][eDummy];

            // SI - only if eof is not reached !
            if (!pDiffCrossSection.eof()) pDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

            pVecm[tDummy].push_back(eDummy);
        }
    }

    //

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
        G4cout << "Born ionisation model is initialized " << G4endl
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

G4double G4DNABornIonisationModel::CrossSectionPerVolume(const G4Material* material,
                                                         const G4ParticleDefinition* particleDefinition,
                                                         G4double ekin,
                                                         G4double,
                                                         G4double)
{
    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNABornIonisationModel" << G4endl;

    if (
            particleDefinition != G4Proton::ProtonDefinition()
            &&
            particleDefinition != G4Electron::ElectronDefinition()
            )

        return 0;

    // Calculate total cross section for model

    G4double lowLim = 0;
    G4double highLim = 0;
    G4double sigma=0;

    G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0)
        //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
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

        if (ekin >= lowLim && ekin < highLim)
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
                G4Exception("G4DNABornIonisationModel::CrossSectionPerVolume","em0002",
                            FatalException,"Model not applicable to particle type.");
            }
        }

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "°°° G4DNABornIonisationModel - XS INFO START" << G4endl;
            G4cout << "°°° Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
            G4cout << "°°° Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "°°° Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
            //      G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
            G4cout << "°°° G4DNABornIonisationModel - XS INFO END" << G4endl;
        }

    } // if (waterMaterial)

    // return sigma*material->GetAtomicNumDensityVector()[1];
    return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                 const G4MaterialCutsCouple* ,//must be set!
                                                 const G4DynamicParticle* particle,
                                                 G4double,
                                                 G4double)
{

    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNABornIonisationModel" << G4endl;

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

        G4int ionizationShell = RandomSelect(k,particleName);

        // sample deexcitation
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

            //	FOR DEBUG ONLY
            //	if (ionizationShell == 4) {
            //
            //	  G4cout << "Z: " << Z << " as: " << as
            //               << " ionizationShell: " << ionizationShell << " bindingEnergy: "<< bindingEnergy/eV << G4endl;
            //        G4cout << "Press <Enter> key to continue..." << G4endl;
            //	  G4cin.ignore();
            //	}

            const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
            secNumberInit = fvect->size();
            fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
            secNumberFinal = fvect->size();
        }

        G4double secondaryKinetic = RandomizeEjectedElectronEnergy(particle->GetDefinition(),k,ionizationShell);

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

        // note thta secondaryKinetic is the nergy of the delta ray, not of all secondaries.
        G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;
        G4double deexSecEnergy = 0;
        for (G4int j=secNumberInit; j < secNumberFinal; j++) {

            deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();

        }

        fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy-secondaryKinetic-deexSecEnergy);

        G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
        fvect->push_back(dp);


        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
                                                               ionizationShell,
                                                               theIncomingTrack);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
                                                                  G4double k, G4int shell)
{
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
        G4int nEnergySteps = 50;

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
        //

        G4double secondaryElectronKineticEnergy=0.;
        do
        {
            secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer-waterStructure.IonisationEnergy(shell));
        } while(G4UniformRand()*crossSectionMaximum >
                DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+waterStructure.IonisationEnergy(shell))/eV,shell));

        return secondaryElectronKineticEnergy;

    }

    else if (particleDefinition == G4Proton::ProtonDefinition())
    {
        G4double maximumKineticEnergyTransfer = 4.* (electron_mass_c2 / proton_mass_c2) * k;

        G4double crossSectionMaximum = 0.;
        for (G4double value = waterStructure.IonisationEnergy(shell);
             value<=4.*waterStructure.IonisationEnergy(shell) ;
             value+=0.1*eV)
        {
            G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
            if (differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
        }

        G4double secondaryElectronKineticEnergy = 0.;
        do
        {
            secondaryElectronKineticEnergy = G4UniformRand() * maximumKineticEnergyTransfer;
        } while(G4UniformRand()*crossSectionMaximum >=
                DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+waterStructure.IonisationEnergy(shell))/eV,shell));

        return secondaryElectronKineticEnergy;
    }

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNABornIonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
                                                                 G4double k,
                                                                 G4double secKinetic,
                                                                 G4double & cosTheta,
                                                                 G4double & phi )
{
    if (particleDefinition == G4Electron::ElectronDefinition())
    {
        phi = twopi * G4UniformRand();
        if (secKinetic < 50.*eV) cosTheta = (2.*G4UniformRand())-1.;
        else if (secKinetic <= 200.*eV)
        {
            if (G4UniformRand() <= 0.1) cosTheta = (2.*G4UniformRand())-1.;
            else cosTheta = G4UniformRand()*(std::sqrt(2.)/2);
        }
        else
        {
            G4double sin2O = (1.-secKinetic/k) / (1.+secKinetic/(2.*electron_mass_c2));
            cosTheta = std::sqrt(1.-sin2O);
        }
    }

    else if (particleDefinition == G4Proton::ProtonDefinition())
    {
        G4double maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
        phi = twopi * G4UniformRand();

        // cosTheta = std::sqrt(secKinetic / maxSecKinetic);

        // Restriction below 100 eV from Emfietzoglou (2000)

        if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
        else cosTheta = (2.*G4UniformRand())-1.;
        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double G4DNABornIonisationModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition, 
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

        if (particleDefinition == G4Proton::ProtonDefinition())
        {
            // k should be in eV and energy transfer eV also
            std::vector<double>::iterator t2 = std::upper_bound(pTdummyVec.begin(),pTdummyVec.end(), k);
            std::vector<double>::iterator t1 = t2-1;

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

            xs11 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE11];
            xs12 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE12];
            xs21 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE21];
            xs22 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE22];

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

G4double G4DNABornIonisationModel::LogLogInterpolate(G4double e1, 
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

G4double G4DNABornIonisationModel::QuadInterpolator(G4double e11, G4double e12, 
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

G4int G4DNABornIonisationModel::RandomSelect(G4double k, const G4String& particle )
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
        G4Exception("G4DNABornIonisationModel::RandomSelect","em0002",
                    FatalException,"Model not applicable to particle type.");
    }

    return level;
}

