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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#include "G4DNAPTBIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"

G4DNAPTBIonisationModel::G4DNAPTBIonisationModel(const G4String& applyToMaterial,
                                                 const G4ParticleDefinition*,
                                                 const G4String& nam, const G4bool isAuger)
    : G4VDNAModel(nam, applyToMaterial)
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
        G4cout << "PTB ionisation model is constructed " << G4endl;
    }

    if(isAuger)
    {
        // create the PTB Auger model
        fDNAPTBAugerModel = new G4DNAPTBAugerModel("e-_G4DNAPTBAugerModel");
    }
    else
    {
        // no PTB Auger model
        fDNAPTBAugerModel = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAPTBIonisationModel::~G4DNAPTBIonisationModel()
{
    // To delete the DNAPTBAugerModel created at initialisation of the ionisation class
    if(fDNAPTBAugerModel) delete fDNAPTBAugerModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                         const G4DataVector& /*cuts*/, G4ParticleChangeForGamma*)
{

    if (verboseLevel > 3)
        G4cout << "Calling G4DNAPTBIonisationModel::Initialise()" << G4endl;

    G4double scaleFactor = 1e-16 * cm*cm;
    G4double scaleFactorBorn = (1.e-22 / 3.343) * m*m;

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
    G4ParticleDefinition* protonDef= G4Proton::ProtonDefinition();

    //*******************************************************
    // Cross section data
    //*******************************************************

    if(particle == electronDef)
    {
        G4String particleName = particle->GetParticleName();

        // Raw materials
        //
        AddCrossSectionData("THF",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_THF",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_THF",
                            scaleFactor);
        SetLowELimit("THF", particleName, 12.*eV);
        SetHighELimit("THF", particleName, 1.*keV);

        AddCrossSectionData("PY",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PY",
                            scaleFactor);
        SetLowELimit("PY", particleName, 12.*eV);
        SetHighELimit("PY", particleName, 1.*keV);

        AddCrossSectionData("PU",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PU",
                            scaleFactor);
        SetLowELimit("PU", particleName, 12.*eV);
        SetHighELimit("PU", particleName, 1.*keV);

        AddCrossSectionData("TMP",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_TMP",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_TMP",
                            scaleFactor);
        SetLowELimit("TMP", particleName, 12.*eV);
        SetHighELimit("TMP", particleName, 1.*keV);

        AddCrossSectionData("G4_WATER",
                            particleName,
                            "dna/sigma_ionisation_e_born",
                            "dna/sigmadiff_ionisation_e_born",
                            scaleFactorBorn);
        SetLowELimit("G4_WATER", particleName, 12.*eV);
        SetHighELimit("G4_WATER", particleName, 1.*keV);

        // DNA materials
        //
        AddCrossSectionData("backbone_THF",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_THF",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_THF",
                            scaleFactor*33./30);
        SetLowELimit("backbone_THF", particleName, 12.*eV);
        SetHighELimit("backbone_THF", particleName, 1.*keV);

        AddCrossSectionData("cytosine_PY",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PY",
                            scaleFactor*42./30);
        SetLowELimit("cytosine_PY", particleName, 12.*eV);
        SetHighELimit("cytosine_PY", particleName, 1.*keV);

        AddCrossSectionData("thymine_PY",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PY",
                            scaleFactor*48./30);
        SetLowELimit("thymine_PY", particleName, 12.*eV);
        SetHighELimit("thymine_PY", particleName, 1.*keV);

        AddCrossSectionData("adenine_PU",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PU",
                            scaleFactor*50./44);
        SetLowELimit("adenine_PU", particleName, 12.*eV);
        SetHighELimit("adenine_PU", particleName, 1.*keV);

        AddCrossSectionData("guanine_PU",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_PU",
                            scaleFactor*56./44);
        SetLowELimit("guanine_PU", particleName, 12.*eV);
        SetHighELimit("guanine_PU", particleName, 1.*keV);

        AddCrossSectionData("backbone_TMP",
                            particleName,
                            "dna/sigma_ionisation_e-_PTB_TMP",
                            "dna/sigmadiff_cumulated_ionisation_e-_PTB_TMP",
                            scaleFactor*33./50);
        SetLowELimit("backbone_TMP", particleName, 12.*eV);
        SetHighELimit("backbone_TMP", particleName, 1.*keV);
    }

    else if (particle == protonDef)
    {
        G4String particleName = particle->GetParticleName();

        // Raw materials
        //
        AddCrossSectionData("THF",
                            particleName,
                            "dna/sigma_ionisation_p_HKS_THF",
                            "dna/sigmadiff_cumulated_ionisation_p_PTB_THF",
                            scaleFactor);
        SetLowELimit("THF", particleName, 70.*keV);
        SetHighELimit("THF", particleName, 10.*MeV);


        AddCrossSectionData("PY",
                            particleName,
                            "dna/sigma_ionisation_p_HKS_PY",
                            "dna/sigmadiff_cumulated_ionisation_p_PTB_PY",
                            scaleFactor);
        SetLowELimit("PY", particleName, 70.*keV);
        SetHighELimit("PY", particleName, 10.*MeV);

        /*
         AddCrossSectionData("PU",
                                    particleName,
                                    "dna/sigma_ionisation_e-_PTB_PU",
                                    "dna/sigmadiff_cumulated_ionisation_e-_PTB_PU",
                                    scaleFactor);
                SetLowELimit("PU", particleName2, 70.*keV);
                SetHighELimit("PU", particleName2, 10.*keV);
*/

        AddCrossSectionData("TMP",
                            particleName,
                            "dna/sigma_ionisation_p_HKS_TMP",
                            "dna/sigmadiff_cumulated_ionisation_p_PTB_TMP",
                            scaleFactor);
        SetLowELimit("TMP", particleName, 70.*keV);
        SetHighELimit("TMP", particleName, 10.*MeV);
    }

    // *******************************************************
    // deal with composite materials
    // *******************************************************

    LoadCrossSectionData(particle->GetParticleName() );

    // *******************************************************
    // Verbose
    // *******************************************************

    // initialise DNAPTBAugerModel
    if(fDNAPTBAugerModel) fDNAPTBAugerModel->Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBIonisationModel::CrossSectionPerVolume(const G4Material* /*material*/,
                                                        const G4String& materialName,
                                                        const G4ParticleDefinition* p,
                                                        G4double ekin,
                                                        G4double /*emin*/,
                                                        G4double /*emax*/)
{
    if(verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNAPTBIonisationModel" << G4endl;

    // initialise the cross section value (output value)
    G4double sigma(0);

    // Get the current particle name
    const G4String& particleName = p->GetParticleName();

    // Set the low and high energy limits
    G4double lowLim = GetLowELimit(materialName, particleName);
    G4double highLim = GetHighELimit(materialName, particleName);

    // Check that we are in the correct energy range
    if (ekin >= lowLim && ekin < highLim)
    {
        // Get the map with all the model data tables
        TableMapData* tableData = GetTableData();

        // Retrieve the cross section value for the current material, particle and energy values
        sigma = (*tableData)[materialName][particleName]->FindValue(ekin);

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "°°° G4DNAPTBIonisationModel - XS INFO START" << G4endl;
            G4cout << "°°° Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
            G4cout << "°°° Cross section per "<< materialName <<" molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "°°° G4DNAPTBIonisationModel - XS INFO END" << G4endl;
        }
    }

    // Return the cross section value
    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                const G4MaterialCutsCouple* /*couple*/,
                                                const G4String& materialName,
                                                const G4DynamicParticle* aDynamicParticle,
                                                G4ParticleChangeForGamma* particleChangeForGamma,
                                                G4double /*tmin*/,
                                                G4double /*tmax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNAPTBIonisationModel" << G4endl;

    // Get the current particle energy
    G4double k = aDynamicParticle->GetKineticEnergy();

    // Get the current particle name
    const G4String& particleName = aDynamicParticle->GetDefinition()->GetParticleName();

    // Get the energy limits
    G4double lowLim = GetLowELimit(materialName, particleName);
    G4double highLim = GetHighELimit(materialName, particleName);

    // Check if we are in the correct energy range
    if (k >= lowLim && k < highLim)
    {
        G4ParticleMomentum primaryDirection = aDynamicParticle->GetMomentumDirection();
        G4double particleMass = aDynamicParticle->GetDefinition()->GetPDGMass();
        G4double totalEnergy = k + particleMass;
        G4double pSquare = k * (totalEnergy + particleMass);
        G4double totalMomentum = std::sqrt(pSquare);

        // Get the ionisation shell from a random sampling
        G4int ionizationShell = RandomSelectShell(k, particleName, materialName);

        // Get the binding energy from the ptbStructure class
        G4double bindingEnergy = ptbStructure.IonisationEnergy(ionizationShell, materialName);

        // Initialize the secondary kinetic energy to a negative value.
        G4double secondaryKinetic (-1000*eV);

        if(materialName!="G4_WATER")
        {
            // Get the energy of the secondary particle
            secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulated(aDynamicParticle->GetDefinition(),k/eV,ionizationShell, materialName);
        }
        else
        {
            secondaryKinetic = RandomizeEjectedElectronEnergy(aDynamicParticle->GetDefinition(),k,ionizationShell, materialName);
        }

        if(secondaryKinetic<=0)
        {
            G4cout<<"Fatal error *************************************** "<<secondaryKinetic/eV<<G4endl;
            G4cout<<"secondaryKinetic: "<<secondaryKinetic/eV<<G4endl;
            G4cout<<"k: "<<k/eV<<G4endl;
            G4cout<<"shell: "<<ionizationShell<<G4endl;
            G4cout<<"material:"<<materialName<<G4endl;
            exit(EXIT_FAILURE);
        }

        G4double cosTheta = 0.;
        G4double phi = 0.;
        RandomizeEjectedElectronDirection(aDynamicParticle->GetDefinition(), k, secondaryKinetic, cosTheta, phi);

        G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
        G4double dirX = sinTheta*std::cos(phi);
        G4double dirY = sinTheta*std::sin(phi);
        G4double dirZ = cosTheta;
        G4ThreeVector deltaDirection(dirX,dirY,dirZ);
        deltaDirection.rotateUz(primaryDirection);

        // The model is written only for electron  and thus we want the change the direction of the incident electron
        // after each ionization. However, if other particle are going to be introduced within this model the following should be added:
        //
        // Check if the particle is an electron
        if(aDynamicParticle->GetDefinition() == G4Electron::ElectronDefinition() )
        {
            // If yes do the following code until next commented "else" statement

            G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));
            G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
            G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
            G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
            G4double finalMomentum = std::sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz);
            finalPx /= finalMomentum;
            finalPy /= finalMomentum;
            finalPz /= finalMomentum;

            G4ThreeVector direction(finalPx,finalPy,finalPz);
            if(direction.unit().getX()>1||direction.unit().getY()>1||direction.unit().getZ()>1)
            {
                G4cout<<"Fatal error ****************************"<<G4endl;
                G4cout<<"direction problem "<<direction.unit()<<G4endl;
                exit(EXIT_FAILURE);
            }

            // Give the new direction to the particle
            particleChangeForGamma->ProposeMomentumDirection(direction.unit()) ;
        }
        // If the particle is not an electron
        else particleChangeForGamma->ProposeMomentumDirection(primaryDirection) ;

        // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
        G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;

        if(scatteredEnergy<=0)
        {
            G4cout<<"Fatal error ****************************"<<G4endl;
            G4cout<<"k: "<<k/eV<<G4endl;
            G4cout<<"secondaryKinetic: "<<secondaryKinetic/eV<<G4endl;
            G4cout<<"shell: "<<ionizationShell<<G4endl;
            G4cout<<"bindingEnergy: "<<bindingEnergy/eV<<G4endl;
            G4cout<<"scatteredEnergy: "<<scatteredEnergy/eV<<G4endl;
            G4cout<<"material: "<<materialName<<G4endl;
            exit(EXIT_FAILURE);
        }

        // Set the new energy of the particle
        particleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);

        // Set the energy deposited by the ionization
        particleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy-secondaryKinetic);

        // Create the new particle with its characteristics
        G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
        fvect->push_back(dp);

        // Check if the auger model is activated (ie instanciated)
        if(fDNAPTBAugerModel)
        {
            // run the PTB Auger model
            if(materialName!="G4_WATER") fDNAPTBAugerModel->ComputeAugerEffect(fvect, materialName, bindingEnergy);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBIonisationModel::ReadDiffCSFile(const G4String& materialName,
                                             const G4String& particleName,
                                             const G4String& file,
                                             const G4double scaleFactor)
{
    // To read and save the informations contained within the differential cross section files

    // get the path of the G4LEDATA data folder
    char *path = getenv("G4LEDATA");
    // if it is not found then quit and print error message
    if(!path)
    {
        G4Exception("G4DNAPTBIonisationModel::ReadAllDiffCSFiles","em0006",
                    FatalException,"G4LEDATA environment variable not set.");
        return;
    }

    // build the fullFileName path of the data file
    std::ostringstream fullFileName;
    fullFileName << path <<"/"<< file<<".dat";

    // open the data file
    std::ifstream diffCrossSection (fullFileName.str().c_str());
    // error if file is not there
    std::stringstream endPath;
    if (!diffCrossSection)
    {
        endPath << "Missing data file: "<<file;
        G4Exception("G4DNAPTBIonisationModel::Initialise","em0003",
                    FatalException, endPath.str().c_str());
    }

    // load data from the file
    fTMapWithVec[materialName][particleName].push_back(0.);

    G4String line;

    // read the file until we reach the end of file point
    // fill fTMapWithVec, diffCrossSectionData, fEnergyTransferData, fProbaShellMap and fEMapWithVector
    while(std::getline(diffCrossSection, line))
    {
        // check if the line is comment or empty
        //
        std::istringstream testIss(line);
        G4String test;
        testIss >> test;
        // check first caracter to determine if following information is data or comments
        if(test=="#")
        {
            // skip the line by beginning a new while loop.
            continue;
        }
        // check if line is empty
        else if(line.empty())
        {
            // skip the line by beginning a new while loop.
            continue;
        }
        //
        // end of the check

        // transform the line into a iss
        std::istringstream iss(line);

        // Initialise the variables to be filled
        double T;
        double E;

        // Filled T and E with the first two numbers of each file line
        iss>>T>>E;

        // Fill the fTMapWithVec container with all the different T values contained within the file.
        // Duplicate must be avoided and this is the purpose of the if statement
        if (T != fTMapWithVec[materialName][particleName].back()) fTMapWithVec[materialName][particleName].push_back(T);

        // iterate on each shell of the corresponding material
        for (int shell=0, eshell=ptbStructure.NumberOfLevels(materialName); shell<eshell; ++shell)
        {
            // map[material][particle][shell][T][E]=diffCrossSectionValue
            // Fill the map with the informations of the input file
            iss>>diffCrossSectionData[materialName][particleName][shell][T][E];

            if(materialName!="G4_WATER")
            {
                // map[material][particle][shell][T][CS]=E
                // Fill the map
                fEnergySecondaryData[materialName][particleName][shell][T][diffCrossSectionData[materialName][particleName][shell][T][E] ]=E;

                // map[material][particle][shell][T]=CS_vector
                // Fill the vector within the map
                fProbaShellMap[materialName][particleName][shell][T].push_back(diffCrossSectionData[materialName][particleName][shell][T][E]);
            }
            else
            {
                diffCrossSectionData[materialName][particleName][shell][T][E]*=scaleFactor;

                fEMapWithVector[materialName][particleName][T].push_back(E);
            }
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBIonisationModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition,
                                                                 G4double k, G4int shell, const G4String& materialName)
{
    if (particleDefinition == G4Electron::ElectronDefinition())
    {
        //G4double Tcut=25.0E-6;
        G4double maximumEnergyTransfer=0.;
        if ((k+ptbStructure.IonisationEnergy(shell, materialName))/2. > k) maximumEnergyTransfer=k;
        else maximumEnergyTransfer = (k+ptbStructure.IonisationEnergy(shell,materialName))/2.;

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

        //if (k > Tcut)
        //{
        G4double crossSectionMaximum = 0.;

        G4double minEnergy = ptbStructure.IonisationEnergy(shell, materialName);
        G4double maxEnergy = maximumEnergyTransfer;
        G4int nEnergySteps = 50;
        G4double value(minEnergy);
        G4double stpEnergy(std::pow(maxEnergy/value, 1./static_cast<G4double>(nEnergySteps-1)));
        G4int step(nEnergySteps);
        while (step>0)
        {
            step--;
            G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell, materialName);
            if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
            value *= stpEnergy;

        }
        //


        G4double secondaryElectronKineticEnergy=0.;

        do
        {
            secondaryElectronKineticEnergy = G4UniformRand() * (maximumEnergyTransfer-ptbStructure.IonisationEnergy(shell, materialName));

        } while(G4UniformRand()*crossSectionMaximum >
                DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+ptbStructure.IonisationEnergy(shell, materialName))/eV,shell, materialName));

        return secondaryElectronKineticEnergy;

        //        }

        //        else if (k < Tcut)
        //        {

        //            G4double bindingEnergy = ptbStructure.IonisationEnergy(shell, materialName);
        //            G4double maxEnergy = ((k-bindingEnergy)/2.);

        //            G4double secondaryElectronKineticEnergy = G4UniformRand()*maxEnergy;
        //            return secondaryElectronKineticEnergy;
        //        }
    }


    else if (particleDefinition == G4Proton::ProtonDefinition())
    {
        G4double maximumKineticEnergyTransfer = 4.* (electron_mass_c2 / proton_mass_c2) * k;

        G4double crossSectionMaximum = 0.;
        for (G4double value = ptbStructure.IonisationEnergy(shell, materialName);
             value<=4.*ptbStructure.IonisationEnergy(shell, materialName) ;
             value+=0.1*eV)
        {
            G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell, materialName);
            if (differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
        }

        G4double secondaryElectronKineticEnergy = 0.;
        do
        {
            secondaryElectronKineticEnergy = G4UniformRand() * maximumKineticEnergyTransfer;
        } while(G4UniformRand()*crossSectionMaximum >=
                DifferentialCrossSection(particleDefinition, k/eV,(secondaryElectronKineticEnergy+ptbStructure.IonisationEnergy(shell, materialName))/eV,shell, materialName));

        return secondaryElectronKineticEnergy;
    }

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAPTBIonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition,
                                                                G4double k,
                                                                G4double secKinetic,
                                                                G4double & cosTheta,
                                                                G4double & phi)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

double G4DNAPTBIonisationModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition,
                                                         G4double k,
                                                         G4double energyTransfer,
                                                         G4int ionizationLevelIndex,
                                                         const G4String& materialName)
{
    G4double sigma = 0.;
    const G4String& particleName = particleDefinition->GetParticleName();

    G4double shellEnergy (ptbStructure.IonisationEnergy(ionizationLevelIndex, materialName));
    G4double kSE (energyTransfer-shellEnergy);

    if (energyTransfer >= shellEnergy)
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
            std::vector<double>::iterator t2 = std::upper_bound(fTMapWithVec[materialName][particleName].begin(),fTMapWithVec[materialName][particleName].end(), k);
            std::vector<double>::iterator t1 = t2-1;

            // SI : the following condition avoids situations where energyTransfer >last vector element
            if (kSE <= fEMapWithVector[materialName][particleName][(*t1)].back() && kSE <= fEMapWithVector[materialName][particleName][(*t2)].back() )
            {
                std::vector<double>::iterator e12 = std::upper_bound(fEMapWithVector[materialName][particleName][(*t1)].begin(),fEMapWithVector[materialName][particleName][(*t1)].end(), kSE);
                std::vector<double>::iterator e11 = e12-1;

                std::vector<double>::iterator e22 = std::upper_bound(fEMapWithVector[materialName][particleName][(*t2)].begin(),fEMapWithVector[materialName][particleName][(*t2)].end(), kSE);
                std::vector<double>::iterator e21 = e22-1;

                valueT1  =*t1;
                valueT2  =*t2;
                valueE21 =*e21;
                valueE22 =*e22;
                valueE12 =*e12;
                valueE11 =*e11;

                xs11 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT1][valueE11];
                xs12 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT1][valueE12];
                xs21 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT2][valueE21];
                xs22 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT2][valueE22];
            }
        }

        if (particleDefinition == G4Proton::ProtonDefinition())
        {
            // k should be in eV and energy transfer eV also
            std::vector<double>::iterator t2 = std::upper_bound(fTMapWithVec[materialName][particleName].begin(),fTMapWithVec[materialName][particleName].end(), k);
            std::vector<double>::iterator t1 = t2-1;

            std::vector<double>::iterator e12 = std::upper_bound(fEMapWithVector[materialName][particleName][(*t1)].begin(),fEMapWithVector[materialName][particleName][(*t1)].end(), kSE);
            std::vector<double>::iterator e11 = e12-1;

            std::vector<double>::iterator e22 = std::upper_bound(fEMapWithVector[materialName][particleName][(*t2)].begin(),fEMapWithVector[materialName][particleName][(*t2)].end(), kSE);
            std::vector<double>::iterator e21 = e22-1;

            valueT1  =*t1;
            valueT2  =*t2;
            valueE21 =*e21;
            valueE22 =*e22;
            valueE12 =*e12;
            valueE11 =*e11;

            xs11 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT1][valueE11];
            xs12 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT1][valueE12];
            xs21 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT2][valueE21];
            xs22 = diffCrossSectionData[materialName][particleName][ionizationLevelIndex][valueT2][valueE22];
        }

        G4double xsProduct = xs11 * xs12 * xs21 * xs22;

        if (xsProduct != 0.)
        {
            sigma = QuadInterpolator(valueE11, valueE12,
                                     valueE21, valueE22,
                                     xs11, xs12,
                                     xs21, xs22,
                                     valueT1, valueT2,
                                     k, kSE);
        }
    }


    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBIonisationModel::RandomizeEjectedElectronEnergyFromCumulated(G4ParticleDefinition* particleDefinition,G4double k, G4int ionizationLevelIndex, const G4String& materialName)
{
    // k should be in eV

    // Schematic explanation.
    // We will do an interpolation to get a final E value (ejected electron energy).
    // 1/ We choose a random number between 0 and 1 (ie we select a cumulated cross section).
    // 2/ We look for T_lower and T_upper.
    // 3/ We look for the cumulated corresponding cross sections and their associated E values.
    //
    // T_low | CS_low_1 -> E_low_1
    //       | CS_low_2 -> E_low_2
    // T_up  | CS_up_1 -> E_up_1
    //       | CS_up_2 -> E_up_2
    //
    // 4/ We interpolate to get our E value.
    //
    // T_low | CS_low_1 -> E_low_1 -----
    //       |                          |----> E_low --
    //       | CS_low_2 -> E_low_2 -----               |
    //                                                 | ---> E_final
    // T_up  | CS_up_1 -> E_up_1 -------               |
    //       |                          |----> E_up ---
    //       | CS_up_2 -> E_up_2 -------

    // Initialize some values
    //
    G4double ejectedElectronEnergy = 0.;
    G4double valueK1 = 0;
    G4double valueK2 = 0;
    G4double valueCumulCS21 = 0;
    G4double valueCumulCS22 = 0;
    G4double valueCumulCS12 = 0;
    G4double valueCumulCS11 = 0;
    G4double secElecE11 = 0;
    G4double secElecE12 = 0;
    G4double secElecE21 = 0;
    G4double secElecE22 = 0;
    G4String particleName = particleDefinition->GetParticleName();

    // ***************************************************************************
    // Get a random number between 0 and 1 to compare with the cumulated CS
    // ***************************************************************************
    //
    // It will allow us to choose an ejected electron energy with respect to the CS.
    G4double random = G4UniformRand();

    // **********************************************
    // Take the input from the data tables
    // **********************************************

    // Cumulated tables are like this: T E cumulatedCS1 cumulatedCS2 cumulatedCS3
    // We have two sets of loaded data: fTMapWithVec which contains data about T (incident particle energy)
    // and fProbaShellMap which contains cumulated cross section data.
    // Since we already have a specific T energy value which could not be explicitly in the table, we must interpolate all the values.

    // First, we select the upper and lower T data values surrounding our T value (ie "k").
    std::vector<double>::iterator k2 = std::upper_bound(fTMapWithVec[materialName][particleName].begin(),fTMapWithVec[materialName][particleName].end(), k);
    std::vector<double>::iterator k1 = k2-1;

    // Check if we have found a k2 value (0 if we did not found it).
    // A missing k2 value can be caused by a energy to high for the data table,
    // Ex : table done for 12*eV -> 1000*eV and k=2000*eV
    // then k2 = 0 and k1 = max of the table.
    // To detect this, we check that k1 is not superior to k2.
    if(*k1 > *k2)
    {
        // Error
        G4cerr<<"**************** Fatal error ******************"<<G4endl;
        G4cerr<<"G4DNAPTBIonisationModel::RandomizeEjectedElectronEnergyFromCumulated"<<G4endl;
        G4cerr<<"You have *k1 > *k2 with k1 "<<*k1<<" and k2 "<<*k2<<G4endl;
        G4cerr<<"This may be because the energy of the incident particle is to high for the data table."<<G4endl;
        G4cerr<<"Particle energy (eV): "<<k<<G4endl;
        exit(EXIT_FAILURE);
    }


    // We have a random number and we select the cumulated cross section data values surrounding our random number.
    // But we need to do that for each T value (ie two T values) previously selected.
    //
    // First one.
    std::vector<double>::iterator cumulCS12 = std::upper_bound(fProbaShellMap[materialName][particleName][ionizationLevelIndex][(*k1)].begin(),
            fProbaShellMap[materialName][particleName][ionizationLevelIndex][(*k1)].end(), random);
    std::vector<double>::iterator cumulCS11 = cumulCS12-1;
    // Second one.
    std::vector<double>::iterator cumulCS22 = std::upper_bound(fProbaShellMap[materialName][particleName][ionizationLevelIndex][(*k2)].begin(),
            fProbaShellMap[materialName][particleName][ionizationLevelIndex][(*k2)].end(), random);
    std::vector<double>::iterator cumulCS21 = cumulCS22-1;

    // Now that we have the "values" through pointers, we access them.
    valueK1  = *k1;
    valueK2  = *k2;
    valueCumulCS11 = *cumulCS11;
    valueCumulCS12 = *cumulCS12;
    valueCumulCS21 = *cumulCS21;
    valueCumulCS22 = *cumulCS22;

    // *************************************************************
    // Do the interpolation to get the ejected electron energy
    // *************************************************************

    // Here we will get four E values corresponding to our four cumulated cross section values previously selected.
    // But we need to take into account a specific case: we have selected a shell by using the ionisation cross section table
    // and, since we get two T values, we could have differential cross sections (or cumulated) equal to 0 for the lower T
    // and not for the upper T. When looking for the cumulated cross section values which surround the selected random number (for the lower T),
    // the upper_bound method will only found 0 values. Thus, the upper_bound method will return the last E value present in the table for the
    // selected T. The last E value being the highest, we will later perform an interpolation between a high E value (for the lower T) and
    // a small E value (for the upper T). This is inconsistent because if the cross section are equal to zero for the lower T then it
    // means it is not possible to ionize and, thus, to have a secondary electron. But, in our situation, it is possible to ionize for the upper T
    // AND for an interpolate T value between Tupper Tlower. That's why the final E value should be interpolate between 0 and the E value (upper T).
    //
    if(cumulCS12==fProbaShellMap[materialName][particleName][ionizationLevelIndex][(*k1)].end())
    {
        // Here we are in the special case and we force Elower1 and Elower2 to be equal at 0 for the interpolation.
        secElecE11 = 0;
        secElecE12 = 0;
        secElecE21 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK2][valueCumulCS21];
        secElecE22 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK2][valueCumulCS22];

        valueCumulCS11 = 0;
        valueCumulCS12 = 0;
    }
    else
    {
        // No special case, interpolation will happen as usual.
        secElecE11 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK1][valueCumulCS11];
        secElecE12 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK1][valueCumulCS12];
        secElecE21 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK2][valueCumulCS21];
        secElecE22 = fEnergySecondaryData[materialName][particleName][ionizationLevelIndex][valueK2][valueCumulCS22];
    }

    ejectedElectronEnergy = QuadInterpolator(valueCumulCS11, valueCumulCS12,
                                      valueCumulCS21, valueCumulCS22,
                                      secElecE11, secElecE12,
                                      secElecE21, secElecE22,
                                      valueK1, valueK2,
                                      k, random);

    // **********************************************
    // Some tests for debugging
    // **********************************************

    G4double bindingEnergy (ptbStructure.IonisationEnergy(ionizationLevelIndex, materialName)/eV);
    if(k-ejectedElectronEnergy-bindingEnergy<=0 || ejectedElectronEnergy<=0)
    {
        G4cout<<"k "<<k<<G4endl;
        G4cout<<"material "<<materialName<<G4endl;
        G4cout<<"secondaryKin "<<ejectedElectronEnergy<<G4endl;
        G4cout<<"shell "<<ionizationLevelIndex<<G4endl;
        G4cout<<"bindingEnergy "<<bindingEnergy<<G4endl;
        G4cout<<"scatteredEnergy "<<k-ejectedElectronEnergy-bindingEnergy<<G4endl;
        G4cout<<"rand "<<random<<G4endl;
        G4cout<<"surrounding k values: valueK1 valueK2\n"<<valueK1<<" "<<valueK2<<G4endl;
        G4cout<<"surrounding E values: secElecE11 secElecE12 secElecE21 secElecE22\n"
             <<secElecE11<<" "<<secElecE12<<" "<<secElecE21<<" "<<secElecE22<<" "<<G4endl;
        G4cout<<"surrounding cumulCS values: valueCumulCS11 valueCumulCS12 valueCumulCS21 valueCumulCS22\n"
             <<valueCumulCS11<<" "<<valueCumulCS12<<" "<<valueCumulCS21<<" "<<valueCumulCS22<<" "<<G4endl;
        G4cerr<<"*****************************"<<G4endl;
        G4cerr<<"Fatal error, EXIT."<<G4endl;
        exit(EXIT_FAILURE);
    }

    return ejectedElectronEnergy*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBIonisationModel::LogLogInterpolate(G4double e1,
                                       G4double e2,
                                       G4double e,
                                       G4double xs1,
                                       G4double xs2)
{
    G4double value (0);

    // Switch to log-lin interpolation for faster code

    if ((e2-e1)!=0 && xs1 !=0 && xs2 !=0)
    {
        G4double d1 = std::log10(xs1);
        G4double d2 = std::log10(xs2);
        value = std::pow(10.,(d1 + (d2 - d1)*(e - e1)/ (e2 - e1)) );
    }

    // Switch to lin-lin interpolation for faster code
    // in case one of xs1 or xs2 (=cum proba) value is zero

    if ((e2-e1)!=0 && (xs1 ==0 || xs2 ==0))
    {
        G4double d1 = xs1;
        G4double d2 = xs2;
        value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    }

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBIonisationModel::QuadInterpolator(G4double e11, G4double e12,
                                      G4double e21, G4double e22,
                                      G4double xs11, G4double xs12,
                                      G4double xs21, G4double xs22,
                                      G4double t1, G4double t2,
                                      G4double t, G4double e)
{
    G4double interpolatedvalue1 (-1);
    if(xs11!=xs12) interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
    else interpolatedvalue1 = xs11;

    G4double interpolatedvalue2 (-1);
    if(xs21!=xs22) interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
    else interpolatedvalue2 = xs21;

    G4double value (-1);
    if(interpolatedvalue1!=interpolatedvalue2) value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
    else value = interpolatedvalue1;

    return value;

    //    G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
    //    G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
    //    G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
    //    return value;
}
