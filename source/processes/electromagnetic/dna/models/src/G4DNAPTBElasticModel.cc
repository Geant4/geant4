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

#include "G4DNAPTBElasticModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Proton.hh"

G4DNAPTBElasticModel::G4DNAPTBElasticModel(const G4String& applyToMaterial, const G4ParticleDefinition*,
                                           const G4String& nam)
    : G4VDNAModel(nam, applyToMaterial)
{
    fKillBelowEnergy = 10*eV; // will be override by the limits defined for each material

    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 )
    {
        G4cout << "PTB Elastic model is constructed " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAPTBElasticModel::~G4DNAPTBElasticModel()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::Initialise(const G4ParticleDefinition* particle,
                                      const G4DataVector& /*cuts*/, G4ParticleChangeForGamma*)
{
    if (verboseLevel > 3)
        G4cout << "Calling G4DNAPTBElasticModel::Initialise()" << G4endl;

    G4double scaleFactor = 1e-16*cm*cm;

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();

    //*******************************************************
    // Cross section data
    //*******************************************************

    if(particle == electronDef)
    {
        G4String particleName = particle->GetParticleName();

        AddCrossSectionData("THF",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_THF",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_THF",
                            scaleFactor);
        SetLowELimit("THF", particleName, 10*eV);
        SetHighELimit("THF", particleName, 1*keV);

        AddCrossSectionData("PY",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PY",
                            scaleFactor);
        SetLowELimit("PY", particleName, 10*eV);
        SetHighELimit("PY", particleName, 1*keV);

        AddCrossSectionData("PU",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PU",
                            scaleFactor);
        SetLowELimit("PU", particleName, 10*eV);
        SetHighELimit("PU", particleName, 1*keV);

        AddCrossSectionData("TMP",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_TMP",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_TMP",
                            scaleFactor);
        SetLowELimit("TMP", particleName, 10*eV);
        SetHighELimit("TMP", particleName, 1*keV);

        AddCrossSectionData("G4_WATER",
                            particleName,
                            "dna/sigma_elastic_e_champion",
                            "dna/sigmadiff_cumulated_elastic_e_champion",
                            scaleFactor);
        SetLowELimit("G4_WATER", particleName, 10*eV);
        SetHighELimit("G4_WATER", particleName, 1*keV);

        // DNA materials
        //
        AddCrossSectionData("backbone_THF",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_THF",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_THF",
                            scaleFactor*33./30);
        SetLowELimit("backbone_THF", particleName, 10*eV);
        SetHighELimit("backbone_THF", particleName, 1*keV);

        AddCrossSectionData("cytosine_PY",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PY",
                            scaleFactor*42./30);
        SetLowELimit("cytosine_PY", particleName, 10*eV);
        SetHighELimit("cytosine_PY", particleName, 1*keV);

        AddCrossSectionData("thymine_PY",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PY",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PY",
                            scaleFactor*48./30);
        SetLowELimit("thymine_PY", particleName, 10*eV);
        SetHighELimit("thymine_PY", particleName, 1*keV);

        AddCrossSectionData("adenine_PU",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PU",
                            scaleFactor*50./44);
        SetLowELimit("adenine_PU", particleName, 10*eV);
        SetHighELimit("adenine_PU", particleName, 1*keV);

        AddCrossSectionData("guanine_PU",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_PU",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_PU",
                            scaleFactor*56./44);
        SetLowELimit("guanine_PU", particleName, 10*eV);
        SetHighELimit("guanine_PU", particleName, 1*keV);

        AddCrossSectionData("backbone_TMP",
                            particleName,
                            "dna/sigma_elastic_e-_PTB_TMP",
                            "dna/sigmadiff_cumulated_elastic_e-_PTB_TMP",
                            scaleFactor*33./50);
        SetLowELimit("backbone_TMP", particleName, 10*eV);
        SetHighELimit("backbone_TMP", particleName, 1*keV);
    }

    //*******************************************************
    // Load the data
    //*******************************************************

    LoadCrossSectionData(particle->GetParticleName() );

    //*******************************************************
    // Verbose output
    //*******************************************************

    if (verboseLevel > 2)
        G4cout << "Loaded cross section files for PTB Elastic model" << G4endl;

    if( verboseLevel>0 )
    {
        G4cout << "PTB Elastic model is initialized " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::ReadDiffCSFile(const G4String& materialName,
                                             const G4String& particleName,
                                             const G4String& file,
                                             const G4double)
{
    // Method to read and save the information contained within the differential cross section files.
    // This method is not yet standard.

    // get the path of the G4LEDATA data folder
    char *path = getenv("G4LEDATA");
    // if it is not found then quit and print error message
    if(!path)
    {
        G4Exception("G4DNAPTBElasticModel::ReadAllDiffCSFiles","em0006",
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
        G4Exception("G4DNAPTBElasticModel::Initialise","em0003",
                    FatalException, endPath.str().c_str());
    }

    tValuesVec[materialName][particleName].push_back(0.);

    G4String line;

    // read the file line by line until we reach the end of file point
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

        // Variables to be filled by the input file
        double tDummy;
        double eDummy;

        // fill the variables with the content of the line
        iss>>tDummy>>eDummy;

        // SI : mandatory Vecm initialization

        // Fill two vectors contained in maps of types:
        // [materialName][particleName]=vector
        // [materialName][particleName][T]=vector
        // to list all the incident energies (tValues) and all the output energies (eValues) within the file
        //
        // Check if we already have the current T value in the vector.
        // If not then add it
        if (tDummy != tValuesVec[materialName][particleName].back())
        {
            // Add the current T value
            tValuesVec[materialName][particleName].push_back(tDummy);

            // Make it correspond to a default zero E value
            eValuesVect[materialName][particleName][tDummy].push_back(0.);
        }

        // Put the differential cross section value of the input file within the diffCrossSectionData map
        iss>>diffCrossSectionData[materialName][particleName][tDummy][eDummy];

        // If the current E value (eDummy) is different from the one already registered in the eVector then add it to the vector
        if (eDummy != eValuesVect[materialName][particleName][tDummy].back()) eValuesVect[materialName][particleName][tDummy].push_back(eDummy);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::CrossSectionPerVolume(const G4Material* /*material*/,
                                                     const G4String& materialName,
                                                     const G4ParticleDefinition* p,
                                                     G4double ekin,
                                                     G4double /*emin*/,
                                                     G4double /*emax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNAPTBElasticModel" << G4endl;

    // Get the name of the current particle
    const G4String& particleName = p->GetParticleName();

    // set killBelowEnergy value for current material
    fKillBelowEnergy = GetLowELimit(materialName, particleName);

    // initialise the return value (cross section) to zero
    G4double sigma(0);

    // check if we are below the high energy limit
    if (ekin < GetHighELimit(materialName, particleName) )
    {
        // This is used to kill the particle if its kinetic energy is below fKillBelowEnergy.
        // If the energy is lower then we return a maximum cross section and thus the SampleSecondaries method will be called for sure.
        // SampleSecondaries will remove the particle from the simulation.
        //
        //SI : XS must not be zero otherwise sampling of secondaries method ignored
        if (ekin < fKillBelowEnergy) return DBL_MAX;

        // Get the tables with the cross section data
        TableMapData* tableData = GetTableData();

        // Retrieve the cross section value
        sigma = (*tableData)[materialName][particleName]->FindValue(ekin);
    }

    if (verboseLevel > 2)
    {
        G4cout << "__________________________________" << G4endl;
        G4cout << "°°° G4DNAPTBElasticModel - XS INFO START" << G4endl;
        G4cout << "°°° Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
        G4cout << "°°° Cross section per molecule (cm^2)=" << sigma/cm/cm << G4endl;
        G4cout << "°°° G4DNAPTBElasticModel - XS INFO END" << G4endl;
    }

    // Return the cross section
    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                             const G4MaterialCutsCouple* /*couple*/,
                                             const G4String& materialName,
                                             const G4DynamicParticle* aDynamicElectron,
                                             G4ParticleChangeForGamma* particleChangeForGamma,
                                             G4double /*tmin*/,
                                             G4double /*tmax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNAPTBElasticModel" << G4endl;

    G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();

    const G4String& particleName = aDynamicElectron->GetParticleDefinition()->GetParticleName();

    // set killBelowEnergy value for material
    fKillBelowEnergy = GetLowELimit(materialName, particleName);

    // If the particle (electron here) energy is below the kill limit then we remove it from the simulation
    if (electronEnergy0 < fKillBelowEnergy)
    {
        particleChangeForGamma->SetProposedKineticEnergy(0.);
        particleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        particleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
    }
    // If we are above the kill limite and below the high limit then we proceed
    else if (electronEnergy0>= fKillBelowEnergy && electronEnergy0 < GetHighELimit(materialName, particleName) )
    {
        // Random sampling of the cosTheta
        G4double cosTheta = RandomizeCosTheta(electronEnergy0, materialName);

        // Random sampling of phi
        G4double phi = 2. * pi * G4UniformRand();

        G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
        G4ThreeVector xVers = zVers.orthogonal();
        G4ThreeVector yVers = zVers.cross(xVers);

        G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
        G4double yDir = xDir;
        xDir *= std::cos(phi);
        yDir *= std::sin(phi);

        // Particle direction after ModelInterface
        G4ThreeVector zPrikeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

        // Give the new direction
        particleChangeForGamma->ProposeMomentumDirection(zPrikeVers.unit()) ;

        // Update the energy which does not change here
        particleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::Theta
(G4ParticleDefinition * particleDefinition, G4double k, G4double integrDiff, const G4String& materialName)
{  
    G4double theta = 0.;
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
    G4String particleName = particleDefinition->GetParticleName();

    if (particleDefinition == G4Electron::ElectronDefinition())
    {
        std::vector<double>::iterator t2 = std::upper_bound(tValuesVec[materialName][particleName].begin(),tValuesVec[materialName][particleName].end(), k);
        std::vector<double>::iterator t1 = t2-1;

        std::vector<double>::iterator e12 = std::upper_bound(eValuesVect[materialName][particleName][(*t1)].begin(),eValuesVect[materialName][particleName][(*t1)].end(), integrDiff);
        std::vector<double>::iterator e11 = e12-1;

        std::vector<double>::iterator e22 = std::upper_bound(eValuesVect[materialName][particleName][(*t2)].begin(),eValuesVect[materialName][particleName][(*t2)].end(), integrDiff);
        std::vector<double>::iterator e21 = e22-1;

        valueT1  =*t1;
        valueT2  =*t2;
        valueE21 =*e21;
        valueE22 =*e22;
        valueE12 =*e12;
        valueE11 =*e11;

        xs11 = diffCrossSectionData[materialName][particleName][valueT1][valueE11];
        xs12 = diffCrossSectionData[materialName][particleName][valueT1][valueE12];
        xs21 = diffCrossSectionData[materialName][particleName][valueT2][valueE21];
        xs22 = diffCrossSectionData[materialName][particleName][valueT2][valueE22];
    }

    if (xs11==0 && xs12==0 && xs21==0 && xs22==0) return (0.);

    theta = QuadInterpolator   ( valueE11, valueE12,
                                 valueE21, valueE22,
                                 xs11, xs12,
                                 xs21, xs22,
                                 valueT1, valueT2,
                                 k, integrDiff );

    return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LinLogInterpolate(G4double e1,
                                                 G4double e2,
                                                 G4double e,
                                                 G4double xs1,
                                                 G4double xs2)
{
    G4double d1 = std::log(xs1);
    G4double d2 = std::log(xs2);
    G4double value = std::exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LinLinInterpolate(G4double e1,
                                                 G4double e2,
                                                 G4double e,
                                                 G4double xs1,
                                                 G4double xs2)
{
    G4double d1 = xs1;
    G4double d2 = xs2;
    G4double value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LogLogInterpolate(G4double e1,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::QuadInterpolator(G4double e11, G4double e12,
                                                G4double e21, G4double e22,
                                                G4double xs11, G4double xs12,
                                                G4double xs21, G4double xs22,
                                                G4double t1, G4double t2,
                                                G4double t, G4double e)
{
    // Log-Log
 /*
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);


  // Lin-Log
  G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
*/

    // Lin-Lin
    G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
    G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
    G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::RandomizeCosTheta(G4double k, const G4String& materialName)
{
    G4double integrdiff=0;
    G4double uniformRand=G4UniformRand();
    integrdiff = uniformRand;

    G4double theta=0.;
    G4double cosTheta=0.;
    theta = Theta(G4Electron::ElectronDefinition(),k/eV,integrdiff, materialName);

    cosTheta= std::cos(theta*pi/180);

    return cosTheta;
}




