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
// CPA100 excitation model class for electrons
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

#include "G4DNACPA100ExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100ExcitationModel::G4DNACPA100ExcitationModel(const G4ParticleDefinition*,
                                                       const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
    fpMolWaterDensity = 0;

    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 )
    {
        G4cout << "CPA100 excitation model is constructed " << G4endl;
    }
    fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100ExcitationModel::~G4DNACPA100ExcitationModel()
{ 
    // Cross section

    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    for (pos = tableData.begin(); pos != tableData.end(); ++pos)
    {
        G4DNACrossSectionDataSet* table = pos->second;
        delete table;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                            const G4DataVector& /*cuts*/)
{

    if (verboseLevel > 3)
        G4cout << "Calling G4DNACPA100ExcitationModel::Initialise()" << G4endl;

    G4String fileElectron("dna/sigma_excitation_e_cpa100");

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();

    G4String electron;

    G4double scaleFactor = 1.e-20 *m*m;

    // *** ELECTRON

    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    lowEnergyLimit[electron] = 11 * eV;  // Default low enetgy limit
    //lowEnergyLimit[electron] = 10.481 * eV;  // From sigexi2 file by MCB
    highEnergyLimit[electron] = 255955 * eV; // idem

    // Cross section
    
    G4DNACrossSectionDataSet* tableE 
     = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor );
    
    /*
      G4DNACrossSectionDataSet* tableE = 
       new G4DNACrossSectionDataSet(new G4DNACPA100LogLogInterpolation, eV, scaleFactor );
    */
    
    tableE->LoadData(fileElectron);

    tableData[electron] = tableE;
    
    //

    if (particle==electronDef)
    {
        SetLowEnergyLimit(lowEnergyLimit[electron]);
        SetHighEnergyLimit(highEnergyLimit[electron]);
    }


//    if( verboseLevel>0 )
    {
        G4cout << "CPA100 excitation model is initialized " << G4endl
               << "Energy range: "
               << LowEnergyLimit() / eV << " eV - "
               << HighEnergyLimit() / keV << " keV for "
               << particle->GetParticleName()
               << G4endl;
    }

    // Initialize water density pointer
    fpMolWaterDensity = 
      G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

    if (isInitialised) { return; }
    fParticleChangeForGamma = GetParticleChangeForGamma();
    isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                           const G4ParticleDefinition* particleDefinition,
                                                           G4double ekin,
                                                           G4double,
                                                           G4double)
{

    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNACPA100ExcitationModel" << G4endl;

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
                G4Exception("G4DNACPA100ExcitationModel::CrossSectionPerVolume","em0002",
                            FatalException,"Model not applicable to particle type.");
            }
        }

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "G4DNACPA100ExcitationModel - XS INFO START" << G4endl;
            G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
            G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
            //      G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
            G4cout << "G4DNACPA100ExcitationModel - XS INFO END" << G4endl;
        }

    } // if (waterMaterial)

    return sigma*waterDensity;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* ,
                                                 const G4MaterialCutsCouple*,
                                                 const G4DynamicParticle* aDynamicParticle,
                                                 G4double,
                                                 G4double)
{

    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNACPA100ExcitationModel" << G4endl;

    G4double k = aDynamicParticle->GetKineticEnergy();

    const G4String& particleName = aDynamicParticle->GetDefinition()->GetParticleName();

    G4int level = RandomSelect(k,particleName);
    G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
    G4double newEnergy = k - excitationEnergy;

    if (newEnergy > 0)
    {
        // fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
        
        // We take into account direction change as described page 87 (II.92) in thesis by S. Edel

        G4double cosTheta = 

         (excitationEnergy/k) / (1. + (k/(2*electron_mass_c2))*(1.-excitationEnergy/k) );
  
        cosTheta = std::sqrt(1.-cosTheta);

        G4double phi = 2. * pi * G4UniformRand();

        G4ThreeVector zVers = aDynamicParticle->GetMomentumDirection();
 
        //G4ThreeVector xVers = zVers.orthogonal();
        //G4ThreeVector yVers = zVers.cross(xVers);
        //G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
        //G4double yDir = xDir;
        //xDir *= std::cos(phi);
        //yDir *= std::sin(phi);
        // G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

        // Computation of scattering angles (from Subroutine DIRAN in CPA100)

        G4double CT1, ST1, CF1, SF1, CT2, ST2, CF2, SF2;
        G4double sinTheta = std::sqrt (1-cosTheta*cosTheta);

        CT1=0;
        ST1=0;
        CF1=0;
        SF1=0;
        CT2=0;
        ST2=0;
        CF2=0;
        SF2=0;

        CT1 = zVers.z();
        ST1=std::sqrt(1.-CT1*CT1);

        if (ST1!=0) CF1 = zVers.x()/ST1; else CF1 = std::cos(2. * pi * G4UniformRand());
        if (ST1!=0) SF1 = zVers.y()/ST1; else SF1 = std::sqrt(1.-CF1*CF1);

        G4double A3, A4, A5, A2, A1;
        A3=0;
        A4=0;
        A5=0;
        A2=0;
        A1=0;

        A3 = sinTheta*std::cos(phi);
        A4 = A3*CT1 + ST1*cosTheta;
        A5 = sinTheta * std::sin(phi);
        A2 = A4 * SF1 + A5 * CF1;
        A1 = A4 * CF1 - A5 * SF1;

        CT2 = CT1*cosTheta - ST1*A3;
        ST2 = std::sqrt(1.-CT2*CT2);

        if (ST2==0) ST2=1E-6;
        CF2 = A1/ST2;
        SF2 = A2/ST2;

        /*
        G4cout << "CT1=" << CT1 << G4endl;
        G4cout << "ST1=" << ST1 << G4endl;
        G4cout << "CF1=" << CF1 << G4endl;
        G4cout << "SF1=" << SF1 << G4endl;
        G4cout << "cosTheta=" << cosTheta << G4endl;
        G4cout << "sinTheta=" << sinTheta << G4endl;
        G4cout << "cosPhi=" << std::cos(phi) << G4endl;
        G4cout << "sinPhi=" << std::sin(phi) << G4endl;
        G4cout << "CT2=" << CT2 << G4endl;
        G4cout << "ST2=" << ST2 << G4endl;
        G4cout << "CF2=" << CF2 << G4endl;
        G4cout << "SF2=" << SF2 << G4endl;
        */

        G4ThreeVector zPrimeVers(ST2*CF2,ST2*SF2,CT2);

        //

        fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit()) ;

        //

        if (!statCode) fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
        else fParticleChangeForGamma->SetProposedKineticEnergy(k);
 
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
    }

    // Chemistry
    
    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
                                                           level,
                                                           theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNACPA100ExcitationModel::RandomSelect(G4double k, const G4String& particle)
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
 G4double tmp=10.481*eV;
 G4cout <<  table->GetComponent(0)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(1)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(2)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(3)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  table->GetComponent(4)->FindValue(tmp)/(1e-20*m*m) << G4endl;
 G4cout <<  
 table->GetComponent(0)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(1)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(2)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(3)->FindValue(tmp)/(1e-20*m*m) +
 table->GetComponent(4)->FindValue(tmp)/(1e-20*m*m) 
 << G4endl;
 abort();
 */
 //
 //Dump
 //
 /*
 G4double minEnergy = 10.481  * eV;
        G4double maxEnergy = 255955. * eV;
        G4int nEnergySteps = 1000;
        G4double energy(minEnergy);
        G4double stpEnergy(std::pow(maxEnergy/energy, 1./static_cast<G4double>(nEnergySteps-1)));
        G4int step(nEnergySteps);
        system ("rm -rf excitation-cap100.out");
 FILE* myFile=fopen("excitation-cpa100.out","a");
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
        G4Exception("G4DNACPA100ExcitationModel::RandomSelect","em0002",
                    FatalException,"Model not applicable to particle type.");
    }
    return level;
}


