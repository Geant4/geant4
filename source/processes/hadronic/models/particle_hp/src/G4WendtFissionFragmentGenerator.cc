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
/*
 * File:   G4WendtFissionFragmentGenerator.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 21, 2013, 13:58 MST
 */

#include "G4ParticleHPManager.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4WendtFissionFragmentGenerator.hh"

G4ThreadLocal G4WendtFissionFragmentGenerator* G4WendtFissionFragmentGenerator::instance = NULL;

G4WendtFissionFragmentGenerator::
G4WendtFissionFragmentGenerator()
{
    // Set the default verbosity
    Verbosity_ = G4FFGDefaultValues::Verbosity;
}
/*
G4WendtFissionFragmentGenerator* G4WendtFissionFragmentGenerator::
GetInstance()
{
    //static G4WendtFissionFragmentGenerator newMe;
    //
    //return &newMe;

    if ( instance == NULL) instance = new G4WendtFissionFragmentGenerator();

    return instance;
}
*/
G4HadFinalState* G4WendtFissionFragmentGenerator::
ApplyYourself(const G4HadProjectile& projectile, G4int Z, G4int A)
{
G4FFG_FUNCTIONENTER__

    G4HadFinalState* finalState = NULL;
    G4DynamicParticleVector* finalParticles = NULL;
    G4int isotope;
    std::map< const G4int, G4FissionFragmentGenerator* >::iterator fissionGenerator;

    // Look for the first available isomer since no M is provided for ApplyYourself()
    for(unsigned int M = 0; M < 10; ++M)
    {
        isotope = G4FissionFragmentGenerator::G4MakeIsotopeCode(Z, A, M);
        fissionGenerator = fissionIsotopes.find(isotope);

        if(fissionGenerator != fissionIsotopes.end())
        {
            // Only generate particles if the generator was constructed
            if(fissionGenerator->second)
            {
                finalParticles = fissionGenerator->second->G4GenerateFission(projectile);
            }

            break;
        }
    }

    if(finalParticles)
    {
        finalState = new G4HadFinalState();

        for(unsigned int i = 0; i < finalParticles->size(); ++i)
        {
            finalState->AddSecondary((*finalParticles)[i]);
        }
    }

    //TK modified 131108 add next line 
    //TK 160112 fix for coverity #53481
    if ( finalState != NULL ) finalState->SetStatusChange(stopAndKill);
G4FFG_FUNCTIONLEAVE__
    return finalState;
}

void G4WendtFissionFragmentGenerator::
InitializeANucleus(const G4int A, const G4int Z, const G4int M, const G4String& dataDirectory)
{
//G4FFG_FUNCTIONENTER__

    const G4int isotope = G4FissionFragmentGenerator::G4MakeIsotopeCode(Z, A, M);
    G4FFGEnumerations::MetaState metaState;
    std::pair< std::map< const G4int, G4FissionFragmentGenerator* >::iterator, bool > newIsotope;

    // Check to see if the isotope/isomer alread exists in the table
    newIsotope = fissionIsotopes.insert(std::make_pair(isotope, (G4FissionFragmentGenerator*)NULL));

    if(newIsotope.second || newIsotope.first->second == NULL)
    {
        // Get the data file
        G4bool flag;
        G4ParticleHPDataUsed dataFile = fileNames.GetName(A, Z, M, dataDirectory, "FF", flag);
        G4String dataFileName = dataFile.GetName();

        // Check if the file exists, and do not create a fission object if it doesn't
        // G4cout << "*** Z = " << Z << "\tA = " << A << "\t\t\t Directory: "<< dataDirectory << " DATA FILE: " << dataFileName << G4endl;
        std::istringstream dataStream(std::ios::in);
        G4ParticleHPManager::GetInstance()->GetDataStream(dataFileName, dataStream);
        if(!dataStream)
        {
            //G4FFG_FUNCTIONLEAVE__
            // G4cerr << "*** Stream error" << G4endl;
            return;
        }

        // Check the data file parameters
        if(!flag
           || ( Z < 2.5 && ( (G4double)abs( dataFile.GetZ() - Z ) > 0.001 || (G4double)abs( (G4int)dataFile.GetA() - A ) > 0.0001 ) ) )
        {
            //G4cerr << "*** Something wrong with the data request.\tFlag :" << flag << G4endl;
            //G4FFG_FUNCTIONLEAVE__
            return;
        }

        G4FissionFragmentGenerator* const fissionGenerator = new G4FissionFragmentGenerator();
        newIsotope.first->second = fissionGenerator;

        switch(M)
        {
        case 1:
            metaState = G4FFGEnumerations::META_1;
            break;

        case 2:
            metaState = G4FFGEnumerations::META_2;
            break;

        default:
            // TODO Display a warning message here indicating that an invalid metastate was passed in
            // Fall through to the ground state by default
        case 0:
            metaState = G4FFGEnumerations::GROUND_STATE;
            break;
        }

        fissionGenerator->G4SetIsotope(isotope);
        fissionGenerator->G4SetMetaState(metaState);
        fissionGenerator->G4SetCause(G4FFGEnumerations::NEUTRON_INDUCED);
        // TODO Load all the fission data and use the projectile energy instead
        fissionGenerator->G4SetIncidentEnergy(G4FFGDefaultValues::ThermalNeutronEnergy);
        fissionGenerator->G4SetYieldType(G4FFGEnumerations::INDEPENDENT);
        fissionGenerator->G4SetSamplingScheme(G4FFGEnumerations::NORMAL);


        // TODO Remove the need for forcing a load in the initialization phase,
        //      i.e. remove the ability to dynamically change the fission parameters
        //      that cause reload because a G4FissionFragmentGenerator class for
        //      each isotope should be loaded in the initialization phase
        if(!fissionGenerator->InitializeFissionProductYieldClass(dataStream))
        {
            // Delete if the initialization fails
            delete fissionGenerator;

            fissionIsotopes.erase(newIsotope.first);
        }
    }

//G4FFG_FUNCTIONLEAVE__
}

G4WendtFissionFragmentGenerator::
~G4WendtFissionFragmentGenerator()
{
    std::map< const G4int, G4FissionFragmentGenerator* >::iterator fissionGenerator;

    for(fissionGenerator = fissionIsotopes.begin(); fissionGenerator != fissionIsotopes.end(); ++fissionGenerator)
    {
        delete fissionGenerator->second;
    }
}
