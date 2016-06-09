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
// $Id: G4DNAMolecularDecay.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAMolecularDecay.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4ITManager.hh"
#include "G4ParticleChange.hh"

using namespace std;

G4DNAMolecularDecay::G4DNAMolecularDecay(const G4String& processName,
                                         G4ProcessType type) : G4VITRestProcess(processName, type)
{
    // set Process Sub Type
    SetProcessSubType(59); // DNA sub-type
    enableAlongStepDoIt = false;
    enablePostStepDoIt = false;
    enableAtRestDoIt=true;

    fVerbose = 0 ;

#ifdef G4VERBOSE
    if (verboseLevel>1)
    {
        G4cout << "G4MolecularDecayProcess constructor " << "  Name:" << processName << G4endl;
    }
#endif

    pParticleChange = &aParticleChange;

    fDecayAtFixedTime = true ;
}

G4DNAMolecularDecay::~G4DNAMolecularDecay()
{
    DecayDisplacementMap::iterator it = fDecayDisplacementMap.begin();

    for( ; it != fDecayDisplacementMap.end() ; it++)
    {
        if(it->second)
        {
            delete it->second;
            it->second = 0;
        }
    }
    fDecayDisplacementMap.clear();
}

G4DNAMolecularDecay::G4DNAMolecularDecay(const G4DNAMolecularDecay &right) :
    G4VITRestProcess(right)
{
    fDecayAtFixedTime = right . fDecayAtFixedTime;
    fDecayDisplacementMap = right.fDecayDisplacementMap;
    fVerbose = right.fVerbose ;
}

G4bool G4DNAMolecularDecay::IsApplicable(const G4ParticleDefinition& aParticleType)
{
    if(aParticleType.GetParticleType()=="Molecule")
    {
#ifdef G4VERBOSE
        if(fVerbose>1)
        {
            G4cout<<"G4MolecularDecay::IsApplicable(";
            G4cout<<aParticleType.GetParticleName()<<",";
            G4cout<<aParticleType.GetParticleType()<<")"<<G4endl;
        }
#endif
        return(true);
    }
    else
    {
        return false;
    }
}

G4double G4DNAMolecularDecay::GetMeanLifeTime(const G4Track& track  ,
                                              G4ForceCondition*)
{
    G4double output = GetMolecule(track)-> GetDecayTime() - track.GetProperTime() ;
    return (output > 0 ? output : 0 );
}

G4VParticleChange* G4DNAMolecularDecay::DecayIt(
    const G4Track& track,
    const G4Step& )
{
    // DEBUG
    //    G4cout << "Is calling G4MolecularDecayProcess::DecayIt" << G4endl;

    aParticleChange.Initialize(track);
    const G4Molecule * theMotherMolecule = GetMolecule(track);
    const G4MoleculeDefinition* moleculeDefinition = theMotherMolecule->GetDefinition();

    //    DEBUG
    //        G4cout <<"Calling G4MolecularDecayProcess::DecayIt"<<G4endl;
    //        G4cout << "The mother molecule state : " << G4endl;
    //        theMotherMolecule -> PrintState();

    if(moleculeDefinition-> GetDecayTable())
    {
        const vector<const G4MolecularDecayChannel*>* DecayVector =
                (theMotherMolecule -> GetDecayChannel());

        if(DecayVector == 0)
        {
            G4ExceptionDescription exceptionDescription;
            theMotherMolecule->GetElectronOccupancy()->DumpInfo();
            exceptionDescription << "No decay channel was found for the molecule : " << theMotherMolecule-> GetName() << G4endl;
            G4Exception("G4DNAMolecularDecay::DecayIt", "G4DNAMolecularDecay::NoDecayChannel",FatalException,exceptionDescription);
            return &aParticleChange;
        }

        G4int DecayVectorSize = DecayVector-> size();
        //        DEBUG
        //            G4cout<< "Number of decay channels : " << DecayVectorSize<<G4endl;
        G4double RdmValue = G4UniformRand();

        const G4MolecularDecayChannel* decayChannel(0);
        G4int i=0;
        do
        {
            decayChannel = (*DecayVector)[i];
            if(RdmValue < decayChannel->GetProbability()) break;
            RdmValue -= decayChannel->GetProbability();
            i++;
        }
        while(i< DecayVectorSize);

        //        DEBUG
        //            G4cout<< "Selected Decay channel : " << decayChannel->GetName()<<G4endl;

        G4double decayEnergy = decayChannel->GetEnergy();
        G4int nbProducts = decayChannel->GetNbProducts();

        if(decayEnergy)
        {
            //            DEBUG
            //                G4cout<<"Deposit energy :" <<decayChannel->GetEnergy()/eV << " eV"<<G4endl;

            aParticleChange.ProposeLocalEnergyDeposit(decayChannel->GetEnergy());
        }

        if(nbProducts)
        {

            //            DEBUG
            //                G4cout<<"Number of products :" <<nbProducts<<G4endl;

            vector<G4ThreeVector> ProductsDisplacement(nbProducts);
            G4ThreeVector theMotherMoleculeDisplacement;

            DecayDisplacementMap::iterator it = fDecayDisplacementMap.find(moleculeDefinition);

            if(it!=fDecayDisplacementMap.end())
            {
                G4VMolecularDecayDisplacer* displacer = it->second;
                ProductsDisplacement = displacer->GetProductsDisplacement(decayChannel);
                theMotherMoleculeDisplacement = displacer-> GetMotherMoleculeDisplacement(decayChannel);
            }
            else
            {
                G4ExceptionDescription errMsg;
                errMsg << "No G4MolecularDecayProcess::theDecayDisplacementMap["
                       << theMotherMolecule->GetName() +"]" ;
                G4Exception("G4MolecularDecayProcess::DecayIt","DNAMolecularDecay001",FatalErrorInArgument, errMsg);
            }

            aParticleChange.SetNumberOfSecondaries(nbProducts);

#ifdef G4VERBOSE
            if(fVerbose)
            {
                G4cout<<"Decay Process : "
                     << theMotherMolecule->GetName()
                     << " (trackID :" << track.GetTrackID() << ") "
                     << decayChannel->GetName()
                     << G4endl;
            }
#endif

            for (G4int j=0; j<nbProducts ; j++)
            {
                G4Molecule* product = new G4Molecule(*decayChannel->GetProduct(j));

                // create a new track object
                // Be carefull as this processes is dedicated to be used in the DNA module
                // The molecular decay will happen one picosecond after the start of the simulation.
                // This may be seen as a bug and will be hopefully improve in the next releases
                G4Track* secondary = product->BuildTrack(picosecond,track.GetPosition()
                                                         + theMotherMoleculeDisplacement + ProductsDisplacement[j]);

                secondary-> SetTrackStatus(fAlive);
#ifdef G4VERBOSE
                if(fVerbose)
                {
                    G4cout<<"Product : "<< product->GetName()<<G4endl;
                }
#endif
                // add the secondary track in the List
                aParticleChange.G4VParticleChange::AddSecondary(secondary);
                G4ITManager<G4Molecule>::Instance()->Push(secondary);
            }
#ifdef G4VERBOSE
            if(fVerbose)
                G4cout<<"-------------"<<G4endl;
#endif
        }
        //        DEBUG
        //        else if(decayEnergy)
        //        {
        //            G4cout << "No products for this channel" << G4endl ;
        //            G4cout<<"-------------"<<G4endl;
        //        }
        else if(!decayEnergy && !nbProducts)
        {
            G4ExceptionDescription errMsg;
            errMsg << "There is no products and no energy specified in the molecular decay channel";
            G4Exception("G4MolecularDecayProcess::DecayIt","DNAMolecularDecay002",FatalErrorInArgument, errMsg);
        }
    }

    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

void G4DNAMolecularDecay::SetDecayDisplacer (const G4ParticleDefinition* molDef, G4VMolecularDecayDisplacer* aDisplacer)
{
    fDecayDisplacementMap[molDef] = aDisplacer;
}

G4VMolecularDecayDisplacer* G4DNAMolecularDecay::GetDecayDisplacer(const G4ParticleDefinition* molDef)
{
    return fDecayDisplacementMap[molDef] ;
}
