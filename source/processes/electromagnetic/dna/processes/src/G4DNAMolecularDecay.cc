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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
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
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4ITManager.hh"
#include "G4ParticleChange.hh"

using namespace std;

G4DNAMolecularDecay::G4DNAMolecularDecay(const G4String& processName,
                                                 G4ProcessType) : G4VITRestProcess(processName, fDecay)
{
    // set Process Sub Type
    SetProcessSubType(static_cast<int>(201));
    enableAlongStepDoIt = false;
    enablePostStepDoIt = false;
    enableAtRestDoIt=true;

    fVerbose = 1 ;

#ifdef G4VERBOSE
    if (verboseLevel>1)
    {
        G4cout << "G4MolecularDecayProcess constructor " << "  Name:" << processName << G4endl;
    }
#endif

    pParticleChange = &aParticleChange;

    fDecayAtFixedTime = true ;
    SetProcessSubType(60);
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
        if(verboseLevel>1)
        {
            G4cout<<"G4MolecularDecayProcess::IsApplicable(";
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

G4double G4DNAMolecularDecay::GetMeanLifeTime(const G4Track& /*track*/  ,
                                                  G4ForceCondition*)
{
    /* G4double output = GetMolecule(aTrack)-> GetDecayTime() - track.GetProperTime() ;
        return (output > 0 ? output : 0 );*/
    return 0.; // for water molecules in DNA
}

G4VParticleChange* G4DNAMolecularDecay::DecayIt(
    const G4Track& track,
    const G4Step& /*step*/)
{
    // DEBUG
    //    G4cout << "Is calling G4MolecularDecayProcess::DecayIt" << G4endl;

    aParticleChange.Initialize(track);
    const G4Molecule * theMotherMolecule = GetMolecule(track);
    const G4MoleculeDefinition* moleculeDefinition = theMotherMolecule->GetDefinition();

//    DEBUG
//    if(fVerbose > 1)
//    {
//        G4cout <<"Calling G4MolecularDecayProcess::DecayIt"<<G4endl;
//        G4cout << "The mother molecule state : " << G4endl;
//        theMotherMolecule -> PrintState();
//    }

    if(moleculeDefinition-> GetDecayTable())
    {
        const vector<const G4MolecularDecayChannel*>* DecayVector =
                (theMotherMolecule -> GetDecayChannel());

        G4int DecayVectorSize = DecayVector-> size();
//        DEBUG
//        if(fVerbose > 1)
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
//        if(fVerbose > 1)
//        {
//            G4cout<< "Selected Decay channel : " << decayChannel->GetName()<<G4endl;
//        }

        G4double decayEnergy = decayChannel->GetEnergy();
        G4int nbProducts = decayChannel->GetNbProducts();

        if(decayEnergy)
        {
//            DEBUG
//            if(fVerbose > 1)
//            {
//                G4cout<<"Deposit energy :" <<decayChannel->GetEnergy()/eV << " eV"<<G4endl;
//            }

            aParticleChange.ProposeLocalEnergyDeposit(decayChannel->GetEnergy());
        }

        if(nbProducts)
        {

//            DEBUG
//            if(fVerbose > 1)
//            {
//                G4cout<<"Number of products :" <<productsVectorSize<<G4endl;
//            }

            vector<G4ThreeVector> ProductsDisplacement(nbProducts);
            G4ThreeVector theMotherMoleculeDisplacement = G4ThreeVector(0,0,0);

            DecayDisplacementMap::iterator it = fDecayDisplacementMap.find(moleculeDefinition);

            if(it!=fDecayDisplacementMap.end())
            {
                G4VMolecularDecayDisplacer* displacer = it->second;
                ProductsDisplacement = displacer->GetProductsDisplacement(decayChannel);
                theMotherMoleculeDisplacement = displacer-> GetMotherMoleculeDisplacement(decayChannel);
            }
            else
            {
                G4String errMsg = "No G4MolecularDecayProcess::theDecayDisplacementMap[" ;
                errMsg += theMotherMolecule->GetName() +"]" ;
                G4Exception("G4MolecularDecayProcess::DecayIt","",FatalErrorInArgument, errMsg);
            }
//          DEBUG
//          G4cout<<"ProductsDisplacement size = " << ProductsDisplacement.size()<<G4endl;

            aParticleChange.SetNumberOfSecondaries(nbProducts);

            if(fVerbose)
            {
                G4cout<<"Decay Process : "
                     << theMotherMolecule->GetName()
                     <<" "<< decayChannel->GetName() ;
//                DEBUG
//                if(fVerbose>1)
//                {
//                    G4cout<<" | Track ID : "<< track.GetTrackID()
//                         <<" | Position : "<<track.GetPosition() / nanometer << "[nm]";
//                }
                G4cout<< G4endl;
            }

            for (G4int j=0; j<nbProducts ; j++)
            {
                G4Molecule* product = new G4Molecule(*decayChannel->GetProduct(j));

                // create a new track object
                G4Track* secondary = product->BuildTrack(picosecond/*track.GetGlobalTime()*/,track.GetPosition()
                                                         + theMotherMoleculeDisplacement + ProductsDisplacement[j]);
//                DEBUG
//                if(fVerbose>1)
//                {
//                    G4cout<<"theMotherMoleculeDisplacement : "<<theMotherMoleculeDisplacement<<G4endl;
//                    G4cout<<"ProductsDisplacement["<<j<<"] : "<<ProductsDisplacement[j]<<G4endl;
//                }

                secondary-> SetTrackStatus(fAlive);

                if(fVerbose)
                {
                    G4cout<<"Product : "<< product->GetName()<<G4endl;
                    // DEBUG
                    // <<" | Position : "<< secondary->GetPosition()/ nanometer << "[nm]"
                    // <<G4endl;
                }

                // add the secondary track in the List
                aParticleChange.G4VParticleChange::AddSecondary(secondary);
                G4ITManager<G4Molecule>::Instance()->Push(secondary);
            }

            if(fVerbose)
                G4cout<<"-------------"<<G4endl;
        }
//        DEBUG
//        else if(fVerbose > 1 && decayEnergy)
//        {
//            G4cout << "No products for this channel" << G4endl ;
//            G4cout<<"-------------"<<G4endl;
//        }
        else if(!decayEnergy && !nbProducts)
        {
            G4String errMsg = "There is no products and no energy specified in the molecular decay channel";
            G4Exception("G4MolecularDecayProcess::DecayIt","",FatalErrorInArgument, errMsg);
        }
    }

    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

void G4DNAMolecularDecay::SetDecayDisplacer (const G4MoleculeDefinition* molDef, G4VMolecularDecayDisplacer* aDisplacer)
{
    fDecayDisplacementMap[molDef] = aDisplacer;
}

G4VMolecularDecayDisplacer* G4DNAMolecularDecay::GetDecayDisplacer(const G4MoleculeDefinition* molDef)
{
    return fDecayDisplacementMap[molDef] ;
}
