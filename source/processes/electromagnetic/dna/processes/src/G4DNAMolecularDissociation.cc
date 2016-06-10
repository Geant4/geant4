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
// $Id: G4DNAMolecularDissociation.cc 93936 2015-11-04 09:37:59Z gcosmo $
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

#include "G4DNAMolecularDissociation.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4ParticleChange.hh"
#include "G4ITTransportationManager.hh"
#include "G4ITNavigator.hh"

using namespace std;

//______________________________________________________________________________

G4DNAMolecularDissociation::
G4DNAMolecularDissociation(const G4String& processName,
                           G4ProcessType type) :
    G4VITRestDiscreteProcess(processName, type)
{
  // set Process Sub Type
  SetProcessSubType(59); // DNA sub-type
  enableAlongStepDoIt = false;
  enablePostStepDoIt = true;
  enableAtRestDoIt = true;

  fVerbose = 0;

#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4MolecularDissociationProcess constructor " << "  Name:"
           << processName << G4endl;
  }
#endif

  pParticleChange = &aParticleChange;

  fDecayAtFixedTime = true;
  fProposesTimeStep = true;
}

//______________________________________________________________________________

G4DNAMolecularDissociation::~G4DNAMolecularDissociation()
{
  DecayDisplacementMap::iterator it = fDecayDisplacementMap.begin();

  for(; it != fDecayDisplacementMap.end(); it++)
  {
    if(it->second)
    {
      delete it->second;
      it->second = 0;
    }
  }
  fDecayDisplacementMap.clear();
}

//______________________________________________________________________________

G4DNAMolecularDissociation::
G4DNAMolecularDissociation(const G4DNAMolecularDissociation &right) :
    G4VITRestDiscreteProcess(right)
{
  fDecayAtFixedTime = right.fDecayAtFixedTime;
  fDecayDisplacementMap = right.fDecayDisplacementMap;
  fVerbose = right.fVerbose;
}

//______________________________________________________________________________

G4bool G4DNAMolecularDissociation::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if(aParticleType.GetParticleType() == "Molecule")
  {
#ifdef G4VERBOSE

    if(fVerbose > 1)
    {
      G4cout << "G4MolecularDissociation::IsApplicable(";
      G4cout << aParticleType.GetParticleName() << ",";
      G4cout << aParticleType.GetParticleType() << ")" << G4endl;
    }
#endif
    return (true);
  }
  else
  {
    return false;
  }
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::GetMeanLifeTime(const G4Track& track,
                                                     G4ForceCondition*)
{
  G4double output = GetMolecule(track)->GetDecayTime() - track.GetProperTime();
  return (output > 0 ? output : 0);
}

//______________________________________________________________________________

G4VParticleChange* G4DNAMolecularDissociation::DecayIt(const G4Track& track,
                                                       const G4Step&)
{
  // DEBUG
  // G4cout << "Is calling G4MolecularDecayProcess::DecayIt" << G4endl;

  aParticleChange.Initialize(track);
  const G4Molecule * theMotherMolecule = GetMolecule(track);
  const G4MoleculeDefinition* moleculeDefinition = theMotherMolecule
      ->GetDefinition();

  // DEBUG
  // G4cout <<"Calling G4MolecularDecayProcess::DecayIt"<<G4endl;
  // G4cout << "The mother molecule state : " << G4endl;
  // theMotherMolecule -> PrintState();

  if(moleculeDefinition->GetDecayTable())
  {
    const vector<const G4MolecularDissociationChannel*>* DecayVector =
        (theMotherMolecule->GetDecayChannel());

    if(DecayVector == 0)
    {
      G4ExceptionDescription exceptionDescription;
      theMotherMolecule->PrintState();
      exceptionDescription << "No decay channel was found for the molecule : "
                           << theMotherMolecule->GetName() << G4endl;
      G4Exception("G4DNAMolecularDissociation::DecayIt",
                  "G4DNAMolecularDissociation::NoDecayChannel",
                  FatalException,
                  exceptionDescription);
      return &aParticleChange;
    }

    G4int DecayVectorSize = DecayVector->size();
    // DEBUG
    // G4cout<< "Number of decay channels : " << DecayVectorSize <<G4endl;
    G4double RdmValue = G4UniformRand();

    const G4MolecularDissociationChannel* decayChannel(0);
    G4int i = 0;
    do
    {
      decayChannel = (*DecayVector)[i];
      if(RdmValue < decayChannel->GetProbability()) break;
      RdmValue -= decayChannel->GetProbability();
      i++;
    }
    while(i < DecayVectorSize);

    // DEBUG
    // G4cout << "Selected Decay channel : "
    //        << decayChannel->GetName() << G4endl;

    G4double decayEnergy = decayChannel->GetEnergy();
    G4int nbProducts = decayChannel->GetNbProducts();

    if(decayEnergy)
    {
      // DEBUG
      // G4cout << "Deposit energy :"
      //        << decayChannel->GetEnergy()/eV << " eV" << G4endl;

      aParticleChange.ProposeLocalEnergyDeposit(decayChannel->GetEnergy());
    }

    if(nbProducts)
    {

      // DEBUG
      // G4cout << "Number of products :" << nbProducts << G4endl;

      vector<G4ThreeVector> ProductsDisplacement(nbProducts);
      G4ThreeVector theMotherMoleculeDisplacement;

      DecayDisplacementMap::iterator it =
          fDecayDisplacementMap.find(moleculeDefinition);

      if(it != fDecayDisplacementMap.end())
      {
        G4VMolecularDecayDisplacer* displacer = it->second;
        ProductsDisplacement = displacer->GetProductsDisplacement(decayChannel);
        theMotherMoleculeDisplacement =
            displacer->GetMotherMoleculeDisplacement(decayChannel);
      }
      else
      {
        G4ExceptionDescription errMsg;
        errMsg << "No G4MolecularDecayProcess::theDecayDisplacementMap["
               << theMotherMolecule->GetName() + "]";
        G4Exception("G4MolecularDecayProcess::DecayIt",
                    "DNAMolecularDecay001",
                    FatalErrorInArgument,
                    errMsg);
      }

      aParticleChange.SetNumberOfSecondaries(nbProducts);

//      G4cout << " nbProducts = " << nbProducts << G4endl;
//      theMotherMolecule->PrintState();

#ifdef G4VERBOSE
      if(fVerbose)
      {
        G4cout << "Decay Process : " << theMotherMolecule->GetName()
               << " (trackID :" << track.GetTrackID() << ") "
               << decayChannel->GetName() << G4endl;
      }
#endif

      G4ITNavigator* navigator =
          G4ITTransportationManager::GetTransportationManager()
              ->GetNavigatorForTracking();

      for(G4int j = 0; j < nbProducts; j++)
      {
        G4Molecule* product = new G4Molecule(decayChannel->GetProduct(j));

        G4ThreeVector displacement = theMotherMoleculeDisplacement
            + ProductsDisplacement[j];
        double mag_displacement = displacement.mag();
        G4ThreeVector displacement_direction = displacement / mag_displacement;

        double prNewSafety = DBL_MAX;

        navigator->CheckNextStep(track.GetPosition(),
                                 displacement_direction,
                                 mag_displacement,
                                 prNewSafety); // returns a value

        if(prNewSafety < mag_displacement) mag_displacement = prNewSafety;

//        const G4AffineTransform& transform = navigator
//            ->GetGlobalToLocalTransform();
//
//        G4ThreeVector localPoint =
//            transform.TransformPoint(track.GetPosition());
//
//        if(track.GetTouchable()->GetSolid()->Inside(localPoint) != EInside::kInside)
//        {
//          G4Exception("G4DNAMolecularDissociation::DecayIt",
//                      "OUTSIDE_OF_MOTHER_VOLUME",
//                      FatalException,
//                      "Product has been placed outside of the volume "
//                      "containing the mother molecule");
//        }

        G4Track* secondary =
            product->BuildTrack(track.GetGlobalTime(),
                                track.GetPosition() + displacement_direction
                                    * mag_displacement);

        secondary->SetTrackStatus(fAlive);
#ifdef G4VERBOSE
        if(fVerbose)
        {
          G4cout << "Product : " << product->GetName() << G4endl;
        }
#endif
        // add the secondary track in the List
        aParticleChange.G4VParticleChange::AddSecondary(secondary);
      }
#ifdef G4VERBOSE
      if(fVerbose) G4cout << "-------------" << G4endl;
#endif
    }
    //DEBUG
    else if(fVerbose && decayEnergy)
    {
      G4cout << "No products for this channel" << G4endl;
      G4cout<<"-------------"<<G4endl;
    }
    /*
     else if(!decayEnergy && !nbProducts)
     {
     G4ExceptionDescription errMsg;
     errMsg << "There is no products and no energy specified in the molecular "
     "dissociation channel";
     G4Exception("G4MolecularDissociationProcess::DecayIt",
     "DNAMolecularDissociation002",
     FatalErrorInArgument,
     errMsg);
     }
     */
  }

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  return &aParticleChange;
}

//______________________________________________________________________________

void G4DNAMolecularDissociation::
SetDecayDisplacer(const G4ParticleDefinition* molDef,
                  G4VMolecularDecayDisplacer* aDisplacer)
{
  fDecayDisplacementMap[molDef] = aDisplacer;
}

//______________________________________________________________________________

G4VMolecularDecayDisplacer*
G4DNAMolecularDissociation::
GetDecayDisplacer(const G4ParticleDefinition* molDef)
{
  return fDecayDisplacementMap[molDef];
}

//______________________________________________________________________________

G4double G4DNAMolecularDissociation::
PostStepGetPhysicalInteractionLength(const G4Track&,
                                     G4double,
                                     G4ForceCondition*)
{
  return 0; //-1*(-log(1-G4UniformRand())*100*1e-15*s);
}
