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
// $Id:$
//
// Author: Mathieu Karamitros
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disappear in the next releases.
//
// History:
// -----------
// 13 Nov 2016 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4TransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4Navigator.hh"

//#define MODEL_VERBOSE

//------------------------------------------------------------------------------

template<typename MODEL>
G4TDNAOneStepThermalizationModel<MODEL>::
G4TDNAOneStepThermalizationModel(const G4ParticleDefinition*,
                                const G4String& nam) :
G4VEmModel(nam), fIsInitialised(false)
{
  fVerboseLevel = 0;
  SetLowEnergyLimit(0.);
  G4DNAWaterExcitationStructure exStructure;
  SetHighEnergyLimit(exStructure.ExcitationEnergy(0));
  fParticleChangeForGamma = 0;
  fpWaterDensity = 0;
  fNavigator = 0;
}

//------------------------------------------------------------------------------

template<typename MODEL>
G4TDNAOneStepThermalizationModel<MODEL>::~G4TDNAOneStepThermalizationModel()
{
  if(fNavigator)
  {
    //    if(fNavigator->GetNavigatorState())
    //      delete fNavigator->GetNavigatorState();
    delete fNavigator;
  }
}

//------------------------------------------------------------------------------
template<typename MODEL>
void G4TDNAOneStepThermalizationModel<MODEL>::
Initialise(const G4ParticleDefinition* particleDefinition,
           const G4DataVector&)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel)
    G4cout << "Calling G4DNAOneStepThermalizationModel::Initialise()"
           << G4endl;
#endif
  if (particleDefinition->GetParticleName() != "e-")
  {
    G4ExceptionDescription errMsg;
    errMsg << "G4DNAOneStepThermalizationModel can only be applied "
    "to electrons";
    G4Exception("G4DNAOneStepThermalizationModel::CrossSectionPerVolume",
                "G4DNAOneStepThermalizationModel001",
                FatalErrorInArgument,errMsg);
    return;
  }
  
  if(!fIsInitialised)
  {
    fIsInitialised = true;
    fParticleChangeForGamma = GetParticleChangeForGamma();
  }
  
  G4Navigator* navigator =
  G4TransportationManager::GetTransportationManager()->
  GetNavigatorForTracking();
  
  fNavigator = new G4Navigator();
  
  if(navigator){ // add these checks for testing mode
    auto world=navigator->GetWorldVolume();
    if(world){
      fNavigator->SetWorldVolume(world);
      //fNavigator->NewNavigatorState();
    }
  }
  
  fpWaterDensity =
  G4DNAMolecularMaterial::Instance()->
  GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
}

//------------------------------------------------------------------------------
template<typename MODEL>
G4double G4TDNAOneStepThermalizationModel<MODEL>::
CrossSectionPerVolume(const G4Material* material,
                      const G4ParticleDefinition*,
                      G4double ekin,
                      G4double,
                      G4double)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel > 1)
    G4cout << "Calling CrossSectionPerVolume() of G4DNAOneStepThermalizationModel"
    << G4endl;
#endif
  
  if(ekin > HighEnergyLimit()){
    return 0.0;
  }
  
  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];
  
  if(waterDensity!= 0.0){
    return DBL_MAX;
  }
  return 0.;
}

//------------------------------------------------------------------------------
template<typename MODEL>
double G4TDNAOneStepThermalizationModel<MODEL>::GetRmean(double k){
  return MODEL::GetRmean(k);
}


//------------------------------------------------------------------------------

template<typename MODEL>
void G4TDNAOneStepThermalizationModel<MODEL>::
GetPenetration(G4double k, G4ThreeVector& displacement)
{
  return MODEL::GetPenetration(k, displacement);
}

//------------------------------------------------------------------------------
template<typename MODEL>
void G4TDNAOneStepThermalizationModel<MODEL>::
SampleSecondaries(std::vector<G4DynamicParticle*>*,
                  const G4MaterialCutsCouple*,
                  const G4DynamicParticle* particle,
                  G4double,
                  G4double)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel)
    G4cout << "Calling SampleSecondaries() of G4DNAOneStepThermalizationModel"
    << G4endl;
#endif
  
  G4double k = particle->GetKineticEnergy();
  
  if (k <= HighEnergyLimit())
  {
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    
    if(G4DNAChemistryManager::IsActivated())
    {
      G4ThreeVector displacement(0,0,0);
      GetPenetration(k, displacement);
      
      //______________________________________________________________
      const G4Track * theIncomingTrack =
      fParticleChangeForGamma->GetCurrentTrack();
      G4ThreeVector finalPosition(theIncomingTrack->GetPosition()+displacement);
      
      fNavigator->SetWorldVolume(theIncomingTrack->GetTouchable()->
                                 GetVolume(theIncomingTrack->GetTouchable()->
                                           GetHistoryDepth()));
      
      double displacementMag = displacement.mag();
      double safety = DBL_MAX;
      G4ThreeVector direction = displacement/displacementMag;
      
      //--
      // 6/09/16 - recupere de molecular dissocation
      double mag_displacement = displacement.mag();
      G4ThreeVector displacement_direction = displacement/mag_displacement;
      
      //     double step = DBL_MAX;
      //     step = fNavigator->CheckNextStep(theIncomingTrack->GetPosition(),
      //                                     displacement_direction,
      //                                     mag_displacement,
      //                                     safety);
      //
      //
      //     if(safety < mag_displacement)
      //     {
      ////       mag_displacement = prNewSafety;
      //       finalPosition = theIncomingTrack->GetPosition()
      //       + (displacement/displacementMag)*safety*0.80;
      //     }
      //--
      
      fNavigator->ResetHierarchyAndLocate(theIncomingTrack->GetPosition(),
                                          direction,
                                          *((G4TouchableHistory*)
                                            theIncomingTrack->GetTouchable()));
      
      fNavigator->ComputeStep(theIncomingTrack->GetPosition(),
                              displacement/displacementMag,
                              displacementMag,
                              safety);
      
      if(safety <= displacementMag)
      {
        finalPosition = theIncomingTrack->GetPosition()
        + (displacement/displacementMag)*safety*0.80;
      }
      
      G4DNAChemistryManager::Instance()->CreateSolvatedElectron(theIncomingTrack,
                                                                &finalPosition);
      
      fParticleChangeForGamma->SetProposedKineticEnergy(25.e-3*eV);
    }
  }
}
