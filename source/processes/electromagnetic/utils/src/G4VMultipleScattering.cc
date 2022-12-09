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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VMultipleScattering
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 25.03.2003
//
// Modifications:
//
// 16-07-03 Use G4VMscModel interface (V.Ivanchenko)
// 03-11-03 Fix initialisation problem in RetrievePhysicsTable (V.Ivanchenko)
// 04-11-03 Update PrintInfoDefinition (V.Ivanchenko)
// 01-03-04 SampleCosineTheta signature changed
// 22-04-04 SampleCosineTheta signature changed back to original
// 27-08-04 Add InitialiseForRun method (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 11-03-05 Shift verbose level by 1 (V.Ivantchenko)
// 15-04-05 optimize internal interface (V.Ivanchenko)
// 15-04-05 remove boundary flag (V.Ivanchenko)
// 27-10-05 introduce virtual function MscStepLimitation() (V.Ivanchenko)
// 12-04-07 Add verbosity at destruction (V.Ivanchenko)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
// 11-03-08 Set skin value does not effect step limit type (V.Ivanchenko)
// 24-06-09 Removed hidden bin in G4PhysicsVector (V.Ivanchenko)
// 04-06-13 Adoptation to MT mode (V.Ivanchenko)
//

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VMultipleScattering.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4EmTableUtil.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::G4VMultipleScattering(const G4String&, G4ProcessType)
  : G4VContinuousDiscreteProcess("msc", fElectromagnetic),
  fNewPosition(0.,0.,0.),
  fNewDirection(0.,0.,1.)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);
  SetProcessSubType(fMultipleScattering);

  lowestKinEnergy = 10*CLHEP::eV;

  geomMin = 0.05*CLHEP::nm;
  minDisplacement2 = geomMin*geomMin;

  pParticleChange = &fParticleChange;

  modelManager = new G4EmModelManager();
  emManager = G4LossTableManager::Instance();
  mscModels.reserve(2);
  emManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::~G4VMultipleScattering()
{
  delete modelManager;
  emManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::AddEmModel(G4int order, G4VMscModel* ptr,
                                       const G4Region* region)
{
  if(nullptr == ptr) { return; }
  G4VEmFluctuationModel* fm = nullptr;
  modelManager->AddEmModel(order, ptr, fm, region);
  ptr->SetParticleChange(pParticleChange);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::SetEmModel(G4VMscModel* ptr, G4int)
{
  if(nullptr == ptr) { return; }
  if(!mscModels.empty()) { 
    for(auto & msc : mscModels) { if(msc == ptr) { return; } } 
  }
  mscModels.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VMultipleScattering::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  G4bool master = emManager->IsMaster();
  if(nullptr == firstParticle) { firstParticle = &part; }

  emManager->PreparePhysicsTable(&part, this, master);
  currParticle = nullptr;

  if(firstParticle == &part) {
    baseMat = emManager->GetTableBuilder()->GetBaseMaterialFlag();

    G4EmTableUtil::PrepareMscProcess(this, part, modelManager,
                                     stepLimit, facrange,
                                     latDisplacement, master,
                                     isIon, baseMat);

    numberOfModels = modelManager->NumberOfModels();
    currentModel = GetModelByIndex(0);

    if(nullptr == safetyHelper) {
      safetyHelper = G4TransportationManager::GetTransportationManager()
        ->GetSafetyHelper();
      safetyHelper->InitialiseHelper();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  G4bool master = emManager->IsMaster();

  if(firstParticle == &part) { 
    emManager->BuildPhysicsTable(firstParticle);
  }
  const G4VMultipleScattering* ptr = this;
  if(!master) {
    ptr = static_cast<const G4VMultipleScattering*>(GetMasterProcess());
  }

  G4EmTableUtil::BuildMscProcess(this, ptr, part, firstParticle,
				 numberOfModels, master);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::StreamInfo(std::ostream& outFile, 
                  const G4ParticleDefinition& part, G4bool rst) const
{
  G4String indent = (rst ? "  " : "");
  outFile << G4endl << indent << GetProcessName() << ": ";
  if (!rst) outFile << " for " << part.GetParticleName();
  outFile  << "  SubType= " << GetProcessSubType() << G4endl;
  modelManager->DumpModelList(outFile, verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::StartTracking(G4Track* track)
{
  G4VEnergyLossProcess* eloss = nullptr;
  if(track->GetParticleDefinition() != currParticle) {
    currParticle = track->GetParticleDefinition();
    fIonisation = emManager->GetEnergyLossProcess(currParticle);
    eloss = fIonisation;
  }
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VMscModel* msc = GetModelByIndex(i);
    msc->StartTracking(track);
    if(nullptr != eloss) {
      msc->SetIonisation(eloss, currParticle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMultipleScattering::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double,
                             G4double currentMinimalStep,
                             G4double&,
                             G4GPILSelection* selection)
{
  // get Step limit proposed by the process
  *selection = NotCandidateForSelection;
  physStepLimit = gPathLength = tPathLength = currentMinimalStep;

  G4double ekin = track.GetKineticEnergy();
  /*
  G4cout << "MSC::AlongStepGPIL: Ekin= " << ekin
         << "  " << currParticle->GetParticleName() 
         << " currMod " << currentModel 
         << G4endl;
  */
  // isIon flag is used only to select a model
  if(isIon) { 
    ekin *= proton_mass_c2/track.GetParticleDefinition()->GetPDGMass(); 
  }
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  
  // select new model, static cast is possible in this class
  if(1 < numberOfModels) {
    currentModel =
      static_cast<G4VMscModel*>(SelectModel(ekin,couple->GetIndex()));
  }
  currentModel->SetCurrentCouple(couple);
  // msc is active is model is active, energy above the limit,
  // and step size is above the limit;
  // if it is active msc may limit the step
  if(currentModel->IsActive(ekin) && tPathLength > geomMin
     && ekin >= lowestKinEnergy) {
    isActive = true;
    tPathLength = 
      currentModel->ComputeTruePathLengthLimit(track, gPathLength);
    if (tPathLength < physStepLimit) { 
      *selection = CandidateForSelection; 
    }
  } else { 
    isActive = false;
    gPathLength = DBL_MAX;
  }
  
  //if(currParticle->GetPDGMass() > GeV)    
  /*
  G4cout << "MSC::AlongStepGPIL: Ekin= " << ekin
         << "  " << currParticle->GetParticleName()
         << " gPathLength= " << gPathLength
         << " tPathLength= " << tPathLength
         << " currentMinimalStep= " << currentMinimalStep
         << " isActive " << isActive << G4endl;
  */
  return gPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VMultipleScattering::PostStepGetPhysicalInteractionLength(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4VMultipleScattering::AlongStepDoIt(const G4Track& track, const G4Step& step)
{
  fParticleChange.InitialiseMSC(track, step);
  fNewPosition = fParticleChange.GetProposedPosition();
  fPositionChanged = false;

  G4double geomLength = step.GetStepLength();

  // very small step - no msc
  if(!isActive) {
    tPathLength = geomLength;

    // sample msc
  } else {
    G4double range = 
      currentModel->GetRange(currParticle,track.GetKineticEnergy(),
                             track.GetMaterialCutsCouple());

    tPathLength = currentModel->ComputeTrueStepLength(geomLength);
  
    /*    
    if(currParticle->GetPDGMass() > 0.9*GeV)    
    G4cout << "G4VMsc::AlongStepDoIt: GeomLength= " 
           << geomLength 
           << " tPathLength= " << tPathLength
           << " physStepLimit= " << physStepLimit
           << " dr= " << range - tPathLength
           << " ekin= " << track.GetKineticEnergy() << G4endl;
    */
    // protection against wrong t->g->t conversion
    tPathLength = std::min(tPathLength, physStepLimit);

    // do not sample scattering at the last or at a small step
    if(tPathLength < range && tPathLength > geomMin) {

      static const G4double minSafety = 1.20*CLHEP::nm;
      static const G4double sFact = 0.99;

      G4ThreeVector displacement = currentModel->SampleScattering(
        step.GetPostStepPoint()->GetMomentumDirection(),minSafety);

      G4double r2 = displacement.mag2();
      //G4cout << "    R= " << sqrt(r2) << " Rmin= " << sqrt(minDisplacement2)
      //     << " flag= " << fDispBeyondSafety << G4endl;
      if(r2 > minDisplacement2) {

        fPositionChanged = true;
        G4double dispR = std::sqrt(r2);
        G4double postSafety = 
          sFact*safetyHelper->ComputeSafety(fNewPosition, dispR); 
        //G4cout<<"    R= "<< dispR<<" postSafety= "<<postSafety<<G4endl;

        // far away from geometry boundary
        if(postSafety > 0.0 && dispR <= postSafety) {
          fNewPosition += displacement;

          //near the boundary
        } else {
          // displaced point is definitely within the volume
          //G4cout<<"    R= "<<dispR<<" postSafety= "<<postSafety<<G4endl;
          if(dispR < postSafety) {
            fNewPosition += displacement;

            // reduced displacement
          } else if(postSafety > geomMin) {
            fNewPosition += displacement*(postSafety/dispR); 

            // very small postSafety
          } else {
            fPositionChanged = false;
          }
        }
        if(fPositionChanged) { 
          safetyHelper->ReLocateWithinVolume(fNewPosition);
          fParticleChange.ProposePosition(fNewPosition); 
        }
      }
    }
  }
  fParticleChange.ProposeTrueStepLength(tPathLength);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::GetContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  G4GPILSelection selection = NotCandidateForSelection;
  G4double x = AlongStepGetPhysicalInteractionLength(track,previousStepSize,
                                                     currentMinimalStep,
                                                     currentSafety, 
                                                     &selection);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::ContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  return GetContinuousStepLimit(track,previousStepSize,currentMinimalStep,
                                currentSafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::GetMeanFreePath(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4VMultipleScattering::StorePhysicsTable(const G4ParticleDefinition* part,
                                         const G4String& directory,
                                         G4bool ascii)
{
  G4bool yes = true;
  if(part != firstParticle || !emManager->IsMaster()) { return yes; }

  return G4EmTableUtil::StoreMscTable(this, part, directory,
				      numberOfModels, verboseLevel,
                                      ascii);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4VMultipleScattering::RetrievePhysicsTable(const G4ParticleDefinition*,
                                            const G4String&,
                                            G4bool)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::ProcessDescription(std::ostream& outFile) const
{
  if(nullptr != firstParticle) {
    StreamInfo(outFile, *firstParticle, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

