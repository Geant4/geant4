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
// $Id: G4VMultipleScattering.cc 107364 2017-11-09 10:53:25Z gcosmo $
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
// 13.04.03 Change printout (V.Ivanchenko)
// 04-06-03 Fix compilation warnings (V.Ivanchenko)
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
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles

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
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::G4VMultipleScattering(const G4String& name, G4ProcessType)
  : G4VContinuousDiscreteProcess("msc", fElectromagnetic),
  numberOfModels(0),
  firstParticle(nullptr),
  currParticle(nullptr),
  stepLimit(fUseSafety),
  facrange(0.04),
  latDisplacement(true),
  isIon(false),
  fDispBeyondSafety(false),
  fNewPosition(0.,0.,0.),
  fNewDirection(0.,0.,1.)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);
  SetProcessSubType(fMultipleScattering);
  if("ionmsc" == name) { firstParticle = G4GenericIon::GenericIon(); }

  lowestKinEnergy = 10*CLHEP::eV;

  physStepLimit = gPathLength = tPathLength = 0.0;
  fIonisation = nullptr;

  geomMin   = 0.05*CLHEP::nm;
  minDisplacement2 = geomMin*geomMin;

  pParticleChange = &fParticleChange;
  safetyHelper = nullptr;
  fPositionChanged = false;
  isActive = false;

  currentModel = nullptr;
  modelManager = new G4EmModelManager();
  emManager = G4LossTableManager::Instance();
  emManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::~G4VMultipleScattering()
{
  /*
  if(1 < verboseLevel) {
    G4cout << "G4VMultipleScattering destruct " << GetProcessName() 
          << G4endl;
  }
  */
  delete modelManager;
  emManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::AddEmModel(G4int order, G4VEmModel* p,
                                       const G4Region* region)
{
  G4VEmFluctuationModel* fm = nullptr;
  modelManager->AddEmModel(order, p, fm, region);
  if(p) { p->SetParticleChange(pParticleChange); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::SetEmModel(G4VMscModel* p, size_t)
{
  for(auto & msc : mscModels) { if(msc == p) { return; } }
  mscModels.push_back(p);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMscModel* G4VMultipleScattering::EmModel(size_t index) const
{
  return (index < mscModels.size()) ? mscModels[index] : nullptr; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VMultipleScattering::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
  G4bool master = emManager->IsMaster();

  if(!firstParticle) { firstParticle = &part; }
  if(part.GetParticleType() == "nucleus") {
    stepLimit = fMinimal;
    latDisplacement = false;
    facrange = 0.2;
    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "alpha"    && pname != "He3" &&
       pname != "hydrogen") {

      const G4ParticleDefinition* theGenericIon = 
        G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");
      if(&part == theGenericIon) { isIon = true; }
     
      if(theGenericIon && firstParticle != theGenericIon) {
        G4ProcessManager* pm =  theGenericIon->GetProcessManager();
        G4ProcessVector* v = pm->GetAlongStepProcessVector();
        size_t n = v->size();
        for(size_t j=0; j<n; ++j) {
          if((*v)[j] == this) {
            firstParticle = theGenericIon;
            isIon = true; 
            break; 
          }
        }
      }
    }
  }

  emManager->PreparePhysicsTable(&part, this, master);
  currParticle = nullptr;

  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << " local particle " << firstParticle->GetParticleName()
           << " isIon: " << isIon << " isMaster: " << master
           << G4endl;
  }

  if(firstParticle == &part) {

    // initialise process
    InitialiseProcess(firstParticle);

    // heavy particles and not ions
    if(!isIon) {
      if(part.GetPDGMass() > MeV) {
        stepLimit = theParameters->MscMuHadStepLimitType(); 
        facrange = theParameters->MscMuHadRangeFactor(); 
        latDisplacement = theParameters->MuHadLateralDisplacement();
      } else {
        stepLimit = theParameters->MscStepLimitType(); 
        facrange = theParameters->MscRangeFactor(); 
        latDisplacement = theParameters->LateralDisplacement();
      }
      if(latDisplacement) { 
        fDispBeyondSafety = theParameters->LatDisplacementBeyondSafety();
      }
    }
    if(master) { SetVerboseLevel(theParameters->Verbose()); }
    else {  SetVerboseLevel(theParameters->WorkerVerbose()); }

    // initialisation of models
    numberOfModels = modelManager->NumberOfModels();
    /*
        G4cout << "### G4VMultipleScattering::PreparePhysicsTable() for "
        << GetProcessName()
        << " and particle " << part.GetParticleName()
        << " Nmod= " << numberOfModels << "  " << this
        << G4endl;
    */
    for(G4int i=0; i<numberOfModels; ++i) {
      G4VMscModel* msc = static_cast<G4VMscModel*>(GetModelByIndex(i));
      if(!msc) { continue; }
      msc->SetIonisation(nullptr, firstParticle);
      msc->SetMasterThread(master);
      currentModel = msc; 
      msc->SetStepLimitType(stepLimit);
      msc->SetLateralDisplasmentFlag(latDisplacement);
      msc->SetSkin(theParameters->MscSkin());
      msc->SetRangeFactor(facrange);
      msc->SetGeomFactor(theParameters->MscGeomFactor());
      msc->SetPolarAngleLimit(theParameters->MscThetaLimit());
      G4double emax = 
        std::min(msc->HighEnergyLimit(),theParameters->MaxKinEnergy());
      msc->SetHighEnergyLimit(emax);
    }

    modelManager->Initialise(firstParticle, G4Electron::Electron(), 
                             10.0, verboseLevel);

    if(!safetyHelper) {
      safetyHelper = G4TransportationManager::GetTransportationManager()
        ->GetSafetyHelper();
      safetyHelper->InitialiseHelper();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  G4String num = part.GetParticleName();
  G4bool master = emManager->IsMaster();
  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << num << " isIon: " << isIon
           << " IsMaster: " << master << G4endl;
  }
  const G4VMultipleScattering* masterProcess = 
    static_cast<const G4VMultipleScattering*>(GetMasterProcess());

  if(firstParticle == &part) { 
    /*    
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << num
           << " IsMaster= " << G4LossTableManager::Instance()->IsMaster()
           << "  " << this
           << G4endl;
    */
    emManager->BuildPhysicsTable(firstParticle);

    if(!master) {
      // initialisation of models
      /*
        G4cout << "### G4VMultipleScattering::BuildPhysicsTable() for "
        << GetProcessName()
        << " and particle " << num
        << " Nmod= " << numberOfModels << "  " << this
        << G4endl;
      */
      for(G4int i=0; i<numberOfModels; ++i) {
        G4VMscModel* msc = static_cast<G4VMscModel*>(GetModelByIndex(i));
        if(!msc) { continue; }
        G4VMscModel* msc0= 
	  static_cast<G4VMscModel*>(masterProcess->GetModelByIndex(i));
        msc->SetCrossSectionTable(msc0->GetCrossSectionTable(), false);
        msc->InitialiseLocal(firstParticle, msc0);
      }
    }
  }

  // explicitly defined printout by particle name
  if(1 < verboseLevel || 
     (0 < verboseLevel && (num == "e-" || 
                           num == "e+"    || num == "mu+" || 
                           num == "mu-"   || num == "proton"|| 
                           num == "pi+"   || num == "pi-" || 
                           num == "kaon+" || num == "kaon-" || 
                           num == "alpha" || num == "anti_proton" || 
                           num == "GenericIon" || num == "alpha+" || 
			   num == "alpha++" )))
    { 
      StreamInfo(G4cout, part);
    }

  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << num
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::StreamInfo(std::ostream& outFile, 
                  const G4ParticleDefinition& part, G4String endOfLine) const
{
  outFile << endOfLine << GetProcessName() << ": ";
	if (endOfLine != G4String("<br>\n")) {
		outFile << "  for " << part.GetParticleName();
	}
	outFile  << "    SubType= " << GetProcessSubType() << endOfLine;
  StreamProcessInfo(outFile, endOfLine);
  modelManager->DumpModelList(outFile, verboseLevel, endOfLine);
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
  /*
  G4cout << "G4VMultipleScattering::StartTracking Nmod= " << numberOfModels
         << "  " << currParticle->GetParticleName() 
         << " E(MeV)= " << track->GetKineticEnergy()
         << "  Ion= " << eloss << "  " << fIonisation << " IsMaster= " 
         << G4LossTableManager::Instance()->IsMaster() 
         << G4endl;
  */
  for(G4int i=0; i<numberOfModels; ++i) {
    /*
    G4cout << "Next model " << i << "  " << msc 
	   << " Emin= " << msc->LowEnergyLimit() 
	   << " Emax= " << msc->HighEnergyLimit() 
	   << " Eact= " << msc->LowEnergyActivationLimit() << G4endl;
    */
    G4VEmModel* msc = GetModelByIndex(i);
    msc->StartTracking(track);
    if(eloss) { 
      G4VMscModel* mscmod = static_cast<G4VMscModel*>(msc);
      if(mscmod) { mscmod->SetIonisation(fIonisation, currParticle); } 
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
  
  // select new model
  if(1 < numberOfModels) {
    currentModel = static_cast<G4VMscModel*>(
      SelectModel(ekin,track.GetMaterialCutsCouple()->GetIndex()));
  }
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
  } else { isActive = false; }
  
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
  fParticleChange.ProposeMomentumDirection(
    step.GetPostStepPoint()->GetMomentumDirection());
  fNewPosition = step.GetPostStepPoint()->GetPosition();
  fParticleChange.ProposePosition(fNewPosition);
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

            // optional extra mechanism is applied only if a particle
            // is stopped by the boundary
          } else if(fDispBeyondSafety && 0.0 == postSafety) {
            fNewPosition += displacement;
            G4double maxshift = 
              std::min(2.0*dispR, geomLength*(physStepLimit/tPathLength - 1.0));
            G4double dist = 0.0;
            G4double safety = postSafety + dispR;
            fNewDirection = *(fParticleChange.GetMomentumDirection());
            /*
              G4cout << "##MSC before Recheck maxshift= " << maxshift
                     << " postsafety= " << postSafety
                     << " Ekin= " << track.GetKineticEnergy()
                     << "  " << track.GetDefinition()->GetParticleName()
                     << G4endl; 
            */
            // check if it is possible to shift to the boundary
	    // and the shift is not large
            if(safetyHelper->RecheckDistanceToCurrentBoundary(fNewPosition,
               fNewDirection, maxshift, &dist, &safety) 
	       && std::abs(dist) < maxshift) {
	         /* 
		    G4cout << "##MSC after Recheck dist= " << dist
		    << " postsafety= " << postSafety
		    << " t= " << tPathLength
		    << " g=  " << geomLength
		    << " p=  " << physStepLimit
		    << G4endl; 
	         */
	      // shift is positive
	      if(dist >= 0.0) {
		tPathLength *= (1.0 + dist/geomLength);
		fNewPosition += dist*fNewDirection; 

		// shift is negative cannot be larger than geomLength
	      } else {
		maxshift = std::min(maxshift, geomLength);
		if(0.0 < maxshift + dist) {
		  const G4ThreeVector& postpoint = step.GetPostStepPoint()->GetPosition();
		  G4ThreeVector point = fNewPosition + dist*fNewDirection;
		  G4double R2 = (postpoint - point).mag2();
		  G4double newdist = dist;
		  // check not more than 10 extra boundaries
		  for(G4int i=0; i<10; ++i) {
		    dist = 0.0;
		    if(safetyHelper->RecheckDistanceToCurrentBoundary(
		       point, fNewDirection, maxshift, &dist, &safety) 
		       && std::abs(newdist + dist) < maxshift) {
		      point += dist*fNewDirection;
		      G4double R2new = (postpoint - point).mag2();
		      //G4cout << "Backward  i= " << i << " dist= " << dist 
		      //     << " R2= " << R2new << G4endl; 
		      if(dist >= 0.0 || R2new > R2) { break; }
		      R2 = R2new;
		      fNewPosition = point;
		      newdist += dist;
		    } else { 
		      break; 
		    }
		  }
		  tPathLength *= (1.0 + newdist/geomLength);
		// shift on boundary is not possible for negative disp
		} else {
		  fNewPosition += displacement*(postSafety/dispR - 1.0); 
		}
	      }
	      // shift on boundary is not possible for any disp
	    } else {
	      fNewPosition += displacement*(postSafety/dispR - 1.0); 
	    }
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

G4VParticleChange* 
G4VMultipleScattering::PostStepDoIt(const G4Track& track, const G4Step&)
{
  fParticleChange.Initialize(track);
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
  if(part != firstParticle) { return yes; }
  const G4VMultipleScattering* masterProcess = 
    static_cast<const G4VMultipleScattering*>(GetMasterProcess()); 
  if(masterProcess && masterProcess != this) { return yes; }

  G4int nmod = modelManager->NumberOfModels();
  static const G4String ss[4] = {"1","2","3","4"};
  for(G4int i=0; i<nmod; ++i) {
    G4VEmModel* msc = modelManager->GetModel(i);
    yes = true;
    G4PhysicsTable* table = msc->GetCrossSectionTable();
    if (table) {
      G4int j = std::min(i,3); 
      G4String name = 
        GetPhysicsTableFileName(part,directory,"LambdaMod"+ss[j],ascii);
      yes = table->StorePhysicsTable(name,ascii);

      if ( yes ) {
        if ( verboseLevel>0 ) {
          G4cout << "Physics table are stored for " 
                 << part->GetParticleName()
                 << " and process " << GetProcessName()
                 << " with a name <" << name << "> " << G4endl;
        }
      } else {
        G4cout << "Fail to store Physics Table for " 
               << part->GetParticleName()
               << " and process " << GetProcessName()
               << " in the directory <" << directory
               << "> " << G4endl;
      }
    }
  }
  return yes;
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

void G4VMultipleScattering::SetIonisation(G4VEnergyLossProcess* p)
{
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VMscModel* msc = static_cast<G4VMscModel*>(GetModelByIndex(i, true));
    if(msc) { msc->SetIonisation(p, firstParticle); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::ProcessDescription(std::ostream& outFile) const
{
  if(firstParticle) {
    StreamInfo(outFile, *firstParticle, G4String("<br>\n"));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

