//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEnergyLossSTD
//
// Author:        Vladimir Ivanchenko 
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// Class Description: 
//
// It is the unified energy loss process it calculates the continuous 
// energy loss for charged particles using a set of Energy Loss
// models valid for different energy regions. There are a possibility
// to create and access to dE/dx and range tables, or to calculate
// that information on fly.
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
#include "G4VEnergyLossSTD.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4VSubCutoffProcessor.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VEnergyLossSTD::G4VEnergyLossSTD(const G4String& name, G4ProcessType type):
                 G4VContinuousDiscreteProcess(name, type),
  subCutoffProcessor(0),
  emFluctModel(0),
  theDEDXTable(0),
  theRangeTable(0),
  theSecondaryRangeTable(0),
  theInverseRangeTable(0),
  theLambdaTable(0),
  minKinEnergy(1.0*eV),
  maxKinEnergy(100.0*GeV),
  particle(0),
  baseParticle(0),
  secondaryParticle(0),
  theGamma(G4Gamma::Gamma()),
  theElectron(G4Electron::Electron()),
  lossFluctuationFlag(false),
  subCutoffFlag(false),
  subCutoffIsDesired(false),  
  rndmStepFlag(false),
  hasRestProcess(true),
  tablesAreBuilt(false),
  integral(true),
  linLossLimit(0.05),
  minSubRange(0.1)
{
  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
  
  // default dRoverRange and finalRange
  SetStepLimits(0.2, 200.0*micrometer);

  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossSTD::~G4VEnergyLossSTD()
{
  Clear();
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  if(subCutoffProcessor) delete subCutoffProcessor;  
  delete modelManager;

  if(emFluctModel) delete emFluctModel;
  (G4LossTableManager::Instance())->Clear();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::Clear()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Clear() for " << GetProcessName() << G4endl;
  }
  theDEDXTable = 0;
  theRangeTable = 0;
  theInverseRangeTable = 0;
  theSecondaryRangeTable = 0;
  modelManager->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::Initialise()
{

  Clear();

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Initialise() for " 
           << GetProcessName() << G4endl;
  }

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass  = particle->GetPDGMass();
  chargeSquare = initialCharge*initialCharge/(eplus*eplus);
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;

  if(particle->GetProcessManager()->GetAtRestProcessVector()->size()) 
               hasRestProcess = true;
  else         hasRestProcess = false;

  if (baseParticle) {
    massRatio = (baseParticle->GetPDGMass())/initialMass;
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    reduceFactor = 1.0/(chargeSqRatio*massRatio);
  }

  theCuts = modelManager->Initialise(particle, secondaryParticle, minSubRange, verboseLevel);

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Initialise(): emFluctModel= " << emFluctModel
           << " subCutoffFlag= " << subCutoffFlag
           << G4endl;
  }  

  if(emFluctModel) emFluctModel->Initialise(particle);

  if (subCutoffFlag) {
    subCutoffProcessor->Initialise(particle, secondaryParticle, modelManager);
  }
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Initialise() is done " 
           << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio 
           << " reduceFactor= " << reduceFactor
           << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::BuildPhysicsTable(const G4ParticleDefinition& part)
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildPhysicsTable() for " 
           << GetProcessName() << G4endl;
  }

  tablesAreBuilt = false;

  // Are particle defined?
  if(!particle) {
    particle = &part; 
  }
  if(particle != &part) {
    G4cout << "G4VEnergyLossSTD::BuildPhysicsTable: " 
           << particle->GetParticleName() 
           << " != " << part.GetParticleName() 
           << G4endl;
    G4Exception("G4VEnergyLossSTD::BuildPhysicsTable with wrong G4ParticleDefinition"); 
  }

  if(!baseParticle) baseParticle = DefineBaseParticle(particle);

  if(particle->GetParticleType() == "nucleus" && 
     particle->GetParticleName() != "GenericIon" &&
     theLambdaTable) return;


  // It is responsability of the G4LossTables to build DEDX and range tables
  (G4LossTableManager::Instance())->BuildDEDXTable(particle);

  // Access to materials
  //const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nMaterials = G4Material::GetNumberOfMaterials();

  //!!!! To be checked because finalRange is changed here !!!!
  for (size_t i=0; i<nMaterials; i++) {
      
    G4double lengthCut = (G4Electron::Electron()->GetLengthCuts())[i];
    if (finalRange > lengthCut) finalRange = lengthCut;
  }

  if(theLambdaTable) theLambdaTable->clearAndDestroy();

  theLambdaTable = BuildLambdaTable();

  if (subCutoffFlag) {
    subCutoffProcessor->SetLambdaSubTable(BuildLambdaSubTable());
  }
  tablesAreBuilt = true;

  //  if(particle == G4Electron::Electron() || particle == G4Proton::Proton()
  //   || particle == G4GenericIon::GenericIon()) 
  if(!baseParticle) PrintInfoDefinition();

  if(0 < verboseLevel) {
    G4cout << "Tables are built for " << particle->GetParticleName()
           << " Integral= " <<  integral
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl; 
      if(theDEDXTable) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "RangeTable address= " << theRangeTable << G4endl; 
      if(theRangeTable) G4cout << (*theRangeTable) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable << G4endl; 
      if(theInverseRangeTable) G4cout << (*theInverseRangeTable) << G4endl;
    }
  }
 
  if(2 < verboseLevel) {
    G4cout << "LambdaTable address= " << theLambdaTable << G4endl; 
    if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossSTD::StorePhysicsTable(G4ParticleDefinition* part,
			 	     const G4String& directory, 
				           G4bool ascii)
{
  G4String filename;
  
  if (!baseParticle) {
 
    G4String f1 = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
    G4String f2 = GetPhysicsTableFileName(part,directory,"Range",ascii);
    G4String f3 = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);

    if(!((G4LossTableManager::Instance())->StorePhysicsTable(this,f1,f2,f3,ascii)))
      {
        return false;
      }
  }

  if (theLambdaTable) {
    filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    if ( !theLambdaTable->StorePhysicsTable(filename, ascii) ){
      G4cout << "Fatal error theLambdaTable->StorePhysicsTable in <" 
             << filename << ">"
             << G4endl;
      return false;
    }
  }

  G4PhysicsTable*  theLambdaSubTable = 0;
  if(subCutoffProcessor) theLambdaSubTable = subCutoffProcessor->LambdaSubTable();

  if (theLambdaSubTable) {
    filename = GetPhysicsTableFileName(part,directory,"LambdaSub",ascii);
    if ( !theLambdaSubTable->StorePhysicsTable(filename, ascii) ){
      G4cout << "Fatal error theLambdaSubTable->StorePhysicsTable in <" 
             << filename << ">"
             << G4endl;
      return false;
    }
  }
  G4cout << GetProcessName() << " for " << particle->GetParticleName() 
         << ": Success to store the PhysicsTables in "  
         << directory << G4endl;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEnergyLossSTD::RetrievePhysicsTable(G4ParticleDefinition* part,
				       const G4String& directory, 
				             G4bool ascii)
{

  // Are particle defined?
  if(!particle) {
    particle = part; 
  }
  if(particle != part) {
    G4cout << "G4VEnergyLossSTD::RetrievePhysicsTable: " 
           << particle->GetParticleName() 
           << " != " << part->GetParticleName() 
           << G4endl;
    G4Exception("G4VEnergyLossSTD::RetrievePhysicsTable with wrong G4ParticleDefinition"); 
  }

  Initialise();

  G4PhysicsTable* theTable;
  G4String filename;
  G4String particleName = particle->GetParticleName();
  G4int nMaterials = G4Material::GetNumberOfMaterials();

  filename = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
  theTable = new G4PhysicsTable(nMaterials);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "DEDX table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;  
  } else {
    G4cout << "DEDX table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;  
  }
  theDEDXTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"Range",ascii);
  theTable = new G4PhysicsTable(nMaterials);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Range table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;  
  } else {
    G4cout << "Range table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;  
  }
  theRangeTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
  theTable = new G4PhysicsTable(nMaterials);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Inverse table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;  
    return false;
  }
  theInverseRangeTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  theLambdaTable = new G4PhysicsTable(nMaterials);
  if (theLambdaTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Lambda table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;  
  } else {
    G4cout << "Lambda table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;  
  }

  G4PhysicsTable*  theLambdaSubTable = 0;
  filename = GetPhysicsTableFileName(part,directory,"LambdaSub",ascii);
  theLambdaSubTable = new G4PhysicsTable(nMaterials);
  if (theLambdaSubTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "LambdaSub table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;  
    if(subCutoffProcessor) subCutoffProcessor->SetLambdaSubTable(theLambdaSubTable);
  } else {
    G4cout << "LambdaSub table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;  
  }

  // G4LossTableManager manages DEDX and range tables

  (G4LossTableManager::Instance())->RetrieveDEDXTable(particle, this);

  
  G4cout << GetProcessName() << " for " << particleName 
         << ": end of retrieving PhysicsTables from directory <"
         << directory << ">" << G4endl;
 	 
  return true;
}
  
		
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetParticles(const G4ParticleDefinition* p1,
                                    const G4ParticleDefinition* p2,
                                    const G4ParticleDefinition* p3)
{
  particle = p1;
  baseParticle = p2;
  secondaryParticle = p3;
  G4bool yes = false;
  if(particle) yes = particle->GetPDGCharge();
  if(baseParticle) {
    if(!baseParticle->GetPDGCharge()) yes = false;
    else {
      if((particle->GetPDGMass() < MeV &&
          baseParticle->GetPDGMass() > MeV) ||
         (particle->GetPDGMass() > MeV &&
          baseParticle->GetPDGMass() < MeV)) yes = false;
    } 
  }
  if(!yes) {
    G4cout << "Warning in G4VEnergyLossSTD::SetParticle: "
           << particle->GetParticleName() 
           << " losses cannot be obtained from "
           << baseParticle->GetParticleName()
           << " losses" << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::AddEmModel(G4VEmModel* p, G4int order)
{
  modelManager->AddEmModel(p, order);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::AddEmFluctuationModel(G4VEmFluctuationModel* p)
{
  if(emFluctModel) delete emFluctModel;
  emFluctModel = p;
  if(p) lossFluctuationFlag = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildDEDXTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildDEDXTable() for " 
           << GetProcessName()
           << " and particle " << particle->GetParticleName() 
           << G4endl;
  }

  // vectors to provide continues dE/dx
  G4DataVector factor;
  G4DataVector dedxLow;
  G4DataVector dedxHigh;

  // Access to materials
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nMaterials = G4Material::GetNumberOfMaterials();
  G4PhysicsTable* theTable = new G4PhysicsTable(nMaterials);

  if(0 < verboseLevel) {
    G4cout << nMaterials << " materials"
           << " minKinEnergy= " << minKinEnergy  
           << " maxKinEnergy= " << maxKinEnergy  
           << G4endl;
  }

  for(size_t i=0; i<nMaterials; i++) {

    // create physics vector and fill it  
    const G4Material* material = (*theMaterialTable)[i];
    G4PhysicsVector* aVector = DEDXPhysicsVector(material);
    modelManager->FillDEDXVector(aVector, material);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildDEDXTable(): table is built" << G4endl;
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildLambdaTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaTable() for process " 
           << GetProcessName() << " and particle " 
           << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  G4PhysicsTable* theTable = new G4PhysicsTable(nMaterials);

  for(G4int i=0; i<nMaterials; i++) {

    // create physics vector and fill it  
    const G4Material* material = (*theMaterialTable)[i];
    G4PhysicsVector* aVector = LambdaPhysicsVector(material);
    modelManager->FillLambdaVector(aVector, material);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built" << G4endl;
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildLambdaSubTable()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaSubTable() for process " 
           << GetProcessName() << " and particle " 
           << particle->GetParticleName() << G4endl;
  }

  // Access to materials
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nMaterials = G4Material::GetNumberOfMaterials();
  G4PhysicsTable* theTable = new G4PhysicsTable(nMaterials);

  for(size_t i=0; i<nMaterials; i++) {

    // create physics vector and fill it  
    const G4Material* material = (*theMaterialTable)[i];
    G4PhysicsVector* aVector = SubLambdaPhysicsVector(material);
    modelManager->FillSubLambdaVector(aVector, material);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(1 < verboseLevel) {
    G4cout << "Table is built" << G4endl;
  }

  return theTable;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossSTD::AlongStepDoIt(const G4Track& track,
                                                   const G4Step& step)
{
  aParticleChange.Initialize(track);
  // The process has range table - calculate energy loss
  if(!theRangeTable) return &aParticleChange;

  // Get the actual (true) Step length  
  G4double length = step.GetStepLength();
  G4double eloss  = 0.0;
  G4bool b;
  G4double finalT = preStepKinEnergy;

  /*
  if(-1 < verboseLevel) {
    G4cout << "AlongStepDoIt for " 
           << GetProcessName() << " and particle " 
           << particle->GetParticleName() 
           << "  eScaled(MeV)= " << finalT/MeV
           << "  slim(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << G4endl;
  }
  */

  static const G4double faclow = 1.5;

    // low energy deposit case
  if (finalT < minKinEnergy || length >= fRange) {
    eloss = finalT;

  } else if(finalT < faclow*minKinEnergy) {

    eloss = finalT*sqrt(length/fRange);

  // Short step
  } else if(length <= linLossLimit * fRange || preStepScaledEnergy > maxKinEnergy) {
    eloss = (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio;

  // Long step
  } else {
    G4double x = (fRange-length)/reduceFactor;
    G4PhysicsVector* v = (*theInverseRangeTable)[currentMaterialIndex];
    // G4double esc = v->GetValue(fRange/reduceFactor, b);
    G4double esp = v->GetValue(x, b);
    eloss = (preStepScaledEnergy - esp)/massRatio;
    
    /*    
    if(-1 < verboseLevel) {
      G4cout << "fRange(mm)= " << fRange/mm
             << " xPost(mm)= " << x/mm
             << " ePre(MeV)= " << preStepScaledEnergy/MeV
             << " ePre0(MeV)= " << esc/MeV
             << " ePost(MeV)= " << esp/MeV
             << " eloss(MeV)= " << eloss/MeV 
             << " eloss0(MeV)= " << (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio
             << G4endl;
    }
    */
    
    if(eloss < 0.0) eloss = 0.0; 
  }    

  // 
  if(eloss > finalT) eloss = finalT;

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  
  // Sub cutoff generation
  if(subCutoffFlag) {

    G4std::vector<G4Track*>* newp = subCutoffProcessor->
           SampleSecondary(step, dynParticle, eloss);

    if(newp) {

      G4ThreeVector finalP = dynParticle->GetMomentum();
      G4int n = newp->size();
      if(n > 0) {
        aParticleChange.SetNumberOfSecondaries(n);
        G4Track* t; 
        G4double e;
        for (G4int i=0; i<n; i++) {
          t = (*newp)[i];
          e = t->GetKineticEnergy();
          const G4ParticleDefinition* pd = t->GetDefinition();
          if (pd != theGamma && pd != theElectron ) e += pd->GetPDGMass();
        
          finalT -= e;
          finalP -= t->GetMomentum();
          aParticleChange.AddSecondary(t);
        }
        aParticleChange.SetMomentumChange(finalP.unit());
      }
      delete newp;
    }
  }
  
    // Sample fluctuations

  if (lossFluctuationFlag && eloss + minKinEnergy > finalT) {

    G4double tmax = MaxSecondaryEnergy(dynParticle);
    tmax = G4std::min(tmax,(*theCuts)[currentMaterialIndex]);

    emFluctModel->SampleFluctuations(currentMaterial, dynParticle, 
                                               tmax, length, eloss);
  }

  /*
  if(1 < verboseLevel) {
    G4cout << "eloss(MeV)= " << eloss/MeV 
           << " currentChargeSquare= " << chargeSquare 
           << G4endl;
  }
  */
  finalT -= eloss;

  //  G4cout << "finalT= " << finalT << G4endl;

  if (finalT < minKinEnergy) {

    eloss += finalT;
    finalT = 0.0;

    // !!!! to be checked - problem of Forced processes

    if (hasRestProcess) aParticleChange.SetStatusChange(fStopButAlive);
    else                aParticleChange.SetStatusChange(fStopAndKill); 
  }

  aParticleChange.SetEnergyChange(finalT);
  aParticleChange.SetLocalEnergyDeposit(eloss);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossSTD::PostStepDoIt(const G4Track& track,
                                                 const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double finalT = track.GetKineticEnergy();

  // Integral approach
  if(integral) {  
    G4bool b;
    G4double postStepLambda = 
      (((*theLambdaTable)[currentMaterialIndex])->GetValue(finalT,b))*chargeSqRatio;
  
    if(preStepLambda*G4UniformRand() > postStepLambda) 
      return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = modelManager->SelectModel(finalT);
  G4double tcut = (*theCuts)[currentMaterialIndex];
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  G4std::vector<G4DynamicParticle*>* newp = 
    currentModel->SampleSecondary(currentMaterial, dynParticle, tcut, finalT);

  if (newp) {
    G4int n = newp->size();
    
    G4ThreeVector finalP = dynParticle->GetMomentum();
    G4StepPoint* postPoint = step.GetPostStepPoint();
    G4ThreeVector point = postPoint->GetPosition();
    G4double time = postPoint->GetGlobalTime();
    G4DynamicParticle* dp;
    G4double e;
    aParticleChange.SetNumberOfSecondaries(n);
    for (G4int i=0; i<n; i++) {
       dp = (*newp)[i];
       e = dp->GetKineticEnergy();
       const G4ParticleDefinition* pd = dp->GetDefinition();
       if (pd != theGamma && pd != theElectron ) e += pd->GetPDGMass();
       finalT -= e;
       finalP -= dp->GetMomentum();
       G4Track* t = new G4Track(dp, time, point);
       aParticleChange.AddSecondary(t); 
    }

    if (finalT < minKinEnergy) {
      aParticleChange.SetLocalEnergyDeposit(finalT);
      finalT = 0.0;

      // !!!! to be checked - problem of Forced processes
 
      if (hasRestProcess) aParticleChange.SetStatusChange(fStopButAlive);
      else                aParticleChange.SetStatusChange(fStopAndKill); 
    }

    aParticleChange.SetEnergyChange(finalT);
    aParticleChange.SetMomentumChange(finalP.unit());

    delete newp;
  }

  return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::PrintInfoDefinition() const
{
  G4cout << G4endl << GetProcessName() << ":  " << G4endl
         << "      dE/dx and range tables from "
	 << G4BestUnit(MinKinEnergy(),"Energy")
         << " to " << G4BestUnit(MaxKinEnergy(),"Energy") 
         << " in " << DEDXBinning() << " bins." << G4endl
         << "      Lambda tables from threshold to "
         << G4BestUnit(MaxKinEnergy(),"Energy") 
         << " in " << LambdaBinning() << " bins." 
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetDEDXTable(G4PhysicsTable* p)
{
  if(theDEDXTable) delete theDEDXTable; 
  theDEDXTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetRangeTable(G4PhysicsTable* p)
{
  if(theRangeTable) delete theRangeTable; 
  theRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  if(theSecondaryRangeTable) delete theSecondaryRangeTable; 
  theSecondaryRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetInverseRangeTable(G4PhysicsTable* p)
{
  if(theInverseRangeTable) delete theInverseRangeTable; 
  theInverseRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::DEDXPhysicsVector(const G4Material*)
{
  return new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nDEDXBins);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::LambdaPhysicsVector(const G4Material* material)
{
  G4PhysicsVector* v = 0;
  G4double cut  = (*theCuts)[material->GetIndex()];
  if(cut > maxKinEnergy) {
    v = new G4PhysicsLogVector(maxKinEnergy*0.5, maxKinEnergy, 2);
  } else {
    G4double tmin = G4std::max(MinPrimaryEnergy(particle, material, cut), minKinEnergy);
    if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
    v = new G4PhysicsLogVector(tmin, maxKinEnergy, nLambdaBins);
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::SubLambdaPhysicsVector(const G4Material* material)
{
  G4PhysicsVector* v = 0;
  G4double cut  = (*theCuts)[material->GetIndex()];
  if(cut > maxKinEnergy) {
    v = new G4PhysicsLogVector(maxKinEnergy*0.5, maxKinEnergy, 2);
  } else {
    G4double tmin = G4std::max(MinPrimaryEnergy(particle, material, cut), minKinEnergy);
    if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
    v = new G4PhysicsLogVector(tmin, maxKinEnergy, nLambdaBins);
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetSubCutoffProcessor(G4VSubCutoffProcessor* p)
{
  if(subCutoffProcessor) delete subCutoffProcessor;
  subCutoffProcessor = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossSTD::MicroscopicCrossSection(G4double kineticEnergy, 
                                             const G4Material* material)
{ 
  // Cross section per atom is calculated
  DefineMaterial(material);
  G4double cross = 0.0;
  G4bool b;
  if(theLambdaTable) {
    cross = (((*theLambdaTable)[currentMaterialIndex])->
                           GetValue(kineticEnergy, b));

    cross /= material->GetTotNbOfAtomsPerVolume();
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4double G4VEnergyLossSTD::MeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  return GetMeanFreePath(track, s, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossSTD::ContinuousStepLimit(const G4Track& track,
                                               G4double x, G4double y, G4double& z)
{
  return GetContinuousStepLimit(track, x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

