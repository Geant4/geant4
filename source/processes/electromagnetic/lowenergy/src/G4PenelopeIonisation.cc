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
// $Id: G4PenelopeIonisation.cc,v 1.2 2003-06-19 14:39:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
//
// File name:     G4PenelopeIonisation
//
// Author:        Luciano Pandola
// 
// Creation date: March 2003
//
// Modifications:
// 
// 25.03.03 L.Pandola First implementation 
// 03.06.03 L.Pandola Added continuous part
// --------------------------------------------------------------
// La tabella di sezioni d'urto si potrebbe calcolare analiticamente (v. Compton)
// calcolando S0 anziche' S1 nella CalculateStoppingPower
// Attualmente le CS sono lette dallo stesso file di LowEnergy

#include "G4PenelopeIonisation.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4LogLogInterpolation.hh"
#include "G4EMDataSet.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4ShellVacancy.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"

G4PenelopeIonisation::G4PenelopeIonisation(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
    crossSectionHandler(0),
    theMeanFreePath(0),
    energySpectrum(0),
    shellVacancy(0),
    kineticEnergy1(0.0),
    cosThetaPrimary(1.0),
    energySecondary(0.0),
    cosThetaSecondary(0.0),
    iOsc(-1)
{
  cutForPhotons = 250.0*eV;
  cutForElectrons = 250.0*eV;
  verboseLevel = 0;
  ionizationEnergy = new std::map<G4int,G4DataVector*>;
  resonanceEnergy  = new std::map<G4int,G4DataVector*>;
  occupationNumber = new std::map<G4int,G4DataVector*>;
  shellFlag = new std::map<G4int,G4DataVector*>;
  ReadData(); //Read data from file
 
}


G4PenelopeIonisation::~G4PenelopeIonisation()
{
  delete crossSectionHandler;
  delete energySpectrum;
  delete theMeanFreePath;
  delete shellVacancy; 
  for (G4int Z=1;Z<100;Z++)
    {
      if (ionizationEnergy->count(Z)) delete (ionizationEnergy->find(Z)->second);
      if (resonanceEnergy->count(Z)) delete (resonanceEnergy->find(Z)->second);
      if (occupationNumber->count(Z)) delete (resonanceEnergy->find(Z)->second);
      if (shellFlag->count(Z)) delete (shellFlag->find(Z)->second);
    }
  delete ionizationEnergy;
  delete resonanceEnergy;
  delete occupationNumber;
  delete shellFlag;
}


void G4PenelopeIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType) 
{
  if(verboseLevel > 0) {
    G4cout << "G4PenelopeIonisation::BuildPhysicsTable start" << G4endl;
      }

  cutForDelta.clear();

  // Create and fill IonisationParameters once
  if( energySpectrum != 0 ) delete energySpectrum;
  energySpectrum = new G4eIonisationSpectrum();

  if(verboseLevel > 0) {
    G4cout << "G4VEnergySpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4CrossSectionHandler once

  if ( crossSectionHandler != 0 ) delete crossSectionHandler;
  G4VDataSetAlgorithm* interpolation = new G4SemiLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new G4eIonisationCrossSectionHandler(energySpectrum,
							     interpolation,
							     lowKineticEnergy,
							     highKineticEnergy,
							     totBin);
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");

  if (verboseLevel > 0) {
    G4cout << GetProcessName()
           << " is created; Cross section data: "
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: "
           << G4endl;
    energySpectrum->PrintData();
  }

  // Build loss table for IonisationIV

  BuildLossTable(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "The loss table is built"
           << G4endl;
      }

  if (&aParticleType==G4Electron::Electron()) {

    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable;
    CounterOfElectronProcess++;
    PrintInfoDefinition();  

  } else {

    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable;
    CounterOfPositronProcess++;
  }

  // Build mean free path data using cut values

  if( theMeanFreePath ) delete theMeanFreePath;
  theMeanFreePath = crossSectionHandler->
                    BuildMeanFreePathForMaterials(&cutForDelta);

  if(verboseLevel > 0) {
    G4cout << "The MeanFreePath table is built"
           << G4endl;
    if(verboseLevel > 1) theMeanFreePath->PrintData();
  }

  // Build common DEDX table for all ionisation processes
 
  BuildDEDXTable(aParticleType);

  if (verboseLevel > 0) {
    G4cout << "G4PenelopeIonisation::BuildPhysicsTable end"
           << G4endl;
  }
}


void G4PenelopeIonisation::BuildLossTable(
					  const G4ParticleDefinition&) //aParticleType)
{
  // Build table for energy loss due to soft brems
  // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss

  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  size_t   totBin = GetNbinEloss();
 
  //  create table

  if (theLossTable) { 
      theLossTable->clearAndDestroy();
      delete theLossTable;
  }
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  theLossTable = new G4PhysicsTable(numOfCouples);

  if (shellVacancy != 0) delete shellVacancy;
  shellVacancy = new G4ShellVacancy();
  G4DataVector* ksi = 0;
  G4DataVector* energy = 0;
  size_t binForFluo = totBin/10;

  G4PhysicsLogVector* bVector = new G4PhysicsLogVector(lowKineticEnergy,
		                		       highKineticEnergy,
						       binForFluo);
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  
  // Clean up the vector of cuts

  cutForDelta.clear();

  // Loop for materials

  for (size_t m=0; m<numOfCouples; m++) {

    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
		        				 highKineticEnergy,
							 totBin);

    // get material parameters needed for the energy loss calculation
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
    const G4Material* material= couple->GetMaterial();

    // the cut cannot be below lowest limit
    G4double tCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[m];
    tCut = std::min(tCut,highKineticEnergy);
    cutForDelta.push_back(tCut);

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector =
                    material->GetAtomicNumDensityVector();
    const G4double electronVolumeDensity = 
      material->GetTotNbOfElectPerVolume();  //electron density
    if(verboseLevel > 0) {
      G4cout << "Energy loss for material # " << m
             << " tCut(keV)= " << tCut/keV
             << G4endl;
      }

    // now comes the loop for the kinetic energy values
    for (size_t i = 0; i<totBin; i++) {

      G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
      G4double ionloss = 0.;

      // loop for elements in the material
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {

        G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
	ionloss   += 
	  CalculateContinuous(lowEdgeEnergy,tCut,Z,electronVolumeDensity) * 
	  theAtomicNumDensityVector[iel];

	if(verboseLevel > 1) {
	  G4cout << "Z= " << Z
		 << " E(keV)= " << lowEdgeEnergy/keV
		 << " loss= " << ionloss
		 << " rho= " << theAtomicNumDensityVector[iel]
		 << G4endl;
	}
      }
      //   if (material->GetName() == "Aluminum"){
      //       G4cout << "Material: " << material->GetName() << " Energia: " << lowEdgeEnergy/keV << " keV,"
      //       	       << "Stopping: " << ionloss/(keV/cm) << " kev/cm" << G4endl;
      //       }
      aVector->PutValue(i,ionloss);
    }
    theLossTable->insert(aVector);

    // fill data for fluorescence

    G4VDataSetAlgorithm* interp = new G4LogLogInterpolation();
    G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
    for (size_t iel=0; iel<NumberOfElements; iel++ ) {

      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      energy = new G4DataVector();
      ksi    = new G4DataVector();

      for (size_t j = 0; j<binForFluo; j++) {

        G4double lowEdgeEnergy = bVector->GetLowEdgeEnergy(j);
        G4double cross   = 0.;
        G4double eAverage= 0.;
	G4int nShells = transitionManager->NumberOfShells(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut,
                                                             lowEdgeEnergy, n);
          G4double pro = energySpectrum->Probability(Z, 0.0, tCut,
                                                             lowEdgeEnergy, n);
          G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
          eAverage   += e * cs * theAtomicNumDensityVector[iel];
          cross      += cs * pro * theAtomicNumDensityVector[iel];
          if(verboseLevel > 1) {
            G4cout << "Z= " << Z
                   << " shell= " << n
                   << " E(keV)= " << lowEdgeEnergy/keV
                   << " Eav(keV)= " << e/keV
                   << " pro= " << pro
                   << " cs= " << cs
                   << G4endl;
          }
	}

        G4double coeff = 0.0;
        if(eAverage > 0.) {
          coeff = cross/eAverage;
          eAverage /= cross;
	}

        if(verboseLevel > 1) {
            G4cout << "Ksi Coefficient for Z= " << Z
                   << " E(keV)= " << lowEdgeEnergy/keV
                   << " Eav(keV)= " << eAverage/keV
                   << " coeff= " << coeff
                   << G4endl;
        }

        energy->push_back(lowEdgeEnergy);
        ksi->push_back(coeff);
      }
      interp = new G4LogLogInterpolation();
      G4VEMDataSet* set = new G4EMDataSet(Z,energy,ksi,interp,1.,1.);
      xsis->AddComponent(set);
    }
    if(verboseLevel) xsis->PrintData();
    shellVacancy->AddXsiTable(xsis);
  }
  delete bVector;
}


G4VParticleChange* G4PenelopeIonisation::PostStepDoIt(const G4Track& track,
					               const G4Step&  step)
{
  aParticleChange.Initialize(track);

  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  const G4DynamicParticle* incidentElectron = track.GetDynamicParticle();
  const G4Material* material = couple->GetMaterial();
  const G4double electronVolumeDensity = 
    material->GetTotNbOfElectPerVolume();  //electron density
  G4double kineticEnergy0 = incidentElectron->GetKineticEnergy();
  G4ParticleMomentum electronDirection0 = incidentElectron->GetMomentumDirection();

  //Inizialisation of variables
  kineticEnergy1=kineticEnergy0;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;

  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy0);
  G4int    index  = couple->GetIndex();
  G4double tCut   = cutForDelta[index];

  CalculateDiscrete(kineticEnergy0,tCut,Z,electronVolumeDensity);
  //the method CalculateDiscrete() set the private variables:
  // kineticEnergy1 = energy of the primary electron after the interaction
  // cosThetaPrimary = cos(theta) of the primary after the interaction
  // energySecondary = energy of the secondary electron
  // cosThetaSecondary = cos(theta) of the secondary
 
  if(energySecondary == 0.0)
    {    
      return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
    }

  //Update the primary particle
  G4double sint = sqrt(1. - cosThetaPrimary*cosThetaPrimary);
  G4double phi  = twopi * G4UniformRand();
  G4double dirx = sint * cos(phi);
  G4double diry = sint * sin(phi);
  G4double dirz = cosThetaPrimary;

  G4ThreeVector electronDirection1(dirx,diry,dirz);
  electronDirection1.rotateUz(electronDirection0);
  aParticleChange.SetMomentumChange(electronDirection1) ;

  if (kineticEnergy1 > 0.)
    {
      aParticleChange.SetEnergyChange(kineticEnergy1) ;
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
    }

  //Generate the delta day
  G4int iosc2 = 0;
  G4double ionEnergy = 0.0;
  if (iOsc > 0) {
    ionEnergy=(*(ionizationEnergy->find(Z)->second))[iOsc];
    iosc2 = (ionizationEnergy->find(Z)->second->size()) - iOsc; //Sono in ordine inverso
  }

  // Verificare la scelta del livello da ionizzare

  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4double bindingEnergy = 0.0;
  G4int shellId = 0;
  if (iOsc > 0){
    const G4AtomicShell* shell = transitionManager->Shell(Z,iosc2);
    bindingEnergy = shell->BindingEnergy();
    shellId = shell->ShellId();
  }
  G4cout << "Oscillatore: " << iOsc << G4endl;
  G4cout << "Energie " << ionEnergy/keV << " " << bindingEnergy/keV << G4endl;
  ionEnergy = std::max(bindingEnergy,ionEnergy); //protection against energy non-conservation 
  G4double eKineticEnergy = energySecondary;

  size_t nTotPhotons=0;
  G4int nPhotons=0;

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t indx = couple->GetIndex();
  G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[indx];
  cutg = std::min(cutForPhotons,cutg);

  G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[indx];
  cute = std::min(cutForPhotons,cute);
  
  std::vector<G4DynamicParticle*>* photonVector=0;
  G4DynamicParticle* aPhoton;
  G4AtomicDeexcitation deexcitationManager;

  if (Z>5 && (ionEnergy > cutg || ionEnergy > cute))
    {
      photonVector = deexcitationManager.GenerateParticles(Z,shellId);
      nTotPhotons = photonVector->size();
      for (size_t k=0;k<nTotPhotons;k++){
	aPhoton = (*photonVector)[k];
	if (aPhoton)
	  {
	    G4double itsCut = cutg;
	    if (aPhoton->GetDefinition() == G4Electron::Electron()) itsCut = cute;
	    G4double itsEnergy = aPhoton->GetKineticEnergy();
	    if (itsEnergy > itsCut && itsEnergy <= ionEnergy)
	      {
		nPhotons++;
		ionEnergy -= itsEnergy;
	      }
	    else
	      {
		delete aPhoton;
		(*photonVector)[k]=0;
	      }
	  }
      }
    }
  G4double energyDeposit=ionEnergy; //il deposito locale e' quello che rimane
  G4int nbOfSecondaries=nPhotons;

  // Generate the delta ray 
  G4double sin2 = sqrt(1. - cosThetaSecondary*cosThetaSecondary);
  G4double phi2  = twopi * G4UniformRand();
  G4DynamicParticle* electron = 0;
  
  G4double xEl = sin2 * cos(phi2); 
  G4double yEl = sin2 * sin(phi2);
  G4double zEl = cosThetaSecondary;
  G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
  
  electron = new G4DynamicParticle (G4Electron::Electron(),
				    eDirection,eKineticEnergy) ;
  nbOfSecondaries++;

  aParticleChange.SetNumberOfSecondaries(nbOfSecondaries);
  if (electron) aParticleChange.AddSecondary(electron);
  for (size_t ll=0;ll<nTotPhotons;ll++)
    {
      aPhoton = (*photonVector)[ll];
      if (aPhoton) aParticleChange.AddSecondary(aPhoton);
    }
  delete photonVector;
  if (energyDeposit < 0)
    {
      G4cout << "WARNING-" 
	     << "G4PenelopeIonisaition::PostStepDoIt - Negative energy deposit"
	     << G4endl;
      energyDeposit=0;
    }
  aParticleChange.SetLocalEnergyDeposit(energyDeposit);
  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void G4PenelopeIonisation::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database.";
  comments += "\n      Delta energy sampled from a parametrised formula.";
  comments += "\n      Implementation of the continuous dE/dx part.";
  comments += "\n      At present it can be used for electrons ";
  comments += "in the energy range [250eV,100GeV].";
  comments += "\n      The process must work with G4PenelopeBremsstrahlung.";

  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}

G4bool G4PenelopeIonisation::IsApplicable(const G4ParticleDefinition& particle)
{
   return ( (&particle == G4Electron::Electron()) );
}

std::vector<G4DynamicParticle*>*
G4PenelopeIonisation::DeexciteAtom(const G4MaterialCutsCouple* couple,
			                  G4double incidentEnergy,
			                  G4double eLoss)
{
  // create vector of secondary particles
  const G4Material* material = couple->GetMaterial();

  std::vector<G4DynamicParticle*>* partVector =
                                 new std::vector<G4DynamicParticle*>;

  if(eLoss > cutForPhotons && eLoss > cutForElectrons) {

    const G4AtomicTransitionManager* transitionManager =
                               G4AtomicTransitionManager::Instance();

    size_t nElements = material->GetNumberOfElements();
    const G4ElementVector* theElementVector = material->GetElementVector();

    std::vector<G4DynamicParticle*>* secVector = 0;
    G4DynamicParticle* aSecondary = 0;
    G4ParticleDefinition* type = 0;
    G4double e;
    G4ThreeVector position;
    G4int shell, shellId;

    // sample secondaries

    G4double eTot = 0.0;
    std::vector<G4int> n =
           shellVacancy->GenerateNumberOfIonisations(couple,
                                                     incidentEnergy,eLoss);
    for (size_t i=0; i<nElements; i++) {

      G4int Z = (G4int)((*theElementVector)[i]->GetZ());
      size_t nVacancies = n[i];

      G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

      if (nVacancies && Z > 5 && (maxE>cutForPhotons || maxE>cutForElectrons)) {

	for (size_t j=0; j<nVacancies; j++) {

	  shell = crossSectionHandler->SelectRandomShell(Z, incidentEnergy);
          shellId = transitionManager->Shell(Z, shell)->ShellId();
	  G4double maxEShell =
                     transitionManager->Shell(Z, shell)->BindingEnergy();

          if (maxEShell>cutForPhotons || maxEShell>cutForElectrons ) {

	    secVector = deexcitationManager.GenerateParticles(Z, shellId);

	    if (secVector != 0) {

	      for (size_t l = 0; l<secVector->size(); l++) {

	        aSecondary = (*secVector)[l];
	        if (aSecondary != 0) {

	          e = aSecondary->GetKineticEnergy();
	          type = aSecondary->GetDefinition();
	          if ( eTot + e <= eLoss &&
	             (type == G4Gamma::Gamma() && e>cutForPhotons ) ||
	             (type == G4Electron::Electron() && e>cutForElectrons)) {

			  eTot += e;
                          partVector->push_back(aSecondary);

		  } else {

                           delete aSecondary;

	          }
	        }
	      }
              delete secVector;
	    }
	  }
	}
      }
    }
  }
  return partVector;
}

G4double G4PenelopeIonisation::GetMeanFreePath(const G4Track& track,
					       G4double, // previousStepSize
					       G4ForceCondition* cond)
{
   *cond = NotForced;
   G4int index = (track.GetMaterialCutsCouple())->GetIndex();
   const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath;
}

void G4PenelopeIonisation::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
}

void G4PenelopeIonisation::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForElectrons = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
}

void G4PenelopeIonisation::ActivateAuger(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}


void G4PenelopeIonisation::CalculateDiscrete(G4double ene,G4double cutoff,
					     G4int Z,G4double electronVolumeDensity)
{
  kineticEnergy1=ene;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;
  iOsc=-1;
  //constants
  G4double rb=ene+2.0*electron_mass_c2;
  G4double gamma = 1.0+ene/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double amol = pow((gamma-1.0)/gamma,2);
  G4double cps = ene*rb;
  G4double cp = sqrt(cps);
  
  G4double delta = CalculateDeltaFermi(ene,Z,electronVolumeDensity);
  G4double distantTransvCS0 = std::max(log(gamma2)-beta2-delta,0.0);

  G4double rl,rl1;
  G4DataVector* qm = new G4DataVector();
  G4DataVector* cumulHardCS = new G4DataVector();
  G4DataVector* typeOfInteraction = new G4DataVector();
  G4DataVector* nbOfLevel = new G4DataVector();
 
  if (cutoff > ene) return; //delta rays are not generated

  //Hard close collisions with outer shells
  G4double wmaxc = 0.5*ene;
  G4double closeCS0 = 0.0;
  G4double closeCS = 0.0;
  if (cutoff>0.1*eV) 
    {
      rl=cutoff/ene;
      rl1=1.0-rl;
      if (rl < 0.5)
	closeCS0 = (amol*(0.5-rl)+(1.0/rl)-(1.0/rl1)+(1.0-amol)*log(rl/rl1))/ene;
    }

  // Cross sections for the different oscillators

  // totalHardCS contains the cumulative hard interaction cross section for the different
  // excitable levels and the different interaction channels (close, distant, etc.),
  // i.e.
  // cumulHardCS[0] = 0.0
  // cumulHardCS[1] = 1st excitable level (distant longitudinal only)
  // cumulHardCS[2] = 1st excitable level (distant longitudinal + transverse)
  // cumulHardCS[3] = 1st excitable level (distant longitudinal + transverse + close)
  // cumulHardCS[4] = 1st excitable level (all channels) + 2nd excitable level (distant long only)
  // etc.
  // This is used for sampling the atomic level which is ionised and the channel of the
  // interaction.
  //
  // For each index iFill of the cumulHardCS vector,
  // nbOfLevel[iFill] contains the current excitable atomic level and
  // typeOfInteraction[iFill] contains the current interaction channel, with the legenda:
  //   1 = distant longitudinal interaction
  //   2 = distant transverse interaction
  //   3 = close collision
  //   4 = close collision with outer shells (in this case nbOfLevel < 0 --> no binding energy)


  G4int nOscil = ionizationEnergy->find(Z)->second->size();
  G4double totalHardCS = 0.0;
  G4double involvedElectrons = 0.0;
  for (G4int i=0;i<nOscil;i++){
    G4double wi = (*(resonanceEnergy->find(Z)->second))[i];
    G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
    //Distant excitations
    if (wi>cutoff && wi<ene)
      {
	if (wi>(1e-6*ene)){
	  G4double cpp=sqrt((ene-wi)*(ene-wi+2.0*electron_mass_c2));
	  qm->push_back(sqrt(pow(cp-cpp,2)+electron_mass_c2*electron_mass_c2)-electron_mass_c2);
	}
	else
	  {
	    qm->push_back(pow(wi,2)/(beta2+2.0*electron_mass_c2));
	  }
	//verificare che quando arriva qui il vettore ha SEMPRE l'i-esimo elemento
	if ((*qm)[i] < wi)
	  {
	    
	    G4double distantLongitCS =  occupNb*log(wi*((*qm)[i]+2.0*electron_mass_c2)/
					 ((*qm)[i]*(wi+2.0*electron_mass_c2)))/wi;
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(1.0); //distant longitudinal
	    nbOfLevel->push_back((G4double) i); //only excitable level are counted 
	    totalHardCS += distantLongitCS;
	    
	    G4double distantTransvCS = occupNb*distantTransvCS0/wi;
	    
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(2.0); //distant tranverse
	    nbOfLevel->push_back((G4double) i);
	    totalHardCS += distantTransvCS;
	  }
      }
    else 
      {
	qm->push_back(wi);
      }
    //close collisions
    if(wi < wmaxc){
      if (wi < cutoff) {
	involvedElectrons += occupNb;
      }
      else
	{
	  rl=wi/ene;
	  rl1=1.0-rl;
	  closeCS = occupNb*(amol*(0.5-rl)+(1.0/rl)-(1.0/rl1)+(1.0-amol)*log(rl/rl1))/ene;
	  cumulHardCS->push_back(totalHardCS);
	  typeOfInteraction->push_back(3.0); //close
	  nbOfLevel->push_back((G4double) i);
	  totalHardCS += closeCS;
	}
    }
  } // loop on the levels
  
  cumulHardCS->push_back(totalHardCS);
  typeOfInteraction->push_back(4.0); //close interaction with outer shells
  nbOfLevel->push_back(-1.0);
  totalHardCS += involvedElectrons*closeCS0;
  cumulHardCS->push_back(totalHardCS); //this is the final value of the totalHardCS

  if (totalHardCS < 1e-30) {
    kineticEnergy1=ene;
    cosThetaPrimary=1.0;
    energySecondary=0.0;
    cosThetaSecondary=0.0;
    iOsc=-1;
    return;
  }


  //Selection of the active oscillator on the basis of the cumulative cross sections
  G4double TST = totalHardCS*G4UniformRand();
  G4int is=0;
  G4int js= nbOfLevel->size();
  do{
    G4int it=(is+js)/2;
    if (TST > (*cumulHardCS)[it]) is=it;
    if (TST <= (*cumulHardCS)[it]) js=it;
  }while((js-is) > 1);

  G4double UII=0.0;
  G4double rkc=cutoff/ene;
  G4double dde;
  G4int kks;

  G4double sampledInteraction = (*typeOfInteraction)[is];
  iOsc = (G4int) (*nbOfLevel)[is];

  //Generates the final state according to the sampled level and 
  //interaction channel
  
  if (sampledInteraction == 1.0)  //Hard distant longitudinal collisions
    {
      G4cout << "Hard distant longitudinal collision" << G4endl;
      dde= (*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=ene-dde;
      G4double qs=(*qm)[iOsc]/(1.0+((*qm)[iOsc]/(2.0*electron_mass_c2)));
      G4double q=qs/(pow((qs/dde)*(1.0+(0.5*dde/electron_mass_c2)),G4UniformRand())-(0.5*qs/electron_mass_c2));
      G4double qtrev = q*(q+2.0*electron_mass_c2);
      G4double cpps = kineticEnergy1*(kineticEnergy1+2.0*electron_mass_c2);
      cosThetaPrimary = (cpps+cps-qtrev)/(2.0*cp*sqrt(cpps));
      if (cosThetaPrimary>1.0) cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      if (kks>4) 
	{
	  energySecondary=dde;
	}
      else
	{
	  energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
	}
      cosThetaSecondary = 0.5*(dde*(ene+rb-dde)+qtrev)/sqrt(cps*qtrev);
      if (cosThetaSecondary>1.0) cosThetaSecondary=1.0;
    }

  else if (sampledInteraction == 2.0)  //Hard distant transverse collisions
    {
      G4cout << "Hard distant transverse collision" << G4endl;
      dde=(*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=ene-dde;
      cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      if (kks>4)
	{
	  energySecondary=dde;
	}
      else
	{
	  energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
	}
      cosThetaSecondary = 1.0;
    }

  else if (sampledInteraction == 3.0 || sampledInteraction == 4.0) //Close interaction
    {
      if (sampledInteraction == 4.0) //interaction with inner shells
	{
	  UII=0.0;
	  rkc = cutoff/ene;
	  iOsc = -1;
	}
      else
	{
	  kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
	  if (kks > 4) {
	    UII=0.0;
	  }
	  else
	    {
	      UII = (*(ionizationEnergy->find(Z)->second))[iOsc];
	    }
	  rkc = (*(resonanceEnergy->find(Z)->second))[iOsc]/ene;
	}
    G4double A = 0.5*amol;
    G4double arkc = A*0.5*rkc;
    G4double phi,rk2,rk,rkf;
    do{
      G4double fb = (1.0+arkc)*G4UniformRand();
      if (fb<1.0)
	{
	  rk=rkc/(1.0-fb*(1.0-(rkc*2.0)));
	}
      else{
	rk = rkc+(fb-1.0)*(0.5-rkc)/arkc;
      }
      rk2 = rk*rk;
      rkf = rk/(1.0-rk);
      phi = 1.0+pow(rkf,2)-rkf+amol*(rk2+rkf);
    }while ((G4UniformRand()*(1.0+A*rk2)) > phi);
    //Energy and scattering angle (primary electron);
    kineticEnergy1 = ene*(1.0-rk);
    cosThetaPrimary = sqrt(kineticEnergy1*rb/(ene*(rb-(rk*ene))));
    //Energy and scattering angle of the delta ray
    energySecondary = ene-kineticEnergy1-UII;
    cosThetaSecondary = sqrt(rk*ene*rb/(ene*(rk*ene+2.0*electron_mass_c2)));
    }

  else

    {
      G4String excep = "G4PenelopeIonisation - Error in the calculation of the final state";
      G4Exception(excep);
    }

  delete qm;
  delete cumulHardCS;
  delete typeOfInteraction;
  delete nbOfLevel;

  return;
}

void G4PenelopeIonisation::ReadData()
{
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeIonisation - G4LEDATA environment variable not set!";
      G4Exception(excep);
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/ion-pen.dat";
  std::ifstream file(pathFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (!(lsdp->is_open()))
    {
      G4String excep = "G4PenelopeIonisation - data file " + pathFile + " not found!";
      G4Exception(excep);
    }

  G4int k1,test,test1,k2,k3;
  G4double a1,a2,a3,a4;
  G4int Z=1,nLevels=0;
  G4DataVector* x1;
  G4DataVector* x2;
  G4DataVector* x3;
  G4DataVector* x4;

  do{
    x1 = new G4DataVector;
    x2 = new G4DataVector;
    x3 = new G4DataVector;
    x4 = new G4DataVector;
    file >> Z >> nLevels;
    for (G4int h=0;h<nLevels;h++){
      //index,occup number,ion energy,res energy,fj0,kz,shell flag
      file >> k1 >> a1 >> a2 >> a3 >> a4 >> k2 >> k3;
      x1->push_back(a1); 
      x2->push_back(a2);
      x3->push_back(a3);
      x4->push_back((G4double) k3);
    }
    occupationNumber->insert(make_pair(Z,x1));
    ionizationEnergy->insert(make_pair(Z,x2));
    resonanceEnergy->insert(make_pair(Z,x3));
    shellFlag->insert(make_pair(Z,x4));
    file >> test >> test1; //-1 -1 close the data for each Z
    if (test > 0) {
      G4String excep = "G4PenelopeIonisation - data file corrupted!";
      G4Exception(excep);
    }
  }while (test != -2); //the very last Z is closed with -2 instead of -1
};


G4double G4PenelopeIonisation::CalculateDeltaFermi(G4double ene,G4int Z,
						   G4double electronVolumeDensity)
{
  G4double plasmaEnergyCoefficient = 1.377e-39; //(e*hbar)^2/(epsilon0*electron_mass)
  G4double plasmaEnergySquared = plasmaEnergyCoefficient*(electronVolumeDensity*m3);
  // sqrt(plasmaEnergySquared) is the plasma energy of the solid (MeV)
  G4double gam = 1.0+ene/electron_mass_c2;
  G4double gam2=gam*gam;
  G4double delta = 0.0;

  //Density effect
  G4double TST = ((G4double) Z)/(gam2*plasmaEnergySquared); 
  //C'e' una piccola differenza rispetto a Penelope, che viene da plasmaEnergySquared
  //(probabilmente le densita' elettroniche sono leggermente diverse!)
  G4double wl2 = 0.0;
  G4double fdel=0.0;
  G4double wr=0;
  size_t nbOsc = resonanceEnergy->find(Z)->second->size();
  for(size_t i=0;i<nbOsc;i++)
    {
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
      wr = (*(resonanceEnergy->find(Z)->second))[i];
      fdel += occupNb/(wr*wr+wl2);
    }
  if (fdel < TST) return delta;
  wl2 = pow((*(resonanceEnergy->find(Z)->second))[nbOsc-1],2);
  do{
    wl2=wl2*2.0;
    fdel = 0.0;
    for (i=0;i<nbOsc;i++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
      wr = (*(resonanceEnergy->find(Z)->second))[i];
      fdel += occupNb/(wr*wr+wl2);
    }
  }while (fdel > TST);
  G4double wl2l=0.0;
  G4double wl2u = wl2;
  G4double control = 0.0;
  do{
    wl2=0.5*(wl2l+wl2u);
    fdel = 0.0;
    for (i=0;i<nbOsc;i++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
      wr = (*(resonanceEnergy->find(Z)->second))[i];
      fdel += occupNb/(wr*wr+wl2);
    }
    if (fdel > TST)
      {
	wl2l = wl2;
      }
    else
      {
	wl2u = wl2;
      }
    control = wl2u-wl2l-wl2*1e-12; 
  }while(control>0);

  //Density correction effect
   for (i=0;i<nbOsc;i++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
      wr = (*(resonanceEnergy->find(Z)->second))[i];
      delta += occupNb*log(1.0+wl2/(wr*wr));
    }
   delta = (delta/((G4double) Z))-wl2/(gam2*plasmaEnergySquared); 
   return delta;
}

G4double G4PenelopeIonisation::CalculateContinuous(G4double ene,G4double cutoff,
					     G4int Z,G4double electronVolumeDensity)
{
  //Constants
  G4double gamma = 1.0+ene/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double constant = pi*pow(classic_electr_radius,2)*2.0*electron_mass_c2/beta2;
 
  
  G4double delta = CalculateDeltaFermi(ene,Z,electronVolumeDensity);
  G4int nbOsc = (G4int) resonanceEnergy->find(Z)->second->size();
  G4double S1 = 0.0;
  G4double stoppingPower = 0.0;
  for (G4int i=0;i<nbOsc;i++){
    G4double resEnergy = (*(resonanceEnergy->find(Z)->second))[i];
    S1 = CalculateStoppingPower(ene,resEnergy,delta,cutoff);
    //in the version of Penelope I run, cutoff=1 keV: if so, the results are identical
    G4double occupNb = (*(occupationNumber->find(Z)->second))[i];
    stoppingPower += occupNb*constant*S1;
  }
  
  return stoppingPower;
}

G4double G4PenelopeIonisation::CalculateStoppingPower(G4double ene,G4double resEne,
					     G4double delta,G4double cutoff)
{
  //Calculate constants
  G4double gamma = 1.0+ene/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double cps = ene*(ene+2.0*electron_mass_c2);
  G4double amol = pow((gamma-1.0)/gamma,2);
  G4double sPower = 0.0;
  if (ene < resEne) return sPower;

  //Distant interactions
  G4double cp1s = (ene-resEne)*(ene-resEne+2.0*electron_mass_c2);
  G4double cp1 = sqrt(cp1s);
  G4double cp = sqrt(cps);
  G4double sdLong=0.0, sdTrans = 0.0, sdDist=0.0;

  //Distant longitudinal interactions
  G4double qm = 0.0;

  if (resEne > ene*(1e-6))
    {
      qm = sqrt(pow(cp-cp1,2)+pow(electron_mass_c2,2))-electron_mass_c2;
    }
  else
    {
      qm = resEne*resEne/(beta2+2.0*electron_mass_c2);
      qm = qm*(1.0-0.5*qm/electron_mass_c2);
    }

  if (qm < resEne)
    {
      sdLong = log(resEne*(qm+2.0*electron_mass_c2)/(qm*(resEne+2.0*electron_mass_c2)));
    }
  else
    {
      sdLong = 0.0;
    }
  
  if (sdLong > 0) {
    sdTrans = std::max(log(gamma2)-beta2-delta,0.0);
    sdDist = sdTrans + sdLong;
    if (cutoff > resEne) sPower = sdDist;
  }


  // Close collisions (Moeller's cross section)
  G4double wl = std::max(cutoff,resEne);
  G4double wu = 0.5*ene;
 
  if (wl < (wu-1*eV)) wu=wl;
  wl = resEne;
  if (wl > (wu-1*eV)) return sPower;
  sPower += log(wu/wl)+(ene/(ene-wu))-(ene/(ene-wl))
    + (2.0 - amol)*log((ene-wu)/(ene-wl))
    + amol*(pow(wu,2)-pow(wl,2))/(2.0*pow(ene,2));

  return sPower;
}
