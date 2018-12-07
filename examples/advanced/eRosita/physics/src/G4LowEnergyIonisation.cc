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
// --------------------------------------------------------------
//
// File name:     G4LowEnergyIonisation
//
// Author:        Alessandra Forti, Vladimir Ivanchenko
// 
// Creation date: March 1999
//
// Modifications:
// - 11.04.2000 VL
//   Changing use of float and G4float casts to G4double casts 
//   because of problems with optimisation (bug ?)
//   10.04.2000 VL
// - Correcting Fluorescence transition probabilities in order to take into account 
//   non-radiative transitions. No Auger electron simulated yet: energy is locally deposited.
//   10.04.2000 VL
// - Correction of incident electron final momentum direction
//   07.04.2000 VL+LU
// - First implementation of continuous energy loss
//   22.03.2000 VL
// - 1 bug corrected in SelectRandomAtom method (units)
//   17.02.2000 Veronique Lefebure
// - 5 bugs corrected: 
//   *in Fluorescence, 2 bugs affecting 
//   . localEnergyDeposition and
//   . number of emitted photons that was then always 1 less
//   *in EnergySampling method: 
//   . expon = Parms[13]+1; (instead of uncorrect -1)
//   . rejection /= Parms[6];(instead of uncorrect Parms[7])
//   . Parms[6] is apparently corrupted in the data file (often = 0)  
//     -->Compute normalisation into local variable rejectionMax
//     and use rejectionMax  in stead of Parms[6]
//
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Added EnergySampling method A. Forti
// Modified PostStepDoIt to insert sampling with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// 20.09.00 V.Ivanchenko update fluctuations 
// 24.04.01 V.Ivanchenko remove RogueWave 
// 22.05.01 V.Ivanchenko update calculation of delta-ray kinematic + 
//                       clean up the code 
// 02.08.01 V.Ivanchenko fix energy conservation for small steps 
// 18.08.01 V.Ivanchenko fix energy conservation for pathalogical delta-energy
// 01.10.01 E. Guardincerri Replaced fluorescence generation in PostStepDoIt
//                          according to design iteration
// 04.10.01 MGP             Minor clean-up in the fluo section, removal of
//                          compilation warnings and extra protection to
//                          prevent from accessing a null pointer        
// 29.09.01 V.Ivanchenko    revision based on design iteration
// 10.10.01 MGP             Revision to improve code quality and 
//                          consistency with design
// 18.10.01 V.Ivanchenko    Add fluorescence AlongStepDoIt
// 18.10.01 MGP             Revision to improve code quality and
//                          consistency with design
// 19.10.01 V.Ivanchenko    update according to new design, V.Ivanchenko
// 26.10.01 V.Ivanchenko    clean up deexcitation
// 28.10.01 V.Ivanchenko    update printout
// 29.11.01 V.Ivanchenko    New parametrisation introduced
// 25.03.02 V.Ivanchneko    Fix in fluorescence
// 28.03.02 V.Ivanchenko    Add flag of fluorescence
// 28.05.02 V.Ivanchenko    Remove flag fStopAndKill
// 31.05.02 V.Ivanchenko    Add path of Fluo + Auger cuts to
//                          AtomicDeexcitation
// 03.06.02 MGP             Restore fStopAndKill
// 19.06.02 VI              Additional printout
// 30.07.02 VI              Fix in restricted energy loss
// 20.09.02 VI              Remove ActivateFlurescence from SetCut...
// 21.01.03 VI              Cut per region
// 12.02.03 VI              Change signature for Deexcitation
// 12.04.03 V.Ivanchenko    Cut per region for fluo AlongStep
// 31.08.04 V.Ivanchenko    Add density correction
//
// --------------------------------------------------------------

#include "G4LowEnergyIonisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RDeIonisationSpectrum.hh"
#include "G4RDeIonisationCrossSectionHandler.hh"
#include "G4RDAtomicTransitionManager.hh"
#include "G4RDAtomicShell.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDSemiLogInterpolation.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4RDShellVacancy.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"

G4LowEnergyIonisation::G4LowEnergyIonisation(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0),
  shellVacancy(0)
{
  cutForPhotons = 250.0*eV;
  cutForElectrons = 250.0*eV;
  verboseLevel = 0;
}


G4LowEnergyIonisation::~G4LowEnergyIonisation()
{
  delete crossSectionHandler;
  delete energySpectrum;
  delete theMeanFreePath;
  delete shellVacancy;
}


void G4LowEnergyIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable start"
           << G4endl;
      }

  cutForDelta.clear();

  // Create and fill IonisationParameters once
  if( energySpectrum != 0 ) delete energySpectrum;
  energySpectrum = new G4RDeIonisationSpectrum();

  if(verboseLevel > 0) {
    G4cout << "G4RDVEnergySpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4RDCrossSectionHandler once

  if ( crossSectionHandler != 0 ) delete crossSectionHandler;
  G4RDVDataSetAlgorithm* interpolation = new G4RDSemiLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new G4RDeIonisationCrossSectionHandler(energySpectrum,
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
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable end"
           << G4endl;
  }
}


void G4LowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& )
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
  shellVacancy = new G4RDShellVacancy();
  G4DataVector* ksi = 0;
  G4DataVector* energy = 0;
  size_t binForFluo = totBin/10;

  G4PhysicsLogVector* bVector = new G4PhysicsLogVector(lowKineticEnergy,
		                		       highKineticEnergy,
						       binForFluo);
  const G4RDAtomicTransitionManager* transitionManager = G4RDAtomicTransitionManager::Instance();
  
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
    if(tCut > highKineticEnergy) tCut = highKineticEnergy;
    cutForDelta.push_back(tCut);
    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector =
                    material->GetAtomicNumDensityVector();
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

	G4int nShells = transitionManager->NumberOfShells(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut,
                                                             lowEdgeEnergy, n);
          G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
          ionloss   += e * cs * theAtomicNumDensityVector[iel];

          if(verboseLevel > 1 || (Z == 14 && lowEdgeEnergy>1. && lowEdgeEnergy<0.)) {
            G4cout << "Z= " << Z
                   << " shell= " << n
                   << " E(keV)= " << lowEdgeEnergy/keV
                   << " Eav(keV)= " << e/keV
                   << " cs= " << cs
	           << " loss= " << ionloss
	           << " rho= " << theAtomicNumDensityVector[iel]
                   << G4endl;
          }
        }
        G4double esp = energySpectrum->Excitation(Z, lowEdgeEnergy);
        ionloss   += esp * theAtomicNumDensityVector[iel];

      }
      if(verboseLevel > 1 || (m == 0 && lowEdgeEnergy>=1. && lowEdgeEnergy<=0.)) {
            G4cout << "Sum: "
                   << " E(keV)= " << lowEdgeEnergy/keV
	           << " loss(MeV/mm)= " << ionloss*mm/MeV
                   << G4endl;
      }
      aVector->PutValue(i,ionloss);
    }
    theLossTable->insert(aVector);

    // fill data for fluorescence

    G4RDVDataSetAlgorithm* interp = new G4RDLogLogInterpolation();
    G4RDVEMDataSet* xsis = new G4RDCompositeEMDataSet(interp, 1., 1.);
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
      interp = new G4RDLogLogInterpolation();
      G4RDVEMDataSet* set = new G4RDEMDataSet(Z,energy,ksi,interp,1.,1.);
      xsis->AddComponent(set);
    }
    if(verboseLevel) xsis->PrintData();
    shellVacancy->AddXsiTable(xsis);
  }
  delete bVector;
}


G4VParticleChange* G4LowEnergyIonisation::PostStepDoIt(const G4Track& track,
					               const G4Step&  step)
{
  // Delta electron production mechanism on base of the model
  // J. Stepanek " A program to determine the radiation spectra due
  // to a single atomic subshell ionisation by a particle or due to
  // deexcitation or decay of radionuclides",
  // Comp. Phys. Comm. 1206 pp 1-19 (1997)

  aParticleChange.Initialize(track);

  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  G4double kineticEnergy = track.GetKineticEnergy();

  // Select atom and shell

  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
  G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  const G4RDAtomicShell* atomicShell =
                (G4RDAtomicTransitionManager::Instance())->Shell(Z, shell);
  G4double bindingEnergy = atomicShell->BindingEnergy();
  G4int shellId = atomicShell->ShellId();

  // Sample delta energy

  G4int    index  = couple->GetIndex();
  G4double tCut   = cutForDelta[index];
  G4double tmax   = energySpectrum->MaxEnergyOfSecondaries(kineticEnergy);
  G4double tDelta = energySpectrum->SampleEnergy(Z, tCut, tmax,
                                                 kineticEnergy, shell);

  if(tDelta == 0.0)
    return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  // Transform to shell potential
  G4double deltaKinE = tDelta + 2.0*bindingEnergy;
  G4double primaryKinE = kineticEnergy + 2.0*bindingEnergy;

  // sampling of scattering angle neglecting atomic motion
  G4double deltaMom = std::sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
  G4double primaryMom = std::sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));

  G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
                            / (deltaMom * primaryMom);

  if (cost > 1.) cost = 1.;
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = twopi * G4UniformRand();
  G4double dirx = sint * std::cos(phi);
  G4double diry = sint * std::sin(phi);
  G4double dirz = cost;

  // Rotate to incident electron direction
  G4ThreeVector primaryDirection = track.GetMomentumDirection();
  G4ThreeVector deltaDir(dirx,diry,dirz);
  deltaDir.rotateUz(primaryDirection);
  dirx = deltaDir.x();
  diry = deltaDir.y();
  dirz = deltaDir.z();


  // Take into account atomic motion del is relative momentum of the motion
  // kinetic energy of the motion == bindingEnergy in V.Ivanchenko model

  cost = 2.0*G4UniformRand() - 1.0;
  sint = std::sqrt(1. - cost*cost);
  phi  = twopi * G4UniformRand();
  G4double del = std::sqrt(bindingEnergy *(bindingEnergy + 2.0*electron_mass_c2))
               / deltaMom;
  dirx += del* sint * std::cos(phi);
  diry += del* sint * std::sin(phi);
  dirz += del* cost;

  // Find out new primary electron direction
  G4double finalPx = primaryMom*primaryDirection.x() - deltaMom*dirx;
  G4double finalPy = primaryMom*primaryDirection.y() - deltaMom*diry;
  G4double finalPz = primaryMom*primaryDirection.z() - deltaMom*dirz;

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
  theDeltaRay->SetKineticEnergy(tDelta);
  G4double norm = 1.0/std::sqrt(dirx*dirx + diry*diry + dirz*dirz);
  dirx *= norm;
  diry *= norm;
  dirz *= norm;
  theDeltaRay->SetMomentumDirection(dirx, diry, dirz);
  theDeltaRay->SetDefinition(G4Electron::Electron());

  G4double theEnergyDeposit = bindingEnergy;

  // fill ParticleChange
  // changed energy and momentum of the actual particle

  G4double finalKinEnergy = kineticEnergy - tDelta - theEnergyDeposit;
  if(finalKinEnergy < 0.0) {
    theEnergyDeposit += finalKinEnergy;
    finalKinEnergy    = 0.0;
    aParticleChange.ProposeTrackStatus(fStopAndKill);

  } else {

    G4double norm = 1.0/std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
    finalPx *= norm;
    finalPy *= norm;
    finalPz *= norm;
    aParticleChange.ProposeMomentumDirection(finalPx, finalPy, finalPz);
  }

  aParticleChange.ProposeEnergy(finalKinEnergy);

  // Generation of Fluorescence and Auger
  size_t nSecondaries = 0;
  size_t totalNumber  = 1;
  std::vector<G4DynamicParticle*>* secondaryVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;

  // Fluorescence data start from element 6

  if (Fluorescence() && Z > 5 && (bindingEnergy >= cutForPhotons
            ||  bindingEnergy >= cutForElectrons)) {

    secondaryVector = deexcitationManager.GenerateParticles(Z, shellId);

    if (secondaryVector != 0) {

      nSecondaries = secondaryVector->size();
      for (size_t i = 0; i<nSecondaries; i++) {

        aSecondary = (*secondaryVector)[i];
        if (aSecondary) {

          G4double e = aSecondary->GetKineticEnergy();
          type = aSecondary->GetDefinition();
          if (e < theEnergyDeposit &&
                ((type == G4Gamma::Gamma() && e > cutForPhotons ) ||
                 (type == G4Electron::Electron() && e > cutForElectrons ))) {

             theEnergyDeposit -= e;
             totalNumber++;

	  } else {

             delete aSecondary;
             (*secondaryVector)[i] = 0;
	  }
	}
      }
    }
  }

  // Save delta-electrons

  aParticleChange.SetNumberOfSecondaries(totalNumber);
  aParticleChange.AddSecondary(theDeltaRay);

  // Save Fluorescence and Auger

  if (secondaryVector) {

    for (size_t l = 0; l < nSecondaries; l++) {

      aSecondary = (*secondaryVector)[l];

      if(aSecondary) {

        aParticleChange.AddSecondary(aSecondary);
      }
    }
    delete secondaryVector;
  }

  if(theEnergyDeposit < 0.) {
    G4cout << "G4LowEnergyIonisation: Negative energy deposit: "
           << theEnergyDeposit/eV << " eV" << G4endl;
    theEnergyDeposit = 0.0;
  }
  aParticleChange.ProposeLocalEnergyDeposit(theEnergyDeposit);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void G4LowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database.";
  comments += "\n      Gamma energy sampled from a parametrised formula.";
  comments += "\n      Implementation of the continuous dE/dx part.";
  comments += "\n      At present it can be used for electrons ";
  comments += "in the energy range [250eV,100GeV].";
  comments += "\n      The process must work with G4LowEnergyBremsstrahlung.";

  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}

G4bool G4LowEnergyIonisation::IsApplicable(const G4ParticleDefinition& particle)
{
   return ( (&particle == G4Electron::Electron()) );
}

std::vector<G4DynamicParticle*>*
G4LowEnergyIonisation::DeexciteAtom(const G4MaterialCutsCouple* couple,
			                  G4double incidentEnergy,
			                  G4double eLoss)
{
  // create vector of secondary particles
  const G4Material* material = couple->GetMaterial();

  std::vector<G4DynamicParticle*>* partVector =
                                 new std::vector<G4DynamicParticle*>;

  if(eLoss > cutForPhotons && eLoss > cutForElectrons) {

    const G4RDAtomicTransitionManager* transitionManager =
                               G4RDAtomicTransitionManager::Instance();

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
	             ((type == G4Gamma::Gamma() && e>cutForPhotons ) ||
	             (type == G4Electron::Electron() && e>cutForElectrons))) {

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

G4double G4LowEnergyIonisation::GetMeanFreePath(const G4Track& track,
						G4double , // previousStepSize
						G4ForceCondition* cond)
{
   *cond = NotForced;
   G4int index = (track.GetMaterialCutsCouple())->GetIndex();
   const G4RDVEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath;
}

void G4LowEnergyIonisation::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
}

void G4LowEnergyIonisation::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForElectrons = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
}

void G4LowEnergyIonisation::ActivateAuger(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}

