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
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisationMA physics process -------
//                by Vladimir Ivanchenko, 14 July 1999 
//                was made on the base of G4hIonisation class
//                developed by Laszlo Urban
// ************************************************************
// It is the extention of the ionisation process for the slow 
// charged hadrons.
// ************************************************************
// 28 July   1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli added ICRU parametrisations for protons
// 20 August 1999 G.Mancinelli added ICRU tables for alpha
// 31 August 1999 V.Ivanchenko update and cleen up
// 30 Sept.  1999 V.Ivanchenko minor upgrade
// 12 Dec.   1999 S. Chauvie added Barkas correction
// 19 Jan.   2000 V.Ivanchenko minor changing in Barkas corrections
// 02 April  2000 S. Chauvie linearization of Barkas effect
// 03 April  2000 V.Ivanchenko Nuclear Stopping power for antiprotons
// 23 May    2000 MG Pia  Clean up for QAO model 
// 24 May    2000 MG Pia  Code properly indented to improve legibility
// 17 July   2000 V.Ivanchenko Bug in scaling AlongStepDoIt method
// 25 July   2000 V.Ivanchenko New design iteration
// 17 August 2000 V.Ivanchenko Add ion fluctuation models
// 18 August 2000 V.Ivanchenko Bug fixed in GetConstrain
// 22 August 2000 V.Ivanchenko Insert paramStepLimit and
//                reorganise access to Barkas and Bloch terms  
// 04 Sept.  2000 V.Ivanchenko rename fluctuations
// 05 Sept.  2000 V.Ivanchenko clean up
// 03 Oct.   2000 V.Ivanchenko CodeWizard clean up
// 03 Nov.   2000 V.Ivanchenko MinKineticEnergy=LowestKineticEnergy=10eV
// 05 Nov.   2000 MG Pia - Removed const cast previously introduced to get
//                the code compiled (const G4Material* now introduced in 
//                electromagnetic/utils utils-V02-00-03 tag)
//                (this is going back and forth, to cope with Michel's
//                utils tag not being accepted yet by system testing)
// 21 Nov.  2000 V.Ivanchenko Fix a problem in fluctuations
// 23 Nov.  2000 V.Ivanchenko Ion type fluctuations only for charge>0
// 10 May   2001 V.Ivanchenko Clean up againist Linux compilation with -Wall
// 23 May   2001 V.Ivanchenko Minor fix in PostStepDoIt
// 07 June  2001 V.Ivanchenko Clean up AntiProtonDEDX + add print out
// 18 June  2001 V.Ivanchenko Cleanup print out
// 18 Oct.  2001 V.Ivanchenko Add fluorescence
// 30 Oct.  2001 V.Ivanchenko Add minGammaEnergy and minElectronEnergy
// 07 Dec   2001 V.Ivanchenko Add SetFluorescence method
// 15 Feb   2002 V.Ivanchenko Fix problem of Generic Ions
// 25 Mar   2002 V.Ivanchenko Fix problem of fluorescence below threshold
// 28 Mar   2002 V.Ivanchenko Set fluorescence off by default
// 09 Apr   2002 V.Ivanchenko Fix table problem of GenericIons
// 28 May   2002 V.Ivanchenko Remove flag fStopAndKill
// 31 May   2002 V.Ivanchenko Add path of Fluo + Auger cuts to
//                            AtomicDeexcitation
// 03 Jun   2002 MGP          Restore fStopAndKill
// 10 Jun   2002 V.Ivanchenko Restore fStopButAlive
// 12 Jun   2002 V.Ivanchenko Fix in fluctuations - if tmax<2*Ipot Gaussian
//                            fluctuations enables
// 20 Sept  2002 V.Ivanchenko Clean up energy ranges for models
// 07 Oct   2002 V.Ivanchenko Clean up initialisation of fluorescence
// 28 Oct   2002 V.Ivanchenko Optimal binning for dE/dx
// 10 Dec   2002 V.Ivanchenko antiProtonLowEnergy -> 25 keV, QEG model below
// 21 Jan   2003 V.Ivanchenko Cut per region
// 10 Mar   2003 V.Ivanchenko Use SubTypes for ions
// 12 Apr   2003 V.Ivanchenko Cut per region for fluo AlongStep
// 18 Apr   2003 V.Ivanchenko finalRange redefinition
// 26 Apr   2003 V.Ivanchenko fix for stepLimit
// 06 May   2004 V.Ivanchenko Migrate to G4VEnergyLossProcess interface

// -----------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hLowEnergyIonisationMA.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Poisson.hh"
#include "G4UnitsTable.hh"
#include "G4hParametrisedLossModel.hh"
#include "G4LowEnergyBraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4Material.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4ShellVacancy.hh"
#include "G4hShellCrossSection.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4BohrFluctuations.hh"
#include "G4IonFluctuations.hh"
#include "G4Region.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisationMA::G4hLowEnergyIonisationMA(const G4String& processName)
  : G4VEnergyLossProcess(processName),
    theTable("ICRU_R49p"),
    flucModel(0),
    shellVacancy(0),
    shellCS(0),
    theMaterial(0),
    currentParticle(0),
    theParticle(0),
    theBaseParticle(0),
    fluobins(20),
    theBarkas(true),
    theFluo(false),
    fluoIsInitialised(false)
{
  highEnergy         = 2.*MeV;
  minGammaEnergy     = 25.*keV;
  minElectronEnergy  = 25.*keV;
  verboseLevel       = 0;
  SetDEDXBinning(360);
  SetLambdaBinning(360);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*GeV);
  shellCS = new G4hShellCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisationMA::~G4hLowEnergyIonisationMA()
{
  if(shellVacancy) delete shellVacancy;
  if(shellCS) delete shellCS;
  G4int length = zFluoDataVector.size();
  if(length) {
    for(G4int i=0; i<length; i++) {
      delete zFluoDataVector[i];
    }
    zFluoDataVector.clear();
  }
  regionsWithFluo.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4hLowEnergyIonisationMA::DefineBaseParticle(
                      const G4ParticleDefinition* p)
{
  const G4ParticleDefinition* bp = 0;
  if(p) {
    theParticle = p;
    G4double q = p->GetPDGCharge();
    if(q > 0.0 && p != G4Proton::Proton())              bp = G4Proton::Proton();
    else if(q < 0.0 && p != G4AntiProton::AntiProton()) bp = G4AntiProton::AntiProton();
    if(p->GetPDGCharge() < 0.0) theTable = "QAO";
  }
  theBaseParticle = bp;
  if(!isInitialised) InitialiseProcess();
  return bp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::InitialiseProcess()
{
  if(theParticle->GetPDGCharge()/eplus > 1.5 ||
     theParticle->GetParticleName() == "GenericIon")
                  flucModel = new G4IonFluctuations();
  else            flucModel = new G4BohrFluctuations();

  SetSecondaryParticle(G4Electron::Electron());

  G4hParametrisedLossModel* param = new G4hParametrisedLossModel(theTable);
  G4double lowEnergy = MinKinEnergy();
  G4double maxEnergy = MaxKinEnergy();

  G4LowEnergyBraggModel* em = new G4LowEnergyBraggModel();
  em->SetLowEnergyLimit(lowEnergy);
  em->SetHighEnergyLimit(highEnergy);
  em->SetParameterization(param);
  AddEmModel(1, em, flucModel);

  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(highEnergy);
  em1->SetHighEnergyLimit(maxEnergy);
  AddEmModel(2, em1, flucModel);

  SetIntegral(true);
  chargeLowLimit = 0.1;
  energyLowLimit = 250.*MeV;
  SetLinearLossLimit(0.15);
  SetStepLimits(0.1, 0.1*mm);

  isInitialised = true;

  if(verboseLevel > 0) {
    G4cout << "G4hLowEnergyIonisationMA::InitiliseProcess done for "
           << theParticle->GetParticleName();
    if(theBaseParticle)	G4cout << "  base particle " << theBaseParticle->GetParticleName();
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::BuildDataForFluorescence()
{
  fluoIsInitialised = true;
  if(verboseLevel > 1) {
    G4cout << "G4hLowEnergyIonisationMA::BuildDataForFluorescence for "
           << theParticle->GetParticleName() << " is started" << G4endl;
  }

  G4double lowEnergy = MinKinEnergy();
  G4double maxEnergy = MaxKinEnergy();

  // fill data for fluorescence

  deexcitationManager.SetCutForSecondaryPhotons(minGammaEnergy);
  deexcitationManager.SetCutForAugerElectrons(minElectronEnergy);

  G4double mass = theParticle->GetPDGMass();
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if (shellVacancy != 0) delete shellVacancy;
  shellVacancy = new G4ShellVacancy();
  G4DataVector* ksi = 0;
  G4DataVector* ksi1 = 0;
  G4DataVector* energy = 0;
  G4DataVector* energy1 = 0;
  G4int length = zFluoDataVector.size();
  if(length > 0) {
    for(G4int i=0; i<length; i++) {
      G4VEMDataSet* x = zFluoDataVector[i];
      delete x;
    }
    zFluoDataVector.clear();
  }

  G4PhysicsLogVector* bVector = new G4PhysicsLogVector(lowEnergy,
		                		       maxEnergy,
						       fluobins);
  const G4AtomicTransitionManager* transitionManager =
                             G4AtomicTransitionManager::Instance();

  G4double bindingEnergy;

  //  loop for materials
  for (size_t j=0; j<numOfCouples; j++) {

    // get material parameters needed for the energy loss calculation
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector =
                    material->GetAtomicNumDensityVector();
    G4VDataSetAlgorithm* interp = new G4SemiLogInterpolation();
    G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
    G4VDataSetAlgorithm* interp1 = new G4SemiLogInterpolation();
    G4VEMDataSet* xsis1 = new G4CompositeEMDataSet(interp1, 1., 1.);

    G4double tCut = (*theCoupleTable->GetEnergyCutsVector(1))[j];
    G4double elDensity = 1.;

    for (size_t iel=0; iel<NumberOfElements; iel++ ) {

      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      G4int nShells = transitionManager->NumberOfShells(Z);
      energy = new G4DataVector();
      ksi    = new G4DataVector();
      energy1= new G4DataVector();
      ksi1   = new G4DataVector();
      //if(NumberOfElements > 1)
      elDensity = theAtomicNumDensityVector[iel]/((G4double)nShells);

      for (size_t j = 0; j<fluobins; j++) {

        G4double tkin  = bVector->GetLowEdgeEnergy(j);
        G4double gamma = tkin/mass + 1.;
        G4double beta2 = 1.0 - 1.0/(gamma*gamma);
        G4double r     = electron_mass_c2/mass;
        G4double tmax  = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*r + r*r);
        G4double cross   = 0.;
        G4double cross1  = 0.;
        G4double eAverage= 0.;
        G4double tmin = std::min(tCut,tmax);
        G4double rel;

        for (G4int n=0; n<nShells; n++) {

          bindingEnergy = transitionManager->Shell(Z, n)->BindingEnergy();
          if (tmin > bindingEnergy) {
            rel = log(tmin/bindingEnergy);
            eAverage   += rel - beta2*(tmin - bindingEnergy)/tmax;
            cross      += 1.0/bindingEnergy - 1.0/tmin - beta2*rel/tmax;
	  }
          if (tmax > tmin) {
	    cross1     += 1.0/tmin - 1.0/tmax - beta2*log(tmax/tmin)/tmax;
	  }
	}

        cross1 *= elDensity;
        energy1->push_back(tkin);
        ksi1->push_back(cross1);

        if(eAverage > 0.) cross /= eAverage;
        else              cross  = 0.;

        energy->push_back(tkin);
        ksi->push_back(cross);
      }
      G4VDataSetAlgorithm* algo = interp->Clone();
      G4VEMDataSet* set = new G4EMDataSet(Z,energy,ksi,algo,1.,1.);
      xsis->AddComponent(set);
      G4VDataSetAlgorithm* algo1 = interp1->Clone();
      G4VEMDataSet* set1 = new G4EMDataSet(Z,energy1,ksi1,algo1,1.,1.);
      xsis1->AddComponent(set1);
    }
    if(verboseLevel > 1) {
      G4cout << "### Shell inverse cross sections for "
             << material->GetName() << G4endl;
      xsis->PrintData();
      G4cout << "### Atom cross sections for "
             << material->GetName() << G4endl;
      xsis1->PrintData();
    }
    shellVacancy->AddXsiTable(xsis);
    zFluoDataVector.push_back(xsis1);
  }
  delete bVector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4Track*>* G4hLowEnergyIonisationMA::SecondariesAlongStep(
                             const G4Step& stepData,
                                   G4double& tmax,
                                   G4double& eloss,
                                   G4double& kineticEnergy)

{
  std::vector<G4DynamicParticle*>* newpart = 0;
  std::vector<G4Track*>* tracks = 0;
  G4DynamicParticle* part = 0;

  const G4Track* track = stepData.GetTrack();
  const G4MaterialCutsCouple* couple = track->GetMaterialCutsCouple();
  G4double hMass = track->GetDynamicParticle()->GetMass();

  if(theFluo) newpart = DeexciteAtom(couple, kineticEnergy, tmax, hMass, eloss);

  if(newpart) {

    size_t nSecondaries = newpart->size();
    fParticleChange.SetNumberOfSecondaries(nSecondaries);
    G4Track* newtrack = 0;
    const G4StepPoint* preStep = stepData.GetPreStepPoint();
    const G4StepPoint* postStep = stepData.GetPostStepPoint();
    G4ThreeVector r = preStep->GetPosition();
    G4ThreeVector deltaR = postStep->GetPosition();
    deltaR -= r;
    G4double t = preStep->GetGlobalTime();
    G4double deltaT = postStep->GetGlobalTime();
    deltaT -= t;
    G4double time, q, e;
    G4ThreeVector position;

    for(size_t i=0; i<nSecondaries; i++) {

      part = (*newpart)[i];
      if(part) {

        e = part->GetKineticEnergy();
        q = G4UniformRand();
        time = deltaT*q + t;
        position  = deltaR*q;
        position += r;
        newtrack = new G4Track(part, time, position);
        tracks->push_back(newtrack);

      }
    }
    delete newpart;
  }
  return tracks;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::SecondariesPostStep(
                                   G4VEmModel* model,
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dynParticle,
                                   G4double& tcut,
                                   G4double& kinEnergy)

{

  G4DynamicParticle* delta = model->SampleSecondary(couple, dynParticle, tcut, kinEnergy);

  if (delta) {

    G4ThreeVector finalP = dynParticle->GetMomentum();
    G4double deltaKinEnergy = delta->GetKineticEnergy();
    G4double mass = dynParticle->GetMass();
    kinEnergy -= deltaKinEnergy;
    finalP -= delta->GetMomentum();
    finalP  = finalP.unit();
    fParticleChange.SetProposedMomentumDirection(finalP);

    // Generation of Fluorescence and Auger
    size_t nSecondaries = 0;
    size_t totalNumber  = 1;
    std::vector<G4DynamicParticle*>* secondaryVector = 0;
    G4DynamicParticle* aSecondary = 0;
    G4ParticleDefinition* type = 0;

    if (theFluo && (kinEnergy > minGammaEnergy || kinEnergy > minElectronEnergy)) {

      if(!fluoIsInitialised) BuildDataForFluorescence();
      // Select atom and shell
      G4int Z = SelectRandomAtom(couple, kinEnergy);

      if(Z > 5) {
        G4int shell = shellCS->SelectRandomShell(Z, kinEnergy,
                                             mass,deltaKinEnergy);
        const G4AtomicShell* atomicShell =
                (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
        G4double bindingEnergy = atomicShell->BindingEnergy();

        if(verboseLevel > 1) {
          G4cout << "PostStep Z= " << Z << " shell= " << shell
		 << " bindingE(keV)= " << bindingEnergy/keV
		 << " deltaE(keV)= " << deltaKinEnergy/keV
		 << " kinE(keV)= " << kinEnergy/keV
		 << G4endl;
	}

	// Fluorescence data start from element 6

	if (kinEnergy >= bindingEnergy &&
           (bindingEnergy >= minGammaEnergy || bindingEnergy >= minElectronEnergy) ) {

	  G4int shellId = atomicShell->ShellId();
	  secondaryVector = deexcitationManager.GenerateParticles(Z, shellId);

	  if (secondaryVector != 0) {

	    nSecondaries = secondaryVector->size();
	    for (size_t i = 0; i<nSecondaries; i++) {

	      aSecondary = (*secondaryVector)[i];
	      if (aSecondary) {

		G4double e = aSecondary->GetKineticEnergy();
		type = aSecondary->GetDefinition();
		if (e < kinEnergy &&
		    ((type == G4Gamma::Gamma() && e > minGammaEnergy ) ||
		     (type == G4Electron::Electron() && e > minElectronEnergy ))) {

		  kinEnergy -= e;
		  totalNumber++;

		} else {

		  delete aSecondary;
		  (*secondaryVector)[i] = 0;
		}
	      }
	    }
	  }
	}
      }
    }
    fParticleChange.SetNumberOfSecondaries(totalNumber);
    fParticleChange.AddSecondary(delta);

    // Save Fluorescence and Auger

    if (secondaryVector) {

      for (size_t l = 0; l < nSecondaries; l++) {

        aSecondary = (*secondaryVector)[l];
        if(aSecondary) fParticleChange.AddSecondary(aSecondary);
      }
    }
    delete secondaryVector;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4hLowEnergyIonisationMA::DeexciteAtom(
                                       const G4MaterialCutsCouple* couple,
	 			             G4double incidentEnergy,
				             G4double tmax,
				             G4double hMass,
				             G4double eLoss)
{

  if (verboseLevel > 1) {
	G4cout << "DeexciteAtom: cutForPhotons(keV)= " << minGammaEnergy/keV
               << "  cutForElectrons(keV)= " << minElectronEnergy/keV
               << "  eLoss(MeV)= " << eLoss
               << G4endl;
  }

  if(!fluoIsInitialised) BuildDataForFluorescence();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  const G4Material* material = couple->GetMaterial();
  G4int index    = couple->GetIndex();
  minGammaEnergy    = (*theCoupleTable->GetEnergyCutsVector(0))[index];
  minElectronEnergy = (*theCoupleTable->GetEnergyCutsVector(1))[index];
  G4double tcut = std::min(minGammaEnergy, minElectronEnergy);

  if(eLoss < tcut) return 0;

  const G4AtomicTransitionManager* transitionManager =
                             G4AtomicTransitionManager::Instance();

  size_t nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4bool stop = true;

  for (size_t j=0; j<nElements; j++) {

    G4int Z = (G4int)((*theElementVector)[j]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (Z > 5 && maxE < eLoss ) {
      stop = false;
      break;
    }
  }

  if(stop) return 0;

  G4double tau   = incidentEnergy/hMass;
  G4double tau1  = tau + 1.0;
  G4double beta2 = tau*(tau + 2.0)/(tau1*tau1);

  // create vector of tracks of secondary particles
  std::vector<G4DynamicParticle*>* partVector =
         new std::vector<G4DynamicParticle*>;
  std::vector<G4DynamicParticle*>* secVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;
  G4double e, tkin, grej;
  G4ThreeVector position;
  G4int shell, shellId;

  // sample secondaries

  G4double etot = 0.0;
  std::vector<G4int> n = shellVacancy->GenerateNumberOfIonisations(couple,
                                         incidentEnergy, eLoss);

  for (size_t i=0; i<nElements; i++) {

    size_t nVacancies = n[i];
    G4int Z = (G4int)((*theElementVector)[i]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (nVacancies && Z>5 && (maxE>minGammaEnergy || maxE>minElectronEnergy )) {
      for(size_t j=0; j<nVacancies; j++) {

        // sampling follows
        do {
	  tkin = tcut/(1.0 + (tcut/maxE - 1.0)*G4UniformRand());
          grej = 1.0 - beta2 * tkin/tmax;

        } while( G4UniformRand() > grej );

        shell = shellCS->SelectRandomShell(Z,incidentEnergy,hMass,tkin);

        shellId = transitionManager->Shell(Z, shell)->ShellId();
        G4double maxE = transitionManager->Shell(Z, shell)->BindingEnergy();

        if (maxE>minGammaEnergy || maxE>minElectronEnergy ) {
          secVector = deexcitationManager.GenerateParticles(Z, shellId);
	} else {
          secVector = 0;
	}

        if (secVector) {

          for (size_t l = 0; l<secVector->size(); l++) {

            aSecondary = (*secVector)[l];
            if(aSecondary) {

              e = aSecondary->GetKineticEnergy();
              type = aSecondary->GetDefinition();
              if ( etot + e <= eLoss &&
                   (type == G4Gamma::Gamma() && e > minGammaEnergy ) ||
                   (type == G4Electron::Electron() && e > minElectronEnergy)) {

                     etot += e;
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

  if(partVector->empty()) {
    delete partVector;
    return 0;
  }

  return partVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4hLowEnergyIonisationMA::SelectRandomAtom(const G4MaterialCutsCouple* couple,
                                                       G4double kineticEnergy) const
{
  const G4Material* material = couple->GetMaterial();
  G4int nElements = material->GetNumberOfElements();
  G4int Z = 0;

  if(nElements == 1) {
    Z = (G4int)(material->GetZ());
    return Z;
  }

  const G4ElementVector* theElementVector = material->GetElementVector();
  std::vector<G4double> p;
  G4int index = couple->GetIndex();

  G4double norm = 0.0;
  for (G4int j=0; j<nElements; j++) {

    const G4VEMDataSet* set = (zFluoDataVector[index])->GetComponent(j);
    G4double cross    = set->FindValue(kineticEnergy);

    p.push_back(cross);
    norm += cross;
  }

  if(norm == 0.0) return 0;

  G4double q = norm*G4UniformRand();

  for (G4int i=0; i<nElements; i++) {

    if(p[i] > q) {
       Z = (G4int)((*theElementVector)[i]->GetZ());
       break;
    }
    q -= p[i];
  }

  return Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisationMA::BarkasTerm(const G4Material* material,
  				                    G4double kineticEnergy) const
//Function to compute the Barkas term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
//
{
  static double FTable[47][2] = {
   { 0.02, 21.5},
   { 0.03, 20.0},
   { 0.04, 18.0},
   { 0.05, 15.6},
   { 0.06, 15.0},
   { 0.07, 14.0},
   { 0.08, 13.5},
   { 0.09, 13.},
   { 0.1,  12.2},
   { 0.2,  9.25},
   { 0.3,  7.0},
   { 0.4,  6.0},
   { 0.5,  4.5},
   { 0.6,  3.5},
   { 0.7,  3.0},
   { 0.8,  2.5},
   { 0.9,  2.0},
   { 1.0,  1.7},
   { 1.2,  1.2},
   { 1.3,  1.0},
   { 1.4,  0.86},
   { 1.5,  0.7},
   { 1.6,  0.61},
   { 1.7,  0.52},
   { 1.8,  0.5},
   { 1.9,  0.43},
   { 2.0,  0.42},
   { 2.1,  0.3},
   { 2.4,  0.2},
   { 3.0,  0.13},
   { 3.08, 0.1},
   { 3.1,  0.09},
   { 3.3,  0.08},
   { 3.5,  0.07},
   { 3.8,  0.06},
   { 4.0,  0.051},
   { 4.1,  0.04},
   { 4.8,  0.03},
   { 5.0,  0.024},
   { 5.1,  0.02},
   { 6.0,  0.013},
   { 6.5,  0.01},
   { 7.0,  0.009},
   { 7.1,  0.008},
   { 8.0,  0.006},
   { 9.0,  0.0032},
   { 10.0, 0.0025} };

  // Information on particle and material
  G4double kinE  = kineticEnergy ;
  if(0.5*MeV > kinE) kinE = 0.5*MeV ;
  G4double gamma = 1.0 + kinE / proton_mass_c2 ;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma) ;
  if(0.0 >= beta2) return 0.0;

  G4double BarkasTerm = 0.0;
  G4double AMaterial = 0.0;
  G4double ZMaterial = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int numberOfElements = material->GetNumberOfElements();

  for (G4int i = 0; i<numberOfElements; i++) {

    AMaterial = (*theElementVector)[i]->GetA()*mole/g;
    ZMaterial = (*theElementVector)[i]->GetZ();

    G4double X = 137.0 * 137.0 * beta2 / ZMaterial;

    // Variables to compute L_1
    G4double Eta0Chi = 0.8;
    G4double EtaChi = Eta0Chi * ( 1.0 + 6.02*pow( ZMaterial,-1.19 ) );
    G4double W = ( EtaChi * pow( ZMaterial,1.0/6.0 ) ) / sqrt(X);
    G4double FunctionOfW = FTable[46][1]*FTable[46][0]/W ;

    for(G4int j=0; j<47; j++) {

      if( W < FTable[j][0] ) {

        if(0 == j) {
     	  FunctionOfW = FTable[0][1] ;

        } else {
          FunctionOfW = (FTable[j][1] - FTable[j-1][1]) * (W - FTable[j-1][0])
	              / (FTable[j][0] - FTable[j-1][0])
	              +  FTable[j-1][1] ;
	}

        break;
      }

    }

    BarkasTerm += FunctionOfW /( sqrt(ZMaterial * X) * X);
  }

  BarkasTerm *= twopi_mc2_rcl2 * (material->GetElectronDensity()) / beta2 ;

  return BarkasTerm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisationMA::BlochTerm(const G4Material* material,
                                                 G4double kineticEnergy,
                                                 G4double cSquare) const
//Function to compute the Bloch term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
//
{
  G4double eloss = 0.0 ;
  G4double gamma = 1.0 + kineticEnergy / proton_mass_c2 ;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma) ;
  G4double y = cSquare / (137.0*137.0*beta2) ;

  if(y < 0.05) {
    eloss = 1.202 ;

  } else {
    eloss = 1.0 / (1.0 + y) ;
    G4double de = eloss ;

    for(G4int i=2; de>eloss*0.01; i++) {
      de = 1.0/( i * (i*i + y)) ;
      eloss += de ;
    }
  }
  eloss *= -1.0 * y * cSquare * twopi_mc2_rcl2 *
            (material->GetElectronDensity()) / beta2 ;

  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::ActivateFluorescence(
                               G4bool val, const G4Region*)
{
  theFluo = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::ActivateAugerElectronProduction(
                               G4bool val, const G4Region*)
{
  if(val) theFluo = val;
  deexcitationManager.ActivateAugerElectronProduction(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisationMA::EffectiveCharge(
                                          const G4ParticleDefinition* p,
                                          const G4Material* material,
			                        G4double kineticEnergy)
{
  G4double mass   = p->GetPDGMass();
  G4double charge = p->GetPDGCharge();
  G4double Zi     = charge/eplus;

  chargeCorrection = 1.0;

  // The aproximation of ion effective charge from:
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Fast ions or hadrons
  G4double reducedEnergy = kineticEnergy * proton_mass_c2/mass ;
  if( reducedEnergy > energyLowLimit || Zi < 1.5 ) return charge ;

  static G4double vFermi[92] = {
    1.0309,  0.15976, 0.59782, 1.0781,  1.0486,  1.0,     1.058,   0.93942, 0.74562, 0.3424,
    0.45259, 0.71074, 0.90519, 0.97411, 0.97184, 0.89852, 0.70827, 0.39816, 0.36552, 0.62712,
    0.81707, 0.9943,  1.1423,  1.2381,  1.1222,  0.92705, 1.0047,  1.2,     1.0661,  0.97411,
    0.84912, 0.95,    1.0903,  1.0429,  0.49715, 0.37755, 0.35211, 0.57801, 0.77773, 1.0207,
    1.029,   1.2542,  1.122,   1.1241,  1.0882,  1.2709,  1.2542,  0.90094, 0.74093, 0.86054,
    0.93155, 1.0047,  0.55379, 0.43289, 0.32636, 0.5131,  0.695,   0.72591, 0.71202, 0.67413,
    0.71418, 0.71453, 0.5911,  0.70263, 0.68049, 0.68203, 0.68121, 0.68532, 0.68715, 0.61884,
    0.71801, 0.83048, 1.1222,  1.2381,  1.045,   1.0733,  1.0953,  1.2381,  1.2879,  0.78654,
    0.66401, 0.84912, 0.88433, 0.80746, 0.43357, 0.41923, 0.43638, 0.51464, 0.73087, 0.81065,
    1.9578,  1.0257} ;

  static G4double lFactor[92] = {
    1.0,  1.0,  1.1,  1.06, 1.01, 1.03, 1.04, 0.99, 0.95, 0.9,
    0.82, 0.81, 0.83, 0.88, 1.0,  0.95, 0.97, 0.99, 0.98, 0.97,
    0.98, 0.97, 0.96, 0.93, 0.91, 0.9,  0.88, 0.9,  0.9,  0.9,
    0.9,  0.85, 0.9,  0.9,  0.91, 0.92, 0.9,  0.9,  0.9,  0.9,
    0.9,  0.88, 0.9,  0.88, 0.88, 0.9,  0.9,  0.88, 0.9,  0.9,
    0.9,  0.9,  0.96, 1.2,  0.9,  0.88, 0.88, 0.85, 0.9,  0.9,
    0.92, 0.95, 0.99, 1.03, 1.05, 1.07, 1.08, 1.1,  1.08, 1.08,
    1.08, 1.08, 1.09, 1.09, 1.1,  1.11, 1.12, 1.13, 1.14, 1.15,
    1.17, 1.2,  1.18, 1.17, 1.17, 1.16, 1.16, 1.16, 1.16, 1.16,
    1.16, 1.16} ;

  static G4double c[6] = {0.2865,  0.1266, -0.001429,
                          0.02402,-0.01135, 0.001475} ;

  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = material->GetNumberOfElements() ;

  //  loop for the elements in the material
  //  to find out average values Z, vF, lF
  G4double z = 0.0, vF = 0.0, lF = 0.0, norm = 0.0 ;

  if( 1 == NumberOfElements ) {
    z = material->GetZ() ;
    G4int iz = G4int(z) - 1 ;
    if(iz < 0) iz = 0 ;
    else if(iz > 91) iz = 91 ;
    vF   = vFermi[iz] ;
    lF   = lFactor[iz] ;

  } else {
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        const G4Element* element = (*theElementVector)[iel] ;
        G4double z2 = element->GetZ() ;
        const G4double weight = theAtomicNumDensityVector[iel] ;
        norm += weight ;
        z    += z2 * weight ;
        G4int iz = G4int(z2) - 1 ;
        if(iz < 0) iz = 0 ;
        else if(iz > 91) iz =91 ;
        vF   += vFermi[iz] * weight ;
        lF   += lFactor[iz] * weight ;
      }
    z  /= norm ;
    vF /= norm ;
    lF /= norm ;
  }

  G4double q;
  // Helium ion case
  if( Zi < 2.5 ) {

    // Normalise to He4 mass
    G4double e = log(std::max(1.0, kineticEnergy / (keV*4.0026) ) );
    G4double x = c[0] ;
    G4double y = 1.0 ;
    for (G4int i=1; i<6; i++) {
      y *= e ;
      x += y * c[i] ;
    }
    q = 7.6 -  e ;
    q = 1.0 + ( 0.007 + 0.00005 * z ) * exp( -q*q ) * sqrt(1.0 - exp(-x)) ;
    if( q < chargeLowLimit ) q = chargeLowLimit ;

    // Heavy ion case
  } else {
    // v1 is ion velocity in vF unit
    G4double v1 = sqrt( reducedEnergy / (25.0 * keV) )/ vF ;
    G4double y ;
    G4double z13 = pow(Zi, 0.3333) ;

    // Faster than Fermi velocity
    if ( v1 > 1.0 ) {
      y = vF * v1 * ( 1.0 + 0.2 / (v1*v1) ) / (z13*z13) ;

      // Slower than Fermi velocity
    } else {
      y = 0.6923 * vF * (1.0 + 2.0*v1*v1/3.0 + v1*v1*v1*v1/15.0) / (z13*z13) ;
    }

    G4double y3 = pow(y, 0.3) ;
    //    G4cout << "y= " << y << " y3= " << y3 << " v1= " << v1 << " vF= " << vF << G4endl; 
    q = 1.0 - exp( 0.803*y3 - 1.3167*y3*y3 - 0.38157*y - 0.008983*y*y ) ;

    if( q < chargeLowLimit ) q = chargeLowLimit ;

    G4double s = 7.6 -  log(std::max(1.0, reducedEnergy/keV)) ;
    s = 1.0 + ( 0.18 + 0.0015 * z ) * exp( -s*s )/ (Zi*Zi) ;

    // Screen length according to
    // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
    // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

    G4double lambda = 10.0 * vF * pow(1.0-q, 0.6667) / (z13 * (6.0 + q)) ;
    chargeCorrection = s * (1.0 + 0.5*(1.0/q - 1.0)*log(1.0 + lambda*lambda)/(vF*vF) );
  }
  //  G4cout << "G4ionIonisation: charge= " << charge << " q= " << q
  //       << " chargeCor= " << chargeCorrection << G4endl; 
  return charge*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisationMA::PrintInfoDefinition()
{
  G4VEnergyLossProcess::PrintInfoDefinition();
  G4cout << "      Bether-Bloch model for Escaled > " << highEnergy << " MeV, "
         << "parametrisation of Bragg peak below, "
         << "Integral mode " << IsIntegral()
         << G4endl;
  if(theBarkas){
    G4cout << "        Parametrization of the Barkas effect is switched on."
           << G4endl ;
  }

  G4bool printHead = true;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  // loop for materials
  for (size_t j=0 ; j < numOfCouples; j++) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();
    G4double cut  = (*theCoupleTable->GetEnergyCutsVector(1))[couple->GetIndex()];
    G4double eexc = material->GetIonisation()->GetMeanExcitationEnergy();

    if(eexc > cut) {
      if(printHead) {
        printHead = false ;

        G4cout << "           material       min.delta energy(keV) " << G4endl;
        G4cout << G4endl;
      }

      G4cout << std::setw(20) << material->GetName()
	     << std::setw(15) << eexc/keV << G4endl;
    }
  }
}
