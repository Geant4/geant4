//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                            *
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
// Geant4 Class file
//  
// Authors: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// Created 22 April 2010 from old G4UAtomicDeexcitation class 
//
// Modified:
// ---------
// 20 Oct 2011  Alf  modified to take into account ECPSSR form Form Factor
// 03 Nov 2011  Alf  Extended Empirical and Form Factor ionisation XS models
//                   out thei ranges with Analytical one.
// 07 Nov 2011  Alf  Restored original ioniation XS for alphas, 
//                   letting scaled ones for other ions.   
// 20 Mar 2012  LP   Register G4PenelopeIonisationCrossSection
//
// -------------------------------------------------------------------
//
// Class description:
// Implementation of atomic deexcitation 
//
// -------------------------------------------------------------------

#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4FluoTransition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"

#include "G4teoCrossSection.hh"
#include "G4empCrossSection.hh"
#include "G4PenelopeIonisationCrossSection.hh"
#include "G4LivermoreIonisationCrossSection.hh"
#include "G4EmCorrections.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4Material.hh"
#include "G4AtomicShells.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UAtomicDeexcitation::G4UAtomicDeexcitation():
  G4VAtomDeexcitation("UAtomDeexcitation"),
  minGammaEnergy(DBL_MAX), 
  minElectronEnergy(DBL_MAX),
  newShellId(-1)
{
  anaPIXEshellCS = nullptr;
  PIXEshellCS    = nullptr;
  ePIXEshellCS   = nullptr;
  emcorr = G4LossTableManager::Instance()->EmCorrections();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  transitionManager = G4AtomicTransitionManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UAtomicDeexcitation::~G4UAtomicDeexcitation()
{
  delete anaPIXEshellCS;
  delete PIXEshellCS;
  delete ePIXEshellCS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UAtomicDeexcitation::InitialiseForNewRun()
{
  if(!IsFluoActive()) { return; }
  transitionManager->Initialise();
  if(!IsPIXEActive()) { return; }

  if(!anaPIXEshellCS) {
    anaPIXEshellCS = new G4teoCrossSection("ECPSSR_Analytical");
  }
  G4cout << G4endl;
  G4cout << "### === G4UAtomicDeexcitation::InitialiseForNewRun()" << G4endl;

  G4EmParameters* param = G4EmParameters::Instance();
  G4String namePIXExsModel = param->PIXECrossSectionModel();
  G4String namePIXExsElectronModel = param->PIXEElectronCrossSectionModel();

  // Check if old cross section for p/ion should be deleted 
  if(PIXEshellCS && namePIXExsModel != PIXEshellCS->GetName()) 
    {
      delete PIXEshellCS;
      PIXEshellCS = nullptr;
    }

  // Instantiate new proton/ion cross section
  if(!PIXEshellCS) {
    if (namePIXExsModel == "ECPSSR_FormFactor")
      {
	PIXEshellCS = new G4teoCrossSection(namePIXExsModel);
      }
    else if(namePIXExsModel == "ECPSSR_ANSTO")
      {
	PIXEshellCS = new G4teoCrossSection(namePIXExsModel);
      }    
    else if(namePIXExsModel == "Empirical")
      {
	PIXEshellCS = new G4empCrossSection(namePIXExsModel);
      }
  }

  // Check if old cross section for e+- should be deleted 
  if(ePIXEshellCS && namePIXExsElectronModel != ePIXEshellCS->GetName()) 
    {
      delete ePIXEshellCS;
      ePIXEshellCS = nullptr;
    } 

  // Instantiate new e+- cross section
  if(nullptr == ePIXEshellCS) 
    {
      if(namePIXExsElectronModel == "Empirical")
	{
	  ePIXEshellCS = new G4empCrossSection("Empirical");
	}
      else if(namePIXExsElectronModel == "ECPSSR_Analytical") 
	{
	  ePIXEshellCS = new G4teoCrossSection("ECPSSR_Analytical");
	}
      else if (namePIXExsElectronModel == "Penelope")
	{
	  ePIXEshellCS = new G4PenelopeIonisationCrossSection();
	}
      else 
	{
	  ePIXEshellCS = new G4LivermoreIonisationCrossSection();
	}
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UAtomicDeexcitation::InitialiseForExtraAtom(G4int)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4AtomicShell* 
G4UAtomicDeexcitation::GetAtomicShell(G4int Z, G4AtomicShellEnumerator shell)
{
  return transitionManager->Shell(Z, (std::size_t)shell);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UAtomicDeexcitation::GenerateParticles(
		            std::vector<G4DynamicParticle*>* vectorOfParticles,
			    const G4AtomicShell* atomicShell, 
			    G4int Z,
			    G4double gammaCut,
			    G4double eCut)
{
  // Defined initial conditions
  G4int givenShellId = atomicShell->ShellId();
  minGammaEnergy = gammaCut;
  minElectronEnergy = eCut;

  // generation secondaries
  G4DynamicParticle* aParticle=0;
  G4int provShellId = 0;
  
  //ORIGINAL METHOD BY ALFONSO MANTERO
  if (!IsAugerCascadeActive())
    {
      //----------------------------  
      G4int counter = 0;
      
      // limits of the EPDL data
      if (Z>5 && Z<105) {

	// The aim of this loop is to generate more than one fluorecence photon 
	// from the same ionizing event 
	do
	  {
	    if (counter == 0) 
	      // First call to GenerateParticles(...):
	      // givenShellId is given by the process
	      {
		provShellId = SelectTypeOfTransition(Z, givenShellId);
		
		if (provShellId >0) 
		  {
		    aParticle =
		      GenerateFluorescence(Z, givenShellId, provShellId);
		  }
		else if (provShellId == -1)
		  {
		    aParticle = GenerateAuger(Z, givenShellId);
		  }
	      }
	    else 
	      // Following calls to GenerateParticles(...):
	      // newShellId is given by GenerateFluorescence(...)
	      {
		provShellId = SelectTypeOfTransition(Z,newShellId);
		if (provShellId >0)
		  {
		    aParticle = GenerateFluorescence(Z,newShellId,provShellId);
		  }
		else if ( provShellId == -1)
		  {
		    aParticle = GenerateAuger(Z, newShellId);		
		  }
	      }
	    ++counter;
	    if (aParticle != 0) 
	      {
		vectorOfParticles->push_back(aParticle);
	      }
	    else {provShellId = -2;}
	  }  
	while (provShellId > -2); 
      }
    } // Auger cascade is not active

  //END OF ORIGINAL METHOD BY ALFONSO MANTERO
  //----------------------

  // NEW METHOD
  // Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)
  if (IsAugerCascadeActive())
    {
      //----------------------
      vacancyArray.push_back(givenShellId);

      // let's check that 5<Z<100
      if (Z<6 || Z>104){
	return;
      }

      // as long as there is vacancy to be filled by either fluo or auger, stay in the loop.
      while(!vacancyArray.empty()){
	//  prepare to process the last element, and then delete it from the vector.
	givenShellId = vacancyArray[0];
	provShellId = SelectTypeOfTransition(Z,givenShellId);

	//G4cout<<"\n------ Atom Transition with Z: "<<Z<<"\tbetween current:"
	//		<<givenShellId<<" & target:"<<provShellId<<G4endl;
	if(provShellId>0){
	  aParticle = GenerateFluorescence(Z,givenShellId,provShellId);
	}
	else if(provShellId == -1){
	  aParticle = GenerateAuger(Z, givenShellId);
	}
	//  if a particle is created, put it in the vector of new particles
	if(aParticle!=0)
	  vectorOfParticles->push_back(aParticle);

	//  one vacancy has been processed. Erase it.
	vacancyArray.erase(vacancyArray.begin());
      }
      //----------------------
      //End of Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)

    } // Auger cascade is active
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4UAtomicDeexcitation::GetShellIonisationCrossSectionPerAtom(
		       const G4ParticleDefinition* pdef, 
		       G4int Z, 
		       G4AtomicShellEnumerator shellEnum,
		       G4double kineticEnergy,
		       const G4Material* mat)
{
  // we must put a control on the shell that are passed: 
  // some shells should not pass (line "0" or "2")

  // check atomic number
  G4double xsec = 0.0;
  if(Z > 93 || Z < 6 ) { return xsec; } //corrected by alf - Z<6 missing
  G4int idx = G4int(shellEnum);
  if(idx >= G4AtomicShells::GetNumberOfShells(Z)) { return xsec; }

  if(pdef == theElectron || pdef == thePositron) {
    xsec = ePIXEshellCS->CrossSection(Z,shellEnum,kineticEnergy,0.0,mat);
    return xsec;
  }

  G4double mass = pdef->GetPDGMass();
  G4double escaled = kineticEnergy;
  G4double q2 = 0.0;

  // scaling to protons for all particles excluding protons and alpha
  G4int pdg = pdef->GetPDGEncoding();
  if (pdg != 2212 && pdg != 1000020040)
    {
      mass = proton_mass_c2;
      escaled = kineticEnergy*mass/(pdef->GetPDGMass());

      if(mat) {
	q2 = emcorr->EffectiveChargeSquareRatio(pdef,mat,kineticEnergy);
      } else {
	G4double q = pdef->GetPDGCharge()/eplus;
	q2 = q*q;
      }
    }
  
  if(PIXEshellCS) {
    xsec = PIXEshellCS->CrossSection(Z,shellEnum,escaled,mass,mat);
  }
  if(xsec < 1e-100) {     
    xsec = anaPIXEshellCS->CrossSection(Z,shellEnum,escaled,mass,mat); 
  }

  if (q2)  {xsec *= q2;}

  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UAtomicDeexcitation::SetCutForSecondaryPhotons(G4double cut)
{
  minGammaEnergy = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UAtomicDeexcitation::SetCutForAugerElectrons(G4double cut)
{
  minElectronEnergy = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4UAtomicDeexcitation::ComputeShellIonisationCrossSectionPerAtom(
				const G4ParticleDefinition* p, 
				G4int Z, 
				G4AtomicShellEnumerator shell,
				G4double kinE,
				const G4Material* mat)
{
  return GetShellIonisationCrossSectionPerAtom(p,Z,shell,kinE,mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4UAtomicDeexcitation::SelectTypeOfTransition(G4int Z, G4int shellId)
{
  if (shellId <=0 ) {
    return 0;
  }
  
  G4int provShellId = -1;
  G4int shellNum = 0;
  G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);  
  
  const G4FluoTransition* refShell = 
    transitionManager->ReachableShell(Z,maxNumOfShells-1);

  // This loop gives shellNum the value of the index of shellId
  // in the vector storing the list of the shells reachable through
  // a radiative transition
  if ( shellId <= refShell->FinalShellId())
    {
      while (shellId != transitionManager->ReachableShell(Z,shellNum)->FinalShellId())
	{
	  if(shellNum ==maxNumOfShells-1)
	    {
	      break;
	    }
	  shellNum++;
	}
      G4int transProb = 0; //AM change 29/6/07 was 1
   
      G4double partialProb = G4UniformRand();      
      G4double partSum = 0;
      const G4FluoTransition* aShell = transitionManager->ReachableShell(Z,shellNum);
      G4int trSize =  (G4int)(aShell->TransitionProbabilities()).size();
    
      // Loop over the shells wich can provide an electron for a 
      // radiative transition towards shellId:
      // in every loop the partial sum of the first transProb shells
      // is calculated and compared with a random number [0,1].
      // If the partial sum is greater, the shell whose index is transProb
      // is chosen as the starting shell for a radiative transition
      // and its identity is returned
      // Else, terminateded the loop, -1 is returned
      while(transProb < trSize){
	partSum += aShell->TransitionProbability(transProb);

	if(partialProb <= partSum)
	  {
	    provShellId = aShell->OriginatingShellId(transProb);
	    break;
	  }
	++transProb;
      }
      // here provShellId is the right one or is -1.
      // if -1, the control is passed to the Auger generation part of the package 
    }
  else 
    {
      provShellId = -1;
    }
  return provShellId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* 
G4UAtomicDeexcitation::GenerateFluorescence(G4int Z, G4int shellId,
					    G4int provShellId )
{ 
  if (shellId <=0 )
    {
      return nullptr;
    }

  //isotropic angular distribution for the outcoming photon
  G4double newcosTh = 1.-2.*G4UniformRand();
  G4double newsinTh = std::sqrt((1.-newcosTh)*(1. + newcosTh));
  G4double newPhi = twopi*G4UniformRand();
  
  G4double xDir = newsinTh*std::sin(newPhi);
  G4double yDir = newsinTh*std::cos(newPhi);
  G4double zDir = newcosTh;
  
  G4ThreeVector newGammaDirection(xDir,yDir,zDir);
  
  G4int shellNum = 0;
  G4int maxNumOfShells = transitionManager->NumberOfReachableShells(Z);
  
  // find the index of the shell named shellId
  while (shellId != transitionManager->
	 ReachableShell(Z,shellNum)->FinalShellId())
    {
      if(shellNum == maxNumOfShells-1)
	{
	  break;
	}
      ++shellNum;
    }
  // number of shell from wich an electron can reach shellId
  G4int transitionSize = (G4int)transitionManager->
    ReachableShell(Z,shellNum)->OriginatingShellIds().size();
  
  G4int index = 0;
  
  // find the index of the shell named provShellId in the vector
  // storing the shells from which shellId can be reached 
  while (provShellId != transitionManager->
	 ReachableShell(Z,shellNum)->OriginatingShellId(index))
    {
      if(index ==  transitionSize-1)
	{
	  break;
	}
      ++index;
    }
  // energy of the gamma leaving provShellId for shellId
  G4double transitionEnergy = transitionManager->
    ReachableShell(Z,shellNum)->TransitionEnergy(index);
  
  if (transitionEnergy < minGammaEnergy) return nullptr;

  // This is the shell where the new vacancy is: it is the same
  // shell where the electron came from
  newShellId = transitionManager->
    ReachableShell(Z,shellNum)->OriginatingShellId(index);
    
  G4DynamicParticle* newPart = new G4DynamicParticle(G4Gamma::Gamma(), 
						     newGammaDirection,
						     transitionEnergy);

  //Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)
  if (IsAugerCascadeActive()) vacancyArray.push_back(newShellId);

  return newPart;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4UAtomicDeexcitation::GenerateAuger(G4int Z, G4int shellId)
{
  if(!IsAugerActive()) { 
    //    G4cout << "auger inactive!" << G4endl; //debug
    return nullptr; 
  }
  
  if (shellId <=0 ) {
    //G4Exception("G4UAtomicDeexcitation::GenerateAuger()","de0002",
    //		JustWarning, "Energy deposited locally");
    return nullptr;
  }

  G4int maxNumOfShells = transitionManager->NumberOfReachableAugerShells(Z);  
  
  const G4AugerTransition* refAugerTransition = 
    transitionManager->ReachableAugerShell(Z,maxNumOfShells-1);

  // This loop gives to shellNum the value of the index of shellId
  // in the vector storing the list of the vacancies in the variuos shells 
  // that can originate a NON-radiative transition
  G4int shellNum = 0;
    
  if ( shellId <= refAugerTransition->FinalShellId() ) 
    // "FinalShellId" is final from the point of view of the electron 
    // who makes the transition, 
    // being the Id of the shell in which there is a vacancy
    {
      G4int pippo = transitionManager->ReachableAugerShell(Z,shellNum)->FinalShellId();
      if (shellId != pippo ) {
	do { 
	  ++shellNum;
 	  if(shellNum == maxNumOfShells)
 	    {
	      // G4cout << "No Auger transition found" << G4endl; //debug
	      return 0;
 	    }
	}
 	while (shellId != (transitionManager->ReachableAugerShell(Z,shellNum)->FinalShellId()) );
      }

      // Now we have that shellnum is the shellIndex of the shell named ShellId
      //      G4cout << " the index of the shell is: "<<shellNum<<G4endl;
      // But we have now to select two shells: one for the transition, 
      // and another for the auger emission.
      G4int transitionLoopShellIndex = 0;      
      G4double partSum = 0;
      const G4AugerTransition* anAugerTransition = 
	transitionManager->ReachableAugerShell(Z,shellNum);

      G4int transitionSize = (G4int)
	(anAugerTransition->TransitionOriginatingShellIds())->size();
      while (transitionLoopShellIndex < transitionSize) {

        std::vector<G4int>::const_iterator pos = 
	  anAugerTransition->TransitionOriginatingShellIds()->cbegin();

        G4int transitionLoopShellId = *(pos+transitionLoopShellIndex);
        G4int numberOfPossibleAuger = (G4int)
	  (anAugerTransition->AugerTransitionProbabilities(transitionLoopShellId))->size();
        G4int augerIndex = 0;
      
	if (augerIndex < numberOfPossibleAuger) {
	  do 
	    {
	      G4double thisProb = anAugerTransition->AugerTransitionProbability(augerIndex, 
										transitionLoopShellId);
	      partSum += thisProb;
	      augerIndex++;
	      
	    } while (augerIndex < numberOfPossibleAuger);
        }
        ++transitionLoopShellIndex;
      }
     
      G4double totalVacancyAugerProbability = partSum;

      //And now we start to select the right auger transition and emission
      G4int transitionRandomShellIndex = 0;
      G4int transitionRandomShellId = 1;
      G4int augerIndex = 0;
      partSum = 0; 
      G4double partialProb = G4UniformRand();
      
      G4int numberOfPossibleAuger = 0;      
      G4bool foundFlag = false;

      while (transitionRandomShellIndex < transitionSize) {

        std::vector<G4int>::const_iterator pos = 
	  anAugerTransition->TransitionOriginatingShellIds()->begin();

        transitionRandomShellId = *(pos+transitionRandomShellIndex);
        
	augerIndex = 0;
	numberOfPossibleAuger = (G4int)(anAugerTransition-> 
				 AugerTransitionProbabilities(transitionRandomShellId))->size();

        while (augerIndex < numberOfPossibleAuger) {
	  G4double thisProb =anAugerTransition->AugerTransitionProbability(augerIndex, 
									   transitionRandomShellId);

          partSum += thisProb;
          
          if (partSum >= (partialProb*totalVacancyAugerProbability) ) { // was /
	    foundFlag = true;
	    break;
	  }
          augerIndex++;
        }
        if (partSum >= (partialProb*totalVacancyAugerProbability) ) {break;} // was /
        ++transitionRandomShellIndex;
      }

      // Now we have the index of the shell from wich comes the auger electron (augerIndex), 
      // and the id of the shell, from which the transition e- come (transitionRandomShellid)
      // If no Transition has been found, 0 is returned.  
      if (!foundFlag) {
	return nullptr;
      } 
      
      // Isotropic angular distribution for the outcoming e-
      G4double newcosTh = 1.-2.*G4UniformRand();
      G4double newsinTh = std::sqrt(1.-newcosTh*newcosTh);
      G4double newPhi = twopi*G4UniformRand();
      
      G4double xDir = newsinTh*std::sin(newPhi);
      G4double yDir = newsinTh*std::cos(newPhi);
      G4double zDir = newcosTh;
      
      G4ThreeVector newElectronDirection(xDir,yDir,zDir);
      
      // energy of the auger electron emitted            
      G4double transitionEnergy = 
	anAugerTransition->AugerTransitionEnergy(augerIndex, transitionRandomShellId);
      
      if (transitionEnergy < minElectronEnergy) {
	return nullptr;
      }

      // This is the shell where the new vacancy is: it is the same
      // shell where the electron came from
      newShellId = transitionRandomShellId;
      
      //Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)
      if (IsAugerCascadeActive())
	{
	  vacancyArray.push_back(newShellId);
	  vacancyArray.push_back(anAugerTransition->AugerOriginatingShellId(augerIndex,transitionRandomShellId));
	}
     
      return new G4DynamicParticle(G4Electron::Electron(), 
				   newElectronDirection,
				   transitionEnergy);
    }
  else 
    {
      return nullptr;
    }
}
