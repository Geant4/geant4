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
// $Id: G4PenelopeIonisationCrossSection.cc 99415 2016-09-21 09:05:43Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// --------
// 14 Mar 2012   L Pandola    First complete implementation for e-
//
//
//! NOTICE: working only for e- at the moment (no interface available for
//!         e+)
//
#include "G4PenelopeOscillatorManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PenelopeIonisationCrossSection.hh"
#include "G4PenelopeIonisationXSHandler.hh"
#include "G4PenelopeCrossSection.hh"
#include "G4Electron.hh"
#include "G4AtomicTransitionManager.hh"

G4PenelopeIonisationCrossSection::G4PenelopeIonisationCrossSection() : 
  G4VhShellCrossSection("Penelope"),shellIDTable(0),
  theCrossSectionHandler(0)
{ 
  oscManager = G4PenelopeOscillatorManager::GetOscillatorManager();
  nMaxLevels = 9;

  // Verbosity scale:
  // 0 = nothing
  // 1 = calculation of cross sections, file openings, sampling of atoms
  // 2 = entering in methods
  verboseLevel = 0;

  fLowEnergyLimit  = 10.0*eV;
  fHighEnergyLimit = 100.0*GeV;

  transitionManager = G4AtomicTransitionManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4PenelopeIonisationCrossSection::~G4PenelopeIonisationCrossSection()
{
  if (theCrossSectionHandler)
    delete theCrossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationCrossSection::CrossSection(G4int Z,
							G4AtomicShellEnumerator shell,
							G4double incidentEnergy,
							G4double ,
							const G4Material* material)
{
  if (verboseLevel > 1)    
    G4cout << "Entering in method G4PenelopeIonisationCrossSection::CrossSection()" << G4endl;
        
  G4double cross = 0.;

  //Material pointer is not available
  if (!material)
    {
      //CRASH!
      G4ExceptionDescription ed;
      ed << "The method has been called with a null G4Material pointer" << G4endl;
      G4Exception("G4PenelopeIonisationCrossSection::CrossSection()","em2042",
		  FatalException,ed);
      return cross;
    }

  if (!theCrossSectionHandler)
    theCrossSectionHandler = new G4PenelopeIonisationXSHandler();

  theCrossSectionHandler->BuildXSTable(material,0.,G4Electron::Electron());

  G4int nmax = std::min(nMaxLevels,transitionManager->NumberOfShells(Z));

  if(G4int(shell) < nmax && 
     incidentEnergy >= fLowEnergyLimit && incidentEnergy <= fHighEnergyLimit) 
    {
      //The shells in Penelope are organized per *material*, rather than per 
      //element, so given a material one has to find the proper index for the 
      //given Z and shellID. An appropriate lookup table is used to avoid 
      //recalculation.
      G4int index = FindShellIDIndex(material,Z,shell);

      //Index is not available!
      if (index < 0)
	return cross;

      const G4PenelopeCrossSection*  theXS = 
	theCrossSectionHandler->GetCrossSectionTableForCouple(G4Electron::Electron(),
							      material,
							      0.);

      //Cross check that everything is fine:
      G4PenelopeOscillator* theOsc = oscManager->GetOscillatorIonisation(material,index);
      if (theOsc->GetParentZ() != Z || theOsc->GetShellFlag()-1 != G4int(shell))
	{
	  //something went wrong!
	  G4ExceptionDescription ed;
	  ed << "There is something wrong here: it looks like the index is wrong" << G4endl;
	  ed << "Requested: shell " << G4int(shell) << " and Z = " << Z << G4endl;
	  ed << "Retrieved: " << theOsc->GetShellFlag()-1 << " and Z = " << theOsc->GetParentZ() << G4endl;	  
	  G4Exception("G4PenelopeIonisationCrossSection::CrossSection()","em2043",
		      JustWarning,ed);
	  return cross;
	}


      G4double crossPerMolecule = (theXS) ? theXS->GetShellCrossSection(index,incidentEnergy) : 0.;

      //Now it must be converted in cross section per atom. I need the number of 
      //atoms of the given Z per molecule.
      G4double atomsPerMolec = oscManager->GetNumberOfZAtomsPerMolecule(material,Z);
      if (atomsPerMolec)
	cross = crossPerMolecule/atomsPerMolec;
     

      if (verboseLevel > 0)
	{
	  G4cout << "Cross section of shell " << G4int(shell) << " and Z= " << Z;
	  G4cout << " of material: " << material->GetName() << " and energy = " << incidentEnergy/keV << " keV" << G4endl;
	  G4cout << "--> " << cross/barn << " barn" << G4endl;
	  G4cout << "Shell binding energy: " << theOsc->GetIonisationEnergy()/eV << " eV;" ;
	  G4cout << " resonance energy: " << theOsc->GetResonanceEnergy()/eV << "eV" << G4endl;
	   if (verboseLevel > 2)
	     {
	       G4cout << "Cross section per molecule: " << crossPerMolecule/barn << " barn" << G4endl;
	       G4cout << "Atoms " << Z << " per molecule: " << atomsPerMolec << G4endl;
	     }
	}     
    }  
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
std::vector<G4double> 
G4PenelopeIonisationCrossSection::GetCrossSection(G4int Z,
						  G4double kinEnergy,
						  G4double, G4double, 
						  const G4Material* mat)
{
  G4int nmax = std::min(nMaxLevels,transitionManager->NumberOfShells(Z));
  std::vector<G4double> vec(nmax,0.0); 
  for(G4int i=0; i<nmax; ++i) {
    vec[i] = CrossSection(Z, G4AtomicShellEnumerator(i), kinEnergy,0.,mat); 
  }
  return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4double> 
G4PenelopeIonisationCrossSection::Probabilities(G4int Z,
                                                 G4double kinEnergy,
                                                 G4double,
                                                 G4double,
                                                 const G4Material* mat)
{
  std::vector<G4double> vec = GetCrossSection(Z, kinEnergy,0,0,mat);
  size_t n = vec.size();
  size_t i=0;
  G4double sum = 0.0;
  for(i=0; i<n; ++i) { sum += vec[i]; }
  if(sum > 0.0) { 
    sum = 1.0/sum; 
    for(i=0; i<n; ++i) { vec[i] = vec[i]*sum; }
  }
  return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4PenelopeIonisationCrossSection::FindShellIDIndex(const G4Material* mat,
							 G4int Z,
							 G4AtomicShellEnumerator shell)
{
   if (verboseLevel > 1)    
     G4cout << "Entering in method G4PenelopeIonisationCrossSection::FindShellIDIndex()" << G4endl;

  if (!shellIDTable)
    shellIDTable = new std::map< std::pair<const G4Material*,G4int>, G4DataVector*>;
 
  std::pair<const G4Material*,G4int> theKey = std::make_pair(mat,Z);
  G4int result = -1;
  G4int ishell = G4int(shell);
  
  if (shellIDTable->count(theKey)) //table already built, and containing the element
    {   
      if (verboseLevel > 2)
	G4cout << "FindShellIDIndex: Table already built for " << mat->GetName() << G4endl;
      G4DataVector* theVec = shellIDTable->find(theKey)->second;
         
      if (ishell>=0 && ishell < (G4int) theVec->size()) //check we are not off-boundary
	result = (G4int) (*theVec)[ishell];    
      else
	{
	  G4ExceptionDescription ed;
	  ed << "Shell ID: " << ishell << " not available for material " << mat->GetName() << " and Z = " << 
	    Z << G4endl;
	  G4Exception("G4PenelopeIonisationCrossSection::FindShellIDIndex()","em2041",JustWarning,
		      ed);
	  return -1;
	}          
    }
  else
    {
      if (verboseLevel > 2)
	G4cout << "FindShellIDIndex: Table to be built for " << mat->GetName() << G4endl;
      //Not contained: look for it
      G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
      size_t numberOfOscillators = theTable->size();
          
      //oscillator loop
      //initialize everything at -1
      G4DataVector* dat = new G4DataVector(nMaxLevels,-1);
      for (size_t iosc=0;iosc<numberOfOscillators;iosc++)
	{      
	  G4PenelopeOscillator* theOsc = (*theTable)[iosc];
	  //level is found!
	  if (theOsc->GetParentZ() == Z)
	    {
	      //individual shells relative to the given material	    
	      G4int shFlag = theOsc->GetShellFlag(); 	  
	      //Notice: GetShellFlag() starts from 1, the G4AtomicShellEnumerator from 0
	      if (shFlag < 30)
		(*dat)[shFlag-1] = (G4double) iosc; //index of the given shell
	      if ((shFlag-1) == ishell) // this is what we were looking for
		result = (G4int) iosc;
	    }      
	}
      shellIDTable->insert(std::make_pair(theKey,dat));
    }

  if (verboseLevel > 1)    
    G4cout << "Leaving method G4PenelopeIonisationCrossSection::FindShellIDIndex() with index = " << result << G4endl;

  return result;
  
}
