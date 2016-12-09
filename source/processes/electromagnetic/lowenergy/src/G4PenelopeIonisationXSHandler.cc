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
// $Id: G4PenelopeIonisationXSHandler.cc 99415 2016-09-21 09:05:43Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// --------
// 08 Mar 2012   L Pandola    First complete implementation
//

#include "G4PenelopeIonisationXSHandler.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeOscillator.hh"
#include "G4PenelopeCrossSection.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh" 

G4PenelopeIonisationXSHandler::G4PenelopeIonisationXSHandler(size_t nb)
  :XSTableElectron(0),XSTablePositron(0),
   theDeltaTable(0),energyGrid(0)
{
  nBins = nb;
  G4double LowEnergyLimit = 100.0*eV;
  G4double HighEnergyLimit = 100.0*GeV;
  oscManager = G4PenelopeOscillatorManager::GetOscillatorManager();
  XSTableElectron = new 
    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;
  XSTablePositron = new 
    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;

  theDeltaTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
  energyGrid = new G4PhysicsLogVector(LowEnergyLimit,
				      HighEnergyLimit, 
				      nBins-1); //one hidden bin is added

  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4PenelopeIonisationXSHandler::~G4PenelopeIonisationXSHandler()
{
  if (XSTableElectron)
    {
      for (auto& item : (*XSTableElectron))
	{	  
	  //G4PenelopeCrossSection* tab = i->second;
	  delete item.second;
	}
      delete XSTableElectron;
      XSTableElectron = nullptr;
    }

  if (XSTablePositron)
    {
      for (auto& item : (*XSTablePositron))
	{
	  //G4PenelopeCrossSection* tab = i->second;
	  delete item.second;
	}
      delete XSTablePositron;
      XSTablePositron = nullptr;
    }
  if (theDeltaTable)
    {      
      for (auto& item : (*theDeltaTable))
 	delete item.second;     
      delete theDeltaTable;
      theDeltaTable = nullptr;
    } 
  if (energyGrid)
    delete energyGrid;
  
  if (verboseLevel > 2)
    G4cout << "G4PenelopeIonisationXSHandler. Tables have been cleared" 
	   << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4PenelopeCrossSection*  
G4PenelopeIonisationXSHandler::GetCrossSectionTableForCouple(const G4ParticleDefinition* part,
								const G4Material* mat,
								const G4double cut) const
{
  if (part != G4Electron::Electron() && part != G4Positron::Positron())
    {
      G4ExceptionDescription ed;
      ed << "Invalid particle: " << part->GetParticleName() << G4endl;
      G4Exception("G4PenelopeIonisationXSHandler::GetCrossSectionTableForCouple()",
		  "em0001",FatalException,ed);
      return nullptr;
    }

  if (part == G4Electron::Electron())
    {
      if (!XSTableElectron)
	{
	  G4Exception("G4PenelopeIonisationXSHandler::GetCrossSectionTableForCouple()",
		      "em0028",FatalException,  
		      "The Cross Section Table for e- was not initialized correctly!");
	  return nullptr;
	}
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (XSTableElectron->count(theKey)) //table already built	
	return XSTableElectron->find(theKey)->second;
      else	 
        return nullptr;	 
    }
  
  if (part == G4Positron::Positron())
    {
      if (!XSTablePositron)
	{
	  G4Exception("G4PenelopeIonisationXSHandler::GetCrossSectionTableForCouple()",
		      "em0028",FatalException,  
		      "The Cross Section Table for e+ was not initialized correctly!");
	  return nullptr;
	}
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (XSTablePositron->count(theKey)) //table already built	
	return XSTablePositron->find(theKey)->second;
      else
        return nullptr;  
   }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationXSHandler::BuildXSTable(const G4Material* mat,G4double cut,
						 const G4ParticleDefinition* part,
						 G4bool isMaster)
{
  //Just to check
  if (!isMaster)
    G4Exception("G4PenelopeIonisationXSHandler::BuildXSTable()",
                "em0100",FatalException,"Worker thread in this method");
  
  //
  //This method fills the G4PenelopeCrossSection containers for electrons or positrons
  //and for the given material/cut couple. The calculation is done as sum over the 
  //individual shells.
  //Equivalent of subroutines EINaT and PINaT of Penelope
  //
  if (verboseLevel > 2)
    {
      G4cout << "G4PenelopeIonisationXSHandler: going to build cross section table " << G4endl;
      G4cout << "for " << part->GetParticleName() << " in " << mat->GetName() << G4endl;
      G4cout << "Cut= " << cut/keV << " keV" << G4endl;
    }

  std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
  //Check if the table already exists
  if (part == G4Electron::Electron())
    {
      if (XSTableElectron->count(theKey)) //table already built	
	return;
    }
  if (part == G4Positron::Positron())
    {
      if (XSTablePositron->count(theKey)) //table already built	
	return;
    }
  
  //check if the material has been built
  if (!(theDeltaTable->count(mat)))
    BuildDeltaTable(mat);


  //Tables have been already created (checked by GetCrossSectionTableForCouple)
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  size_t numberOfOscillators = theTable->size();

  if (energyGrid->GetVectorLength() != nBins) 
    {
      G4ExceptionDescription ed;
      ed << "Energy Grid looks not initialized" << G4endl;
      ed << nBins << " " << energyGrid->GetVectorLength() << G4endl;
      G4Exception("G4PenelopeIonisationXSHandler::BuildXSTable()",
		  "em2030",FatalException,ed);
    }

  G4PenelopeCrossSection* XSEntry = new G4PenelopeCrossSection(nBins,numberOfOscillators);
 
  //loop on the energy grid
  for (size_t bin=0;bin<nBins;bin++)
    {
       G4double energy = energyGrid->GetLowEdgeEnergy(bin);
       G4double XH0=0, XH1=0, XH2=0;
       G4double XS0=0, XS1=0, XS2=0;
   
       //oscillator loop
       for (size_t iosc=0;iosc<numberOfOscillators;iosc++)
	 {
	   G4DataVector* tempStorage = 0;

	   G4PenelopeOscillator* theOsc = (*theTable)[iosc];
	   G4double delta = GetDensityCorrection(mat,energy);
	   if (part == G4Electron::Electron())	     
	     tempStorage = ComputeShellCrossSectionsElectron(theOsc,energy,cut,delta);
	   else if (part == G4Positron::Positron())
	     tempStorage = ComputeShellCrossSectionsPositron(theOsc,energy,cut,delta);
	   //check results are all right
	   if (!tempStorage)
	     {
	       G4ExceptionDescription ed;
	       ed << "Problem in calculating the shell XS for shell # " 
		  << iosc << G4endl;
	       G4Exception("G4PenelopeIonisationXSHandler::BuildXSTable()",
			   "em2031",FatalException,ed);			   
	       delete XSEntry;
	       return;
	     }
	   if (tempStorage->size() != 6)
	     {
	       G4ExceptionDescription ed;
	       ed << "Problem in calculating the shell XS " << G4endl;
	       ed << "Result has dimension " << tempStorage->size() << " instead of 6" << G4endl;
	       G4Exception("G4PenelopeIonisationXSHandler::BuildXSTable()",
			   "em2031",FatalException,ed);	
	     }
	   G4double stre = theOsc->GetOscillatorStrength();

	   XH0 += stre*(*tempStorage)[0];
	   XH1 += stre*(*tempStorage)[1];
	   XH2 += stre*(*tempStorage)[2];
	   XS0 += stre*(*tempStorage)[3];
	   XS1 += stre*(*tempStorage)[4];
	   XS2 += stre*(*tempStorage)[5];
	   XSEntry->AddShellCrossSectionPoint(bin,iosc,energy,stre*(*tempStorage)[0]);
	   if (tempStorage)
	     {
	       delete tempStorage;
	       tempStorage = 0;
	     }
	 }       
       XSEntry->AddCrossSectionPoint(bin,energy,XH0,XH1,XH2,XS0,XS1,XS2);
    }
  //Do (only once) the final normalization
  XSEntry->NormalizeShellCrossSections();

  //Insert in the appropriate table
  if (part == G4Electron::Electron())      
    XSTableElectron->insert(std::make_pair(theKey,XSEntry));
  else if (part == G4Positron::Positron())
    XSTablePositron->insert(std::make_pair(theKey,XSEntry));
  else
    delete XSEntry;
  
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationXSHandler::GetDensityCorrection(const G4Material* mat,
							     const G4double energy) const
{
  G4double result = 0;
  if (!theDeltaTable)
    {
      G4Exception("G4PenelopeIonisationXSHandler::GetDensityCorrection()",
		  "em2032",FatalException,
		  "Delta Table not initialized. Was Initialise() run?");
      return 0;
    }
  if (energy <= 0*eV)
    {
      G4cout << "G4PenelopeIonisationXSHandler::GetDensityCorrection()" << G4endl;
      G4cout << "Invalid energy " << energy/eV << " eV " << G4endl;
      return 0;
    }
  G4double logene = std::log(energy);
 
  if (theDeltaTable->count(mat))
    {
      const G4PhysicsFreeVector* vec = theDeltaTable->find(mat)->second;
      result = vec->Value(logene); //the table has delta vs. ln(E)      
    }
  else
    {
      G4ExceptionDescription ed;      
      ed << "Unable to build table for " << mat->GetName() << G4endl;
      G4Exception("G4PenelopeIonisationXSHandler::GetDensityCorrection()",
		  "em2033",FatalException,ed);
    }

  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationXSHandler::BuildDeltaTable(const G4Material* mat)
{
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  G4double plasmaSq = oscManager->GetPlasmaEnergySquared(mat);
  G4double totalZ = oscManager->GetTotalZ(mat);
  size_t numberOfOscillators = theTable->size();

  if (energyGrid->GetVectorLength() != nBins) 
    {
      G4ExceptionDescription ed;
      ed << "Energy Grid for Delta table looks not initialized" << G4endl;
      ed << nBins << " " << energyGrid->GetVectorLength() << G4endl;
      G4Exception("G4PenelopeIonisationXSHandler::BuildDeltaTable()",
		  "em2030",FatalException,ed);
    }

  G4PhysicsFreeVector* theVector = new G4PhysicsFreeVector(nBins);

  //loop on the energy grid
  for (size_t bin=0;bin<nBins;bin++)
    {
      G4double delta = 0.;
      G4double energy = energyGrid->GetLowEdgeEnergy(bin);

      //Here calculate delta
      G4double gam = 1.0+(energy/electron_mass_c2);
      G4double gamSq = gam*gam;

      G4double TST = totalZ/(gamSq*plasmaSq);
      G4double wl2 = 0;
      G4double fdel = 0;

      //loop on oscillators
      for (size_t i=0;i<numberOfOscillators;i++)
	{
	  G4PenelopeOscillator* theOsc = (*theTable)[i];
	  G4double wri = theOsc->GetResonanceEnergy();
	  fdel += theOsc->GetOscillatorStrength()/(wri*wri+wl2);
	}      
      if (fdel >= TST) //if fdel < TST, delta = 0
	{
	  //get last oscillator
	  G4PenelopeOscillator* theOsc = (*theTable)[numberOfOscillators-1]; 
	  wl2 = theOsc->GetResonanceEnergy()*theOsc->GetResonanceEnergy();
	  
	  //First iteration
	  G4bool loopAgain = false;
	  do
	    {
	      loopAgain = false;
	      wl2 += wl2;
	      fdel = 0.;
	      for (size_t i=0;i<numberOfOscillators;i++)
		{
		  G4PenelopeOscillator* theOscLocal1 = (*theTable)[i];
		  G4double wri = theOscLocal1->GetResonanceEnergy();
		  fdel += theOscLocal1->GetOscillatorStrength()/(wri*wri+wl2);
		}
	      if (fdel > TST)
		loopAgain = true;
	    }while(loopAgain);

	  G4double wl2l = 0;
	  G4double wl2u = wl2;
	  //second iteration
	  do
	    {	     
	      loopAgain = false;
	      wl2 = 0.5*(wl2l+wl2u);
	      fdel = 0;
	      for (size_t i=0;i<numberOfOscillators;i++)
		{
		  G4PenelopeOscillator* theOscLocal2 = (*theTable)[i];
		  G4double wri = theOscLocal2->GetResonanceEnergy();
		  fdel += theOscLocal2->GetOscillatorStrength()/(wri*wri+wl2);
		}
	      if (fdel > TST)
		wl2l = wl2;
	      else
		wl2u = wl2;
	      if ((wl2u-wl2l)>1e-12*wl2)
		loopAgain = true;
	    }while(loopAgain);
	  
	  //Eventually get density correction
	  delta = 0.;
	  for (size_t i=0;i<numberOfOscillators;i++)
	    {
	      G4PenelopeOscillator* theOscLocal3 = (*theTable)[i];
	      G4double wri = theOscLocal3->GetResonanceEnergy();
	      delta += theOscLocal3->GetOscillatorStrength()*
		std::log(1.0+(wl2/(wri*wri)));	  	     
	    }
	  delta = (delta/totalZ)-wl2/(gamSq*plasmaSq);
	}
      energy = std::max(1e-9*eV,energy); //prevents log(0)
      theVector->PutValue(bin,std::log(energy),delta);
    }
  theDeltaTable->insert(std::make_pair(mat,theVector));
  return;
}
							
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4DataVector* G4PenelopeIonisationXSHandler::ComputeShellCrossSectionsElectron(G4PenelopeOscillator* theOsc,
										  G4double energy,
										  G4double cut,
										  G4double delta)
{
  //
  //This method calculates the hard and soft cross sections (H0-H1-H2-S0-S1-S2) for 
  //the given oscillator/cut and at the given energy.
  //It returns a G4DataVector* with 6 entries (H0-H1-H2-S0-S1-S2)
  //Equivalent of subroutines EINaT1 of Penelope
  //
  // Results are _per target electron_
  //
  G4DataVector* result = new G4DataVector();
  for (size_t i=0;i<6;i++)
    result->push_back(0.);
  G4double ionEnergy = theOsc->GetIonisationEnergy();
  
  //return a set of zero's if the energy it too low to excite the current oscillator
  if (energy < ionEnergy)
    return result;

  G4double H0=0.,H1=0.,H2=0.;
  G4double S0=0.,S1=0.,S2=0.;

  //Define useful constants to be used in the calculation
  G4double gamma = 1.0+energy/electron_mass_c2;
  G4double gammaSq = gamma*gamma;
  G4double beta = (gammaSq-1.0)/gammaSq;
  G4double pielr2 = pi*classic_electr_radius*classic_electr_radius; //pi*re^2
  G4double constant = pielr2*2.0*electron_mass_c2/beta;
  G4double XHDT0 = std::log(gammaSq)-beta;

  G4double cpSq = energy*(energy+2.0*electron_mass_c2);
  G4double cp = std::sqrt(cpSq);
  G4double amol = (energy/(energy+electron_mass_c2))*(energy/(energy+electron_mass_c2));

  //
  // Distant interactions
  //
  G4double resEne = theOsc->GetResonanceEnergy();
  G4double cutoffEne = theOsc->GetCutoffRecoilResonantEnergy();
  if (energy > resEne)
    {
      G4double cp1Sq = (energy-resEne)*(energy-resEne+2.0*electron_mass_c2);
      G4double cp1 = std::sqrt(cp1Sq);
      
      //Distant longitudinal interactions
      G4double QM = 0;
      if (resEne > 1e-6*energy)
	QM = std::sqrt((cp-cp1)*(cp-cp1)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
      else
	{
	  QM = resEne*resEne/(beta*2.0*electron_mass_c2);
	  QM = QM*(1.0-0.5*QM/electron_mass_c2);
	}
      G4double SDL1 = 0;
      if (QM < cutoffEne)
	SDL1 = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)));
      
      //Distant transverse interactions
      if (SDL1)
	{
	  G4double SDT1 = std::max(XHDT0-delta,0.0);
	  G4double SD1 = SDL1+SDT1;
	  if (cut > resEne)
	    {
	      S1 = SD1; //XS1
	      S0 = SD1/resEne; //XS0
	      S2 = SD1*resEne; //XS2
	    }
	  else
	    {
	      H1 = SD1; //XH1
	      H0 = SD1/resEne; //XH0
	      H2 = SD1*resEne; //XH2
	    }
	}
    }
  //
  // Close collisions (Moller's cross section)
  //
  G4double wl = std::max(cut,cutoffEne);
  G4double ee = energy + ionEnergy;
  G4double wu = 0.5*ee;
  if (wl < wu-(1e-5*eV))
    {
      H0 += (1.0/(ee-wu)) - (1.0/(ee-wl)) - (1.0/wu) + (1.0/wl) + 
	(1.0-amol)*std::log(((ee-wu)*wl)/((ee-wl)*wu))/ee + 
	amol*(wu-wl)/(ee*ee);
      H1 += std::log(wu/wl)+(ee/(ee-wu))-(ee/(ee-wl)) + 
	(2.0-amol)*std::log((ee-wu)/(ee-wl)) + 
	amol*(wu*wu-wl*wl)/(2.0*ee*ee);
      H2 += (2.0-amol)*(wu-wl)+(wu*(2.0*ee-wu)/(ee-wu)) - 
	(wl*(2.0*ee-wl)/(ee-wl)) + 
	(3.0-amol)*ee*std::log((ee-wu)/(ee-wl)) + 	
	amol*(wu*wu*wu-wl*wl*wl)/(3.0*ee*ee);
      wu = wl;
    }
  wl = cutoffEne;
  
  if (wl > wu-(1e-5*eV))
    {
      (*result)[0] = constant*H0;
      (*result)[1] = constant*H1;
      (*result)[2] = constant*H2;
      (*result)[3] = constant*S0;
      (*result)[4] = constant*S1;
      (*result)[5] = constant*S2;
      return result;
    }

  S0 += (1.0/(ee-wu))-(1.0/(ee-wl)) - (1.0/wu) + (1.0/wl) + 
    (1.0-amol)*std::log(((ee-wu)*wl)/((ee-wl)*wu))/ee +
    amol*(wu-wl)/(ee*ee);
  S1 += std::log(wu/wl)+(ee/(ee-wu))-(ee/(ee-wl)) + 
    (2.0-amol)*std::log((ee-wu)/(ee-wl)) + 
    amol*(wu*wu-wl*wl)/(2.0*ee*ee);
  S2 += (2.0-amol)*(wu-wl)+(wu*(2.0*ee-wu)/(ee-wu)) - 
    (wl*(2.0*ee-wl)/(ee-wl)) + 
    (3.0-amol)*ee*std::log((ee-wu)/(ee-wl)) + 
    amol*(wu*wu*wu-wl*wl*wl)/(3.0*ee*ee);

  (*result)[0] = constant*H0;
  (*result)[1] = constant*H1;
  (*result)[2] = constant*H2;
  (*result)[3] = constant*S0;
  (*result)[4] = constant*S1;
  (*result)[5] = constant*S2;
  return result;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4DataVector* G4PenelopeIonisationXSHandler::ComputeShellCrossSectionsPositron(G4PenelopeOscillator* theOsc,
										  G4double energy,
										  G4double cut,
										  G4double delta)
{
  //
  //This method calculates the hard and soft cross sections (H0-H1-H2-S0-S1-S2) for 
  //the given oscillator/cut and at the given energy.
  //It returns a G4DataVector* with 6 entries (H0-H1-H2-S0-S1-S2)
  //Equivalent of subroutines PINaT1 of Penelope
  //
  // Results are _per target electron_
  //
  G4DataVector* result = new G4DataVector();
  for (size_t i=0;i<6;i++)
    result->push_back(0.);
  G4double ionEnergy = theOsc->GetIonisationEnergy();
  
  //return a set of zero's if the energy it too low to excite the current oscillator
  if (energy < ionEnergy)
    return result;

  G4double H0=0.,H1=0.,H2=0.;
  G4double S0=0.,S1=0.,S2=0.;

  //Define useful constants to be used in the calculation
  G4double gamma = 1.0+energy/electron_mass_c2;
  G4double gammaSq = gamma*gamma;
  G4double beta = (gammaSq-1.0)/gammaSq;
  G4double pielr2 = pi*classic_electr_radius*classic_electr_radius; //pi*re^2
  G4double constant = pielr2*2.0*electron_mass_c2/beta;
  G4double XHDT0 = std::log(gammaSq)-beta;

  G4double cpSq = energy*(energy+2.0*electron_mass_c2);
  G4double cp = std::sqrt(cpSq);
  G4double amol = (energy/(energy+electron_mass_c2))*(energy/(energy+electron_mass_c2));
  G4double g12 = (gamma+1.0)*(gamma+1.0);
  //Bhabha coefficients
  G4double bha1 = amol*(2.0*g12-1.0)/(gammaSq-1.0);
  G4double bha2 = amol*(3.0+1.0/g12);
  G4double bha3 = amol*2.0*gamma*(gamma-1.0)/g12;
  G4double bha4 = amol*(gamma-1.0)*(gamma-1.0)/g12;

  //
  // Distant interactions
  //
  G4double resEne = theOsc->GetResonanceEnergy();
  G4double cutoffEne = theOsc->GetCutoffRecoilResonantEnergy();
  if (energy > resEne)
    {
      G4double cp1Sq = (energy-resEne)*(energy-resEne+2.0*electron_mass_c2);
      G4double cp1 = std::sqrt(cp1Sq);
      
      //Distant longitudinal interactions
      G4double QM = 0;
      if (resEne > 1e-6*energy)
	QM = std::sqrt((cp-cp1)*(cp-cp1)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
      else
	{
	  QM = resEne*resEne/(beta*2.0*electron_mass_c2);
	  QM = QM*(1.0-0.5*QM/electron_mass_c2);
	}
      G4double SDL1 = 0;
      if (QM < cutoffEne)
	SDL1 = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)));
      
      //Distant transverse interactions
      if (SDL1)
	{
	  G4double SDT1 = std::max(XHDT0-delta,0.0);
	  G4double SD1 = SDL1+SDT1;
	  if (cut > resEne)
	    {
	      S1 = SD1; //XS1
	      S0 = SD1/resEne; //XS0
	      S2 = SD1*resEne; //XS2
	    }
	  else
	    {
	      H1 = SD1; //XH1
	      H0 = SD1/resEne; //XH0
	      H2 = SD1*resEne; //XH2
	    }
	}
    }

  //
  // Close collisions (Bhabha's cross section)
  //
  G4double wl = std::max(cut,cutoffEne);
  G4double wu = energy; 
  G4double energySq = energy*energy;
  if (wl < wu-(1e-5*eV))
    {
      G4double wlSq = wl*wl;
      G4double wuSq = wu*wu;
      H0 += (1.0/wl) - (1.0/wu)- bha1*std::log(wu/wl)/energy  
	+ bha2*(wu-wl)/energySq  
	- bha3*(wuSq-wlSq)/(2.0*energySq*energy)
	+ bha4*(wuSq*wu-wlSq*wl)/(3.0*energySq*energySq);
      H1 += std::log(wu/wl) - bha1*(wu-wl)/energy
	+ bha2*(wuSq-wlSq)/(2.0*energySq)
	- bha3*(wuSq*wu-wlSq*wl)/(3.0*energySq*energy)
	+ bha4*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energySq);
      H2 += wu - wl - bha1*(wuSq-wlSq)/(2.0*energy)
	+ bha2*(wuSq*wu-wlSq*wl)/(3.0*energySq)
	- bha3*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energy)
	+ bha4*(wuSq*wuSq*wu-wlSq*wlSq*wl)/(5.0*energySq*energySq);
      wu = wl;
    }
  wl = cutoffEne;
  
  if (wl > wu-(1e-5*eV))
    {
      (*result)[0] = constant*H0;
      (*result)[1] = constant*H1;
      (*result)[2] = constant*H2;
      (*result)[3] = constant*S0;
      (*result)[4] = constant*S1;
      (*result)[5] = constant*S2;
      return result;
    }

  G4double wlSq = wl*wl;
  G4double wuSq = wu*wu;

  S0 += (1.0/wl) - (1.0/wu) - bha1*std::log(wu/wl)/energy 
    + bha2*(wu-wl)/energySq  
    - bha3*(wuSq-wlSq)/(2.0*energySq*energy)
    + bha4*(wuSq*wu-wlSq*wl)/(3.0*energySq*energySq);

  S1 += std::log(wu/wl) - bha1*(wu-wl)/energy
    + bha2*(wuSq-wlSq)/(2.0*energySq)
    - bha3*(wuSq*wu-wlSq*wl)/(3.0*energySq*energy)
    + bha4*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energySq);

  S2 += wu - wl - bha1*(wuSq-wlSq)/(2.0*energy)
    + bha2*(wuSq*wu-wlSq*wl)/(3.0*energySq)
    - bha3*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energy)
    + bha4*(wuSq*wuSq*wu-wlSq*wlSq*wl)/(5.0*energySq*energySq);

 (*result)[0] = constant*H0;
 (*result)[1] = constant*H1;
 (*result)[2] = constant*H2;
 (*result)[3] = constant*S0;
 (*result)[4] = constant*S1;
 (*result)[5] = constant*S2;

 return result;
}
