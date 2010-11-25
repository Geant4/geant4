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
// $Id: G4PenelopeGammaConversionModel.cc,v 1.7 2010-11-25 09:45:13 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 06 Oct 2008   L Pandola    Migration from process to model 
// 17 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - do not apply production threshold on level of the model
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they might be checked later on
//

#include "G4PenelopeGammaConversionModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopeGammaConversionModel::G4PenelopeGammaConversionModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),fTheScreeningRadii(0),crossSectionHandler(0),isInitialised(false)
{
  fIntrinsicLowEnergyLimit = 2.0*electron_mass_c2;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  fSmallEnergy = 1.1*MeV;

  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeGammaConversionModel::~G4PenelopeGammaConversionModel()
{  
  if (crossSectionHandler) delete crossSectionHandler;
  if (fTheScreeningRadii) delete fTheScreeningRadii;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeGammaConversionModel::Initialise(const G4ParticleDefinition*,
						const G4DataVector& )
{
  if (verboseLevel > 3)
    G4cout << "Calling  G4PenelopeGammaConversionModel::Initialise()" << G4endl;

  //Delete the old cross section handler, if necessary
  if (crossSectionHandler)
    {
      crossSectionHandler->Clear();
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }
  
  //Re-initialize cross section handler
  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Initialise(0,fIntrinsicLowEnergyLimit,HighEnergyLimit(),400);
  crossSectionHandler->Clear();
  G4String crossSectionFile = "penelope/pp-cs-pen-";
  crossSectionHandler->LoadData(crossSectionFile);
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata =
    crossSectionHandler->BuildMeanFreePathForMaterials();
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;

  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for PenelopeGammaConversion" << G4endl;

  if (verboseLevel > 0) { 
    G4cout << "Penelope Gamma Conversion model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeGammaConversionModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  //
  // Penelope model.
  // Cross section (including triplet production) read from database and managed 
  // through the G4CrossSectionHandler utility. Cross section data are from
  // M.J. Berger and J.H. Hubbel (XCOM), Report NBSIR 887-3598
  //
  
  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopePhotoElectricModel" << G4endl;

  G4int iZ = (G4int) Z;
  //  if (!crossSectionHandler) //VI: should not be checked in run time
  //  {
  //    G4cout << "G4PenelopeGammaConversionModel::ComputeCrossSectionPerAtom" << G4endl;
  //    G4cout << "The cross section handler is not correctly initialized" << G4endl;
  //    G4Exception();
  //  }
  G4double cs = crossSectionHandler->FindValue(iZ,energy);

  if (verboseLevel > 2)
    G4cout << "Gamma conversion cross section at " << energy/MeV << " MeV for Z=" << Z << 
      " = " << cs/barn << " barn" << G4endl;
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4PenelopeGammaConversionModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						  const G4MaterialCutsCouple* couple,
						  const G4DynamicParticle* aDynamicGamma,
						  G4double,
						  G4double)
{
  //
  // Penelope model.
  // Final state is sampled according to the Bethe-Heitler model with Coulomb
  // corrections, according to the semi-empirical model of
  //  J. Baro' et al., Radiat. Phys. Chem. 44 (1994) 531.
  //
  // The model uses the high energy Coulomb correction from 
  //  H. Davies et al., Phys. Rev. 93 (1954) 788
  // and atomic screening radii tabulated from 
  //  J.H. Hubbel et al., J. Phys. Chem. Ref. Data 9 (1980) 1023
  // for Z= 1 to 92. This managed in this model by the method
  // GetScreeningRadius(). 
  //
  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeGammaConversionModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();

  // Always kill primary
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->SetProposedKineticEnergy(0.);

  if (photonEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return ;
    }

  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  G4double eps ;
  G4double eki = electron_mass_c2 / photonEnergy ;

  // Do it fast if photon energy < 1.1 MeV
  if (photonEnergy < fSmallEnergy )
    {
      eps = eki + (1-2*eki) * G4UniformRand();
    }
  else
    {
      // Select randomly one element in the current material
      if (verboseLevel > 2)
	G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;
      //use crossSectionHandler instead of G4EmElementSelector because in this case 
      //the dimension of the table is equal to the dimension of the database 
      //(less interpolation errors)
      G4int Z_int = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);
      if (verboseLevel > 2)
	G4cout << "Selected Z = " << Z_int << G4endl;
      
      //Low energy and Coulomb corrections
      G4double Z=(G4double) Z_int;
      G4double ZAlpha = Z*fine_structure_const;
      G4double ScreenRadius = GetScreeningRadius(Z);
      G4double funct1=0,g0=0;
      G4double g1min=0,g2min=0;
      funct1 = 4.0*std::log(ScreenRadius);
      g0 = funct1-4*CoulombCorrection(ZAlpha)+LowEnergyCorrection(ZAlpha,eki); 
      G4double bmin = 2*eki*ScreenRadius;
      std::vector<G4double> ScreenFunctionValues = ScreenFunction(bmin);
      if (ScreenFunctionValues.size() != 2)
	{
	  G4cout << "G4PenelopeGammaConversionModel::SampleSecondaries" << G4endl;
	  G4cout << "ScreenFunction did not return 2 values! Something wrong! " << G4endl;
	  G4Exception();
	}
      g1min=g0+ScreenFunctionValues[0];
      g2min=g0+ScreenFunctionValues[1];
      G4double xr,a1,p1;
      xr=0.5-eki;
      a1=(2.0/3.0)*g1min*xr*xr;
      p1=a1/(a1+g2min);

      //Random sampling of eps
      G4double rand1,rand2,rand3,b;
      G4double g1;
   
      do{
        rand1 = G4UniformRand();
        if (rand1 < p1) {
	  rand2 = 2.0*G4UniformRand()-1.0;
	  if (rand2 < 0) {
	    eps = 0.5 - xr*std::pow(std::abs(rand2),(1./3.));
	  }
	  else
	    {
	      eps = 0.5 + xr*std::pow(rand2,(1./3.));
	    }
	  b = (eki*ScreenRadius)/(2*eps*(1.0-eps));
	  std::vector<G4double> ScreenFunctionSampling = ScreenFunction(b);
	  g1 = g0+ScreenFunctionSampling[0];
	  if (g1 < 0) g1=0;
	  rand3 = G4UniformRand()*g1min;
        }
        else
	  {
	    eps = eki+2.0*xr*G4UniformRand();
	    b = (eki*ScreenRadius)/(2*eps*(1.0-eps));
	    std::vector<G4double> ScreenFunctionSampling = ScreenFunction(b);
	    g1 = g0+ScreenFunctionSampling[1];
	    if (g1 < 0) g1=0; 
	    rand3 = G4UniformRand()*g2min;
	  }	   
      } while (rand3>g1);
    } //End of eps sampling

  G4double electronTotEnergy = eps*photonEnergy;
  G4double positronTotEnergy = (1.0-eps)*photonEnergy;
  
  // Scattered electron (positron) angles. ( Z - axis along the parent photon)

  //electron kinematics
  G4double costheta_el,costheta_po;
  G4double phi_el,phi_po;
  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ; 
  costheta_el = G4UniformRand()*2.0-1.0;
  G4double kk = std::sqrt(electronKineEnergy*(electronKineEnergy+2.*electron_mass_c2));
  costheta_el = (costheta_el*electronTotEnergy+kk)/(electronTotEnergy+costheta_el*kk);
  phi_el  = twopi * G4UniformRand() ;
  G4double dirX_el = std::sqrt(1.-costheta_el*costheta_el) * std::cos(phi_el);
  G4double dirY_el = std::sqrt(1.-costheta_el*costheta_el) * std::sin(phi_el);
  G4double dirZ_el = costheta_el;

  //positron kinematics
  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;
  costheta_po = G4UniformRand()*2.0-1.0;
  kk = std::sqrt(positronKineEnergy*(positronKineEnergy+2.*electron_mass_c2));
  costheta_po = (costheta_po*positronTotEnergy+kk)/(positronTotEnergy+costheta_po*kk);
  phi_po  = twopi * G4UniformRand() ;
  G4double dirX_po = std::sqrt(1.-costheta_po*costheta_po) * std::cos(phi_po);
  G4double dirY_po = std::sqrt(1.-costheta_po*costheta_po) * std::sin(phi_po);
  G4double dirZ_po = costheta_po;

  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular
  // distribution with respect to the Z axis along the parent photon
  G4double localEnergyDeposit = 0. ;

  //Generate explicitely the electron in the pair, only if it is > threshold
  //VI: applying cut here provides inconsistency 

  if (electronKineEnergy > 0.0)
    {
      G4ThreeVector electronDirection ( dirX_el, dirY_el, dirZ_el);
      electronDirection.rotateUz(photonDirection);
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   electronDirection,
							   electronKineEnergy);
      fvect->push_back(electron);
    }
  else
    {
      localEnergyDeposit += electronKineEnergy;
      electronKineEnergy = 0;
    }

  //Generate the positron. Real particle in any case, because it will annihilate. If below
  //threshold, produce it at rest
  // VI: here there was a bug - positron and electron cuts are different
  if (positronKineEnergy < 0.0)
    {
      localEnergyDeposit += positronKineEnergy;
      positronKineEnergy = 0; //produce it at rest
    }
  G4ThreeVector positronDirection(dirX_po,dirY_po,dirZ_po);
  positronDirection.rotateUz(photonDirection);
  G4DynamicParticle* positron = new G4DynamicParticle(G4Positron::Positron(),
						      positronDirection, positronKineEnergy);
  fvect->push_back(positron);

  //Add rest of energy to the local energy deposit
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);
  
  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeGammaConversion" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      if (electronKineEnergy)
	G4cout << "Electron (explicitely produced) " << electronKineEnergy/keV << " keV" 
	       << G4endl;
      if (positronKineEnergy)
	G4cout << "Positron (not at rest) " << positronKineEnergy/keV << " keV" << G4endl;
      G4cout << "Rest masses of e+/- " << 2.0*electron_mass_c2/keV << " keV" << G4endl;
      if (localEnergyDeposit)
	G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (electronKineEnergy+positronKineEnergy+
					  localEnergyDeposit+2.0*electron_mass_c2)/keV <<
        " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
 if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(electronKineEnergy+positronKineEnergy+
				      localEnergyDeposit+2.0*electron_mass_c2-photonEnergy);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeGammaConversion: problem with energy conservation: " 
	       << (electronKineEnergy+positronKineEnergy+
		   localEnergyDeposit+2.0*electron_mass_c2)/keV 
	       << " keV (final) vs. " << photonEnergy/keV << " keV (initial)" << G4endl;
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4double> G4PenelopeGammaConversionModel::ScreenFunction(G4double b)
{
  std::vector<G4double> result;
  result.clear();
  G4double bsquare=b*b;
  G4double a0,f1,f2;
  f1=2.0-2*std::log(1+bsquare);
  f2=f1-(2.0/3.0);
  if (b < 1.0e-10)
    {
      f1=f1-twopi*b;
    }
  else
    {
      a0 = 4*b*std::atan(1.0/b);
      f1 = f1 - a0;
      f2 = f2+2*bsquare*(4.0-a0-3*std::log((1+bsquare)/bsquare));
    }
  result.push_back(0.5*(3*f1-f2));
  result.push_back(0.25*(3*f1+f2));
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeGammaConversionModel::GetScreeningRadius(G4double Z)
{
  G4double result = 0;
  G4bool foundElement = false;
  G4int iZ = (G4int) Z;
  if (!fTheScreeningRadii)
    fTheScreeningRadii = new std::map<G4int,G4double>;
  
  if (fTheScreeningRadii->count(iZ))
    {
      //The element is already loaded: just return it
      result = fTheScreeningRadii->find(iZ)->second;
      return result;
    }
  else //retrieve all from file
    {
      char* path = getenv("G4LEDATA");
      if (!path)
	{
	  G4String excep = "G4PenelopeGammaConversionModel - G4LEDATA environment variable not set!";
	  G4Exception(excep);
	  return result;
	}
      G4String pathString(path);
      G4String pathFile = pathString + "/penelope/pp-pen.dat";
      std::ifstream file(pathFile);
      
      if (!(file.is_open()))
	{
	  G4String excep = "G4PenelopeGammaConversionModel - data file " + pathFile + "not found!";
	  G4Exception(excep);
	}
      G4int k;
      G4double a1,a2;
      while(!file.eof()) {
	file >> k >> a1 >> a2;
	fTheScreeningRadii->insert(std::make_pair(k,a1));
	if ((G4double) k == Z)
	  {
	    result = a1;
	    foundElement = true;
	  }
      }
      file.close();
      if (verboseLevel > 2)
	G4cout << "Read file pp-pen.dat" << G4endl;
      if (foundElement)
	return result;
      else
	{
	  G4String excep = "G4PenelopeGammaConversionModel - Screening Radius for not found in the data file";
	  G4Exception(excep);
	  return 0;
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeGammaConversionModel::CoulombCorrection(G4double a)
{
  G4double fc=0;
  G4double b[7] = {0.202059,-0.03693,0.00835,-0.00201,0.00049,-0.00012,0.00003};
  G4double aSquared = a*a;
  G4double aFourth = aSquared*aSquared;
  G4double aEighth = aFourth*aFourth;

  fc = ((1.0/(1.0+a*a))+b[0]+b[1]*aSquared+b[2]*aFourth+b[3]*(aSquared*aFourth)+
	b[4]*aEighth+b[5]*(aEighth*aSquared)+b[6]*(aEighth*aFourth));
  fc=aSquared*fc;
  return fc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeGammaConversionModel::LowEnergyCorrection(G4double a,G4double eki)
{
  G4double f0=0,t=0;
  G4double b[12] = {-1.744,-12.10,11.18,8.523,73.26,-41.41,-13.52,-121.1,94.41,8.946,62.05,-63.41};
  t=std::sqrt(2.0*eki);
  G4double tSq = t*t;
  f0=(b[0]+b[1]*a+b[2]*a*a)*t+(b[3]+b[4]*a+b[5]*a*a)*(tSq)+(b[6]+b[7]*a+b[8]*a*a)*(tSq*t)+
    (b[9]+b[10]*a+b[11]*a*a)*(tSq*tSq);
  return f0;

}
