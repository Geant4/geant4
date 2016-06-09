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
// --------------------------------------------------------------------
//
//
//
// --------------------------------------------------------------
//
// Author: L.Pandola
// History:
// --------
// 02 Dec 2002 L.Pandola    1st implementation
// 12 Feb 2003   MG Pia     Migration to "cuts per region"
// 10 Mar 2003 V.Ivanchenko Remove CutPerMaterial warning
// 13 Mar 2003 L.Pandola    Code "cleaned"
// 25 Mar 2003 L.Pandola    Changed the name of the database file to read
// 24 Apr 2003 V.Ivanchenko Cut per region mfpt
// 17 Mar 2004 L.Pandola    Removed unnecessary calls to std::pow(a,b)
// --------------------------------------------------------------

#include "G4PenelopeGammaConversion.hh"

#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4Positron.hh"
#include "G4IonisParamElm.hh"
#include "G4Material.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VRangeTest.hh"
#include "G4RangeTest.hh"
#include "G4MaterialCutsCouple.hh"


G4PenelopeGammaConversion::G4PenelopeGammaConversion(const G4String& processName)
  : G4VDiscreteProcess(processName),   
    lowEnergyLimit(1.022000*MeV),
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(1.022000*MeV),
    intrinsicHighEnergyLimit(100*GeV),
    smallEnergy(1.1*MeV)

{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4PenelopeGammaConversion::G4PenelopeGammaConversion - energy limit outside intrinsic process validity range");
    }

  // The following pointer is owned by G4DataHandler  
  crossSectionHandler = new G4CrossSectionHandler();
  // Log log interpolation (default)
  crossSectionHandler->Initialise(0,1.0220*MeV,100.*GeV,400);
  meanFreePathTable = 0;
  rangeTest = new G4RangeTest;

   if (verboseLevel > 0) 
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: " 
	      << lowEnergyLimit / MeV << " MeV - "
	      << highEnergyLimit / GeV << " GeV" 
	      << G4endl;
     }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4PenelopeGammaConversion is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}
 
G4PenelopeGammaConversion::~G4PenelopeGammaConversion()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete rangeTest;
}

void G4PenelopeGammaConversion::BuildPhysicsTable(const G4ParticleDefinition& )
{

  crossSectionHandler->Clear();
  G4String crossSectionFile = "penelope/pp-cs-pen-";
  crossSectionHandler->LoadData(crossSectionFile);
  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4PenelopeGammaConversion::PostStepDoIt(const G4Track& aTrack,
							    const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
 
  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy = incidentPhoton->GetKineticEnergy();
  G4ParticleMomentum photonDirection = incidentPhoton->GetMomentumDirection();

  G4double eps ;
  G4double eki = electron_mass_c2 / photonEnergy ;

  // Do it fast if photon energy < 1.1 MeV
  if (photonEnergy < smallEnergy )
    {
      eps = eki + (1-2*eki) * G4UniformRand();
    }
  else
    {
      // Select randomly one element in the current material
      const G4Element* element = crossSectionHandler->SelectRandomElement(couple,photonEnergy);

      if (element == 0)
	{
	  G4cout << "G4PenelopeGammaConversion::PostStepDoIt - element = 0" << G4endl;
	}
      G4IonisParamElm* ionisation = element->GetIonisation();
      if (ionisation == 0) 
	{
	  G4cout << "G4PenelopeGammaConversion::PostStepDoIt - ionisation = 0" << G4endl;
	}
      
      //Low energy and Coulomb corrections
      G4double Z=ionisation->GetZ(); 
      G4double ZAlpha = Z*fine_structure_const;
      G4double ScreenRadius = GetScreeningRadius(Z);
      G4double funct1=0,g0=0;
      G4double g1min=0,g2min=0;
      funct1 = 4.0*std::log(ScreenRadius);
      g0 = funct1-4*CoulombCorrection(ZAlpha)+LowEnergyCorrection(ZAlpha,eki); 
      G4double bmin = 2*eki*ScreenRadius;
      g1min=g0+ScreenFunction(bmin,1);
      g2min=g0+ScreenFunction(bmin,2);
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
	  g1 = g0+ScreenFunction(b,1);
	  if (g1 < 0) g1=0;
	  rand3 = G4UniformRand()*g1min;
        }
        else
	  {
	    eps = eki+2.0*xr*G4UniformRand();
	    b = (eki*ScreenRadius)/(2*eps*(1.0-eps));
	    g1 = g0+ScreenFunction(b,2); 
	    if (g1 < 0) g1=0; 
	    rand3 = G4UniformRand()*g2min;
	  }	   
      } while (rand3>g1);
    } //End of eps sampling

  G4double electronTotEnergy;
  G4double positronTotEnergy;

  electronTotEnergy = eps*photonEnergy;
  positronTotEnergy = (1.0-eps)*photonEnergy;
  
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

  aParticleChange.SetNumberOfSecondaries(2) ;


  // Generate the electron only if with large enough range w.r.t. cuts and safety

  G4double safety = aStep.GetPostStepPoint()->GetSafety();

  if (rangeTest->Escape(G4Electron::Electron(),couple,electronKineEnergy,safety))
    {
      G4ThreeVector electronDirection ( dirX_el, dirY_el, dirZ_el);
      electronDirection.rotateUz(photonDirection);
      G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
							    electronDirection,
							    electronKineEnergy);
      aParticleChange.AddSecondary(particle1) ;
    }
  else
    {
      localEnergyDeposit += electronKineEnergy ;
    }


  if (! (rangeTest->Escape(G4Positron::Positron(),couple,positronKineEnergy,safety)))
    {
      localEnergyDeposit += positronKineEnergy ;
      positronKineEnergy = 0. ;
    }
  G4ThreeVector positronDirection(dirX_po,dirY_po,dirZ_po);
  positronDirection.rotateUz(photonDirection);

  // Create G4DynamicParticle object for the particle2
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
						       positronDirection, positronKineEnergy);
  aParticleChange.AddSecondary(particle2) ;

  aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposit) ;

  // Kill the incident photon
  aParticleChange.ProposeMomentumDirection(0.,0.,0.) ;
  aParticleChange.ProposeEnergy(0.) ;
  aParticleChange.ProposeTrackStatus(fStopAndKill) ;

  //  Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
}

G4bool G4PenelopeGammaConversion::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() );
}

G4double G4PenelopeGammaConversion::GetMeanFreePath(const G4Track& track,
						    G4double, // previousStepSize
						    G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  size_t materialIndex = couple->GetIndex();

  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
  return meanFreePath;
}

G4double G4PenelopeGammaConversion::ScreenFunction(G4double b,G4int icase)
{
  G4double bsquare=b*b;
  G4double a0,f1,f2,g1,g2;
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
  g1=0.5*(3*f1-f2);
  g2=0.25*(3*f1+f2);
  if (icase==1) {
   return g1;
  }
  else
    {
      return g2;
    }
}
      
G4double G4PenelopeGammaConversion::CoulombCorrection(G4double a)
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

G4double G4PenelopeGammaConversion::LowEnergyCorrection(G4double a,G4double eki)
{
  G4double f0=0,t=0;
  G4double b[12] = {-1.744,-12.10,11.18,8.523,73.26,-41.41,-13.52,-121.1,94.41,8.946,62.05,-63.41};
  t=std::sqrt(2.0*eki);
  G4double tSq = t*t;
  f0=(b[0]+b[1]*a+b[2]*a*a)*t+(b[3]+b[4]*a+b[5]*a*a)*(tSq)+(b[6]+b[7]*a+b[8]*a*a)*(tSq*t)+
    (b[9]+b[10]*a+b[11]*a*a)*(tSq*tSq);
  return f0;
}

G4double G4PenelopeGammaConversion::GetScreeningRadius(G4double Z) 
{
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeGammaConversion - G4LEDATA environment variable not set!";
      G4Exception(excep);
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/pp-pen.dat";
  std::ifstream file(pathFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (!(lsdp->is_open()))
    {
      G4String excep = "G4PenelopeGammaConversion - data file " + pathFile + "not found!";
      G4Exception(excep);
    }
  G4int k;
  G4double a1,a2;
  while(!file.eof()) {
    file >> k >> a1 >> a2;
    if ((G4double) k == Z)
      {
	return a1;
      }
  } 
  G4String excep = "G4PenelopeGammaConversion - Screening Radius for not found in the data file";
  G4Exception(excep);
  return 0;
}
