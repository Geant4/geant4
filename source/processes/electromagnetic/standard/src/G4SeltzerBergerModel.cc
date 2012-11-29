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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SeltzerBergerModel
//
// Author:        Vladimir Ivanchenko use inheritance from Andreas Schaelicke
//                base class implementing ultra relativistic bremsstrahlung
//                model 
//
// Creation date: 04.10.2011
//
// Modifications:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4SeltzerBergerModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4ModifiedTsai.hh"

#include "G4Physics2DVector.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4Physics2DVector* G4SeltzerBergerModel::dataSB[101] = {0};
G4double G4SeltzerBergerModel::ylimit[101] = {0.0};
G4double G4SeltzerBergerModel::expnumlim = -12.;

G4SeltzerBergerModel::G4SeltzerBergerModel(const G4ParticleDefinition* p,
					   const G4String& nam)
  : G4eBremsstrahlungRelModel(p,nam),useBicubicInterpolation(false)
{
  SetLowEnergyLimit(0.0);
  SetLPMFlag(false);
  nwarn = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SeltzerBergerModel::~G4SeltzerBergerModel()
{
  for(size_t i=0; i<101; ++i) { 
    if(dataSB[i]) {
      delete dataSB[i]; 
      dataSB[i] = 0;
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SeltzerBergerModel::Initialise(const G4ParticleDefinition* p,
				      const G4DataVector& cuts)
{
  // check environment variable
  // Build the complete string identifying the file with the data set
  char* path = getenv("G4LEDATA");

  // Access to elements
  const G4ElementTable* theElmTable = G4Element::GetElementTable();
  size_t numOfElm = G4Element::GetNumberOfElements();
  if(numOfElm > 0) {
    for(size_t i=0; i<numOfElm; ++i) {
      G4int Z = G4int(((*theElmTable)[i])->GetZ());
      if(Z < 1)        { Z = 1; }
      else if(Z > 100) { Z = 100; }
      //G4cout << "Z= " << Z << G4endl;
      // Initialisation
      if(!dataSB[Z]) { ReadData(Z, path); }
    }
  }

  G4eBremsstrahlungRelModel::Initialise(p, cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SeltzerBergerModel::ReadData(size_t Z, const char* path)
{
  //  G4cout << "ReadData Z= " << Z << G4endl;
  // G4cout << "Status for Z= " << dataSB[Z] << G4endl;
  //if(path) { G4cout << path << G4endl; }
  if(dataSB[Z]) { return; }
  const char* datadir = path;

  if(!datadir) {
    datadir = getenv("G4LEDATA");
    if(!datadir) {
      G4Exception("G4SeltzerBergerModel::ReadData()","em0006",FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }
  std::ostringstream ost;
  ost << datadir << "/brem_SB/br" << Z;
  std::ifstream fin(ost.str().c_str());
  if( !fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "Bremsstrahlung data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4SeltzerBergerModel::ReadData()","em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.23 or later.");
    return;
  } 
  //G4cout << "G4SeltzerBergerModel read from <" << ost.str().c_str() 
  //	 << ">" << G4endl;
  G4Physics2DVector* v = new G4Physics2DVector();
  const G4double emaxlog = 4*log(10.);
  if(v->Retrieve(fin)) { 
    if(useBicubicInterpolation) { v->SetBicubicInterpolation(true); }
    dataSB[Z] = v; 
    ylimit[Z] = v->Value(0.97, emaxlog);
  } else {
    G4ExceptionDescription ed;
    ed << "Bremsstrahlung data file <" << ost.str().c_str()
       << "> is not retrieved!" << G4endl;
    G4Exception("G4SeltzerBergerModel::ReadData()","em0005",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.23 or later.");
    delete v;
  }
  // G4cout << dataSB[Z] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4SeltzerBergerModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
{

  if(gammaEnergy < 0.0 || kinEnergy <= 0.0) { return 0.0; }
  G4double x = gammaEnergy/kinEnergy;
  G4double y = log(kinEnergy/MeV);
  G4int Z = G4int(currentZ);
  //G4cout << "G4SeltzerBergerModel::ComputeDXSectionPerAtom Z= " << Z
  //	 << " x= " << x << " y= " << y << " " << dataSB[Z] << G4endl;
  if(!dataSB[Z]) { ReadData(Z); }
  G4double invb2 = totalEnergy*totalEnergy/(kinEnergy*(kinEnergy + 2*particleMass));
  G4double cross = dataSB[Z]->Value(x,y)*invb2*millibarn/bremFactor;
  
  if(!isElectron) {
    G4double invbeta1 = sqrt(invb2);
    G4double e2 = kinEnergy - gammaEnergy;
    if(e2 > 0.0) {
      G4double invbeta2 = (e2 + particleMass)/sqrt(e2*(e2 + 2*particleMass));
      G4double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2);
      if(xxx < expnumlim) { cross = 0.0; }
      else { cross *= exp(xxx); }
    } else {
      cross = 0.0;
    }
  }
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4SeltzerBergerModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp, 
					const G4MaterialCutsCouple* couple,
					const G4DynamicParticle* dp,
					G4double cutEnergy,
					G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  G4double emax = std::min(maxEnergy, kineticEnergy);
  if(cut >= emax) { return; }

  SetupForMaterial(particle, couple->GetMaterial(), kineticEnergy);

  const G4Element* elm = 
    SelectRandomAtom(couple,particle,kineticEnergy,cut,emax);
  SetCurrentElement(elm->GetZ());
  G4int Z = G4int(currentZ);

  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  /*
  G4cout << "G4SeltzerBergerModel::SampleSecondaries E(MeV)= " 
	 << kineticEnergy/MeV
	 << " Z= " << Z << " cut(MeV)= " << cut/MeV 
	 << " emax(MeV)= " << emax/MeV << " corr= " << densityCorr << G4endl;
  */
  G4double xmin = log(cut*cut + densityCorr);
  G4double xmax = log(emax*emax  + densityCorr);
  G4double y = log(kineticEnergy/MeV);

  G4double gammaEnergy, v; 

  // majoranta
  G4double x0 = cut/kineticEnergy;
  G4double vmax = dataSB[Z]->Value(x0, y)*1.02;
  //  G4double invbeta1 = 0;

  const G4double epeaklimit= 300*MeV; 
  const G4double elowlimit = 10*keV; 

  // majoranta corrected for e-
  if(isElectron && x0 < 0.97 && 
     ((kineticEnergy > epeaklimit) || (kineticEnergy < elowlimit))) {
    G4double ylim = std::min(ylimit[Z],1.1*dataSB[Z]->Value(0.97, y)); 
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }

  //G4cout<<"y= "<<y<<" xmin= "<<xmin<<" xmax= "<<xmax<<" vmax= "<<vmax<<G4endl;
  //  G4int ncount = 0;
  do {
    //++ncount;
    G4double x = exp(xmin + G4UniformRand()*(xmax - xmin)) - densityCorr;
    if(x < 0.0) { x = 0.0; }
    gammaEnergy = sqrt(x);
    G4double x1 = gammaEnergy/kineticEnergy;
    v = dataSB[Z]->Value(x1, y);

    // correction for positrons        
    if(!isElectron) {
      G4double e1 = kineticEnergy - cut;
      G4double invbeta1 = (e1 + particleMass)/sqrt(e1*(e1 + 2*particleMass));
      G4double e2 = kineticEnergy - gammaEnergy;
      G4double invbeta2 = (e2 + particleMass)/sqrt(e2*(e2 + 2*particleMass));
      G4double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2); 

      if(xxx < expnumlim) { v = 0.0; }
      else { v *= exp(xxx); }
    }
   
    if (v > 1.05*vmax && nwarn < 20) {
      ++nwarn;
      G4cout << "### G4SeltzerBergerModel Warning: Majoranta exceeded! "
	     << v << " > " << vmax << " by " << v/vmax
	     << " Egamma(MeV)= " << gammaEnergy
	     << " Ee(MeV)= " << kineticEnergy
	     << " Z= " << Z << "  " << particle->GetParticleName()
	//<< " ncount= " << ncount
	     << G4endl;
     
      if ( 20 == nwarn ) {
	G4cout << "### G4SeltzerBergerModel Warnings will not be printed anymore"
	       << G4endl;
      }
    }
  } while (v < vmax*G4UniformRand());

  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //

  G4ThreeVector gammaDirection = 
    GetAngularDistribution()->SampleDirection(dp, totalEnergy-gammaEnergy,
					      Z, couple->GetMaterial());

  // create G4DynamicParticle object for the Gamma
  G4DynamicParticle* gamma = 
    new G4DynamicParticle(theGamma,gammaDirection,gammaEnergy);
  vdp->push_back(gamma);
  
  G4ThreeVector direction = (totMomentum*dp->GetMomentumDirection()
			     - gammaEnergy*gammaDirection).unit();

  /*
  G4cout << "### G4SBModel: v= "
	 << " Eg(MeV)= " << gammaEnergy
	 << " Ee(MeV)= " << kineticEnergy
	 << " DirE " << direction << " DirG " << gammaDirection
	 << G4endl;
  */
  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    G4DynamicParticle* el = 
      new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
			    direction, finalE);
    vdp->push_back(el);

    // continue tracking
  } else {
    fParticleChange->SetProposedMomentumDirection(direction);
    fParticleChange->SetProposedKineticEnergy(finalE);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


