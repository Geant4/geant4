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
// Authors: G.Depaola & F.Longo
//
//

#include "G4LivermorePolarizedGammaConversionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4AutoLock.hh"
#include "G4Exp.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;
namespace { G4Mutex LivermorePolarizedGammaConversionModelMutex = G4MUTEX_INITIALIZER; }

G4PhysicsFreeVector* G4LivermorePolarizedGammaConversionModel::data[] = {nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedGammaConversionModel::G4LivermorePolarizedGammaConversionModel(
   const G4ParticleDefinition*, const G4String& nam)
  :G4VEmModel(nam), smallEnergy(2.*MeV), isInitialised(false)
{
  fParticleChange = nullptr;
  lowEnergyLimit = 2*electron_mass_c2;
  
  Phi=0.;
  Psi=0.;
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = calculation of cross sections, file openings, samping of atoms
  // 2 = entering in methods
  
  if(verboseLevel > 0) {
    G4cout << "Livermore Polarized GammaConversion is constructed " 
	   << G4endl;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedGammaConversionModel::~G4LivermorePolarizedGammaConversionModel()
{  
  if(IsMaster()) {
    for(G4int i=0; i<maxZ; ++i) {
      if(data[i]) { 
	delete data[i];
	data[i] = nullptr;
      }
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedGammaConversionModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 1)
    {
      G4cout << "Calling1 G4LivermorePolarizedGammaConversionModel::Initialise()" 
	     << G4endl
             << "Energy range: "
	     << LowEnergyLimit() / MeV << " MeV - "
             << HighEnergyLimit() / GeV << " GeV"
             << G4endl;
    }
  
  if(IsMaster()) 
    {
      // Initialise element selector    
      InitialiseElementSelectors(particle, cuts);
      
      // Access to elements
      const char* path = G4FindDataDir("G4LEDATA");
      
      G4ProductionCutsTable* theCoupleTable =
	G4ProductionCutsTable::GetProductionCutsTable();
      
      G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
      
      for(G4int i=0; i<numOfCouples; ++i) 
	{
	  const G4Material* material = 
	    theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
	  const G4ElementVector* theElementVector = material->GetElementVector();
	  std::size_t nelm = material->GetNumberOfElements();
	  
	  for (std::size_t j=0; j<nelm; ++j) 
	    {
	      G4int Z = (G4int)(*theElementVector)[j]->GetZ();
	      if(Z < 1)          { Z = 1; }
	      else if(Z > maxZ)  { Z = maxZ; }
	      if(!data[Z]) { ReadData(Z, path); }
	    }
	}
    }
  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedGammaConversionModel::InitialiseLocal(
     const G4ParticleDefinition*, G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::MinPrimaryEnergy(const G4Material*,
			     const G4ParticleDefinition*, G4double)
{
  return lowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedGammaConversionModel::ReadData(std::size_t Z, const char* path)
{
  if (verboseLevel > 1) 
    {
      G4cout << "Calling ReadData() of G4LivermorePolarizedGammaConversionModel" 
	     << G4endl;
    }
  
  if(data[Z]) { return; }
  
  const char* datadir = path;
  
  if(!datadir) 
    {
      datadir = G4FindDataDir("G4LEDATA");
      if(!datadir) 
	{
	  G4Exception("G4LivermorePolarizedGammaConversionModel::ReadData()",
		      "em0006",FatalException,
		      "Environment variable G4LEDATA not defined");
	  return;
	}
    }
  //  
  data[Z] = new G4PhysicsFreeVector(0,/*spline=*/true);
  //
  std::ostringstream ost;
  ost << datadir << "/livermore/pair/pp-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  
  if( !fin.is_open()) 
    {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedGammaConversionModel data file <" << ost.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4LivermorePolarizedGammaConversionModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.27 or later.");
      return;
    } 
  else 
    {
      
      if(verboseLevel > 3) { G4cout << "File " << ost.str() 
				    << " is opened by G4LivermorePolarizedGammaConversionModel" << G4endl;}
      
      data[Z]->Retrieve(fin, true);
    } 
  
  // Activation of spline interpolation
  data[Z]->FillSecondDerivatives();  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
				       G4double GammaEnergy,
				       G4double Z, G4double,
				       G4double, G4double)
{
  if (verboseLevel > 1) {
    G4cout << "G4LivermorePolarizedGammaConversionModel::ComputeCrossSectionPerAtom()" 
	   << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit) { return 0.0; } 
  
  G4double xs = 0.0;
  
  G4int intZ=G4int(Z);
  
  if(intZ < 1 || intZ > maxZ) { return xs; }
  
  G4PhysicsFreeVector* pv = data[intZ];
  
  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) 
    {
      InitialiseForElement(0, intZ);
      pv = data[intZ];
      if(!pv) { return xs; }
    }
  // x-section is taken from the table
  xs = pv->Value(GammaEnergy); 
  
  if(verboseLevel > 0)
    {
      G4int n = G4int(pv->GetVectorLength() - 1);
      G4cout  <<  "****** DEBUG: tcs value for Z=" << Z << " at energy (MeV)=" 
	      << GammaEnergy/MeV << G4endl;
      G4cout  <<  "  cs (Geant4 internal unit)=" << xs << G4endl;
      G4cout  <<  "    -> first cs value in EADL data file (iu) =" << (*pv)[0] << G4endl;
      G4cout  <<  "    -> last  cs value in EADL data file (iu) =" << (*pv)[n] << G4endl;
      G4cout  <<  "*********************************************************" << G4endl;
    }
  
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedGammaConversionModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
							    const G4MaterialCutsCouple* couple,
							    const G4DynamicParticle* aDynamicGamma,
							    G4double,
							    G4double)
{

  // Fluorescence generated according to:
  // J. Stepanek ,"A program to determine the radiation spectra due to a single atomic
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides",
  // Comp. Phys. Comm. 1206 pp 1-1-9 (1997)
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedGammaConversionModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();

  if(photonEnergy <= lowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      return;
    }

  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization();
  G4ThreeVector gammaDirection0 = aDynamicGamma->GetMomentumDirection();

  // Make sure that the polarization vector is perpendicular to the
  // gamma direction. If not
  if(!(gammaPolarization0.isOrthogonal(gammaDirection0, 1e-6))||(gammaPolarization0.mag()==0))
    { // only for testing now
      gammaPolarization0 = GetRandomPolarization(gammaDirection0);
    }
  else
    {
      if ( gammaPolarization0.howOrthogonal(gammaDirection0) != 0)
	{
	  gammaPolarization0 = GetPerpendicularPolarization(gammaDirection0, gammaPolarization0);
	}
    }

  // End of Protection

  G4double epsilon ;
  G4double epsilon0Local = electron_mass_c2 / photonEnergy ;

  // Do it fast if photon energy < 2. MeV

  if (photonEnergy < smallEnergy )
    {
      epsilon = epsilon0Local + (0.5 - epsilon0Local) * G4UniformRand();
    }
  else
    {
      // Select randomly one element in the current material 
      const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
      const G4Element* element = SelectRandomAtom(couple,particle,photonEnergy);
      
      if (element == nullptr)
        {
          G4cout << "G4LivermorePolarizedGammaConversionModel::SampleSecondaries - element = 0" << G4endl;
	  return;
        }
      
      
      G4IonisParamElm* ionisation = element->GetIonisation();      
      if (ionisation == nullptr)
        {
          G4cout << "G4LivermorePolarizedGammaConversionModel::SampleSecondaries - ionisation = 0" << G4endl;
	  return;
        }
      
      // Extract Coulomb factor for this Element
      G4double fZ = 8. * (ionisation->GetlogZ3());
      if (photonEnergy > 50. * MeV) fZ += 8. * (element->GetfCoulomb());

      // Limits of the screening variable
      G4double screenFactor = 136. * epsilon0Local / (element->GetIonisation()->GetZ3()) ;
      G4double screenMax = G4Exp ((42.24 - fZ)/8.368) - 0.952 ;
      G4double screenMin = std::min(4.*screenFactor,screenMax) ;

      // Limits of the energy sampling
      G4double epsilon1 = 0.5 - 0.5 * sqrt(1. - screenMin / screenMax) ;
      G4double epsilonMin = std::max(epsilon0Local,epsilon1);
      G4double epsilonRange = 0.5 - epsilonMin ;

      // Sample the energy rate of the created electron (or positron)
      G4double screen;
      G4double gReject ;

      G4double f10 = ScreenFunction1(screenMin) - fZ;
      G4double f20 = ScreenFunction2(screenMin) - fZ;
      G4double normF1 = std::max(f10 * epsilonRange * epsilonRange,0.);
      G4double normF2 = std::max(1.5 * f20,0.);

      do {
        if (normF1 / (normF1 + normF2) > G4UniformRand() )
          {
            epsilon = 0.5 - epsilonRange * pow(G4UniformRand(), 0.3333) ;
            screen = screenFactor / (epsilon * (1. - epsilon));
            gReject = (ScreenFunction1(screen) - fZ) / f10 ;
          }
        else
          {
            epsilon = epsilonMin + epsilonRange * G4UniformRand();
            screen = screenFactor / (epsilon * (1 - epsilon));
            gReject = (ScreenFunction2(screen) - fZ) / f20 ;
	  }
      } while ( gReject < G4UniformRand() );
    }   //  End of epsilon sampling
  
  // Fix charges randomly
  G4double electronTotEnergy;
  G4double positronTotEnergy;

  if (G4UniformRand() > 0.5)  
    {
      electronTotEnergy = (1. - epsilon) * photonEnergy;
      positronTotEnergy = epsilon * photonEnergy;
    }
  else
    {
      positronTotEnergy = (1. - epsilon) * photonEnergy;
      electronTotEnergy = epsilon * photonEnergy;
    }
  
  // Scattered electron (positron) angles. ( Z - axis along the parent photon)
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  // derived from Tsai distribution (Rev. Mod. Phys. 49, 421 (1977)
  G4double Ene = electronTotEnergy/electron_mass_c2; // Normalized energy

  G4double cosTheta = 0.;
  G4double sinTheta = 0.;

  SetTheta(&cosTheta,&sinTheta,Ene);
  G4double phi,psi=0.;

  //corrected e+ e- angular angular distribution //preliminary!
  phi = SetPhi(photonEnergy);
  psi = SetPsi(photonEnergy,phi);
  Psi = psi;
  Phi = phi;

  G4double phie, phip; 
  G4double choice, choice2;
  choice = G4UniformRand();
  choice2 = G4UniformRand();

  if (choice2 <= 0.5)
    {
      // do nothing 
      //  phi = phi;
    }
  else
    {
      phi = -phi;
    }
  
  if (choice <= 0.5)
    {
      phie = psi; //azimuthal angle for the electron
      phip = phie+phi; //azimuthal angle for the positron
    }
  else
    {
      // opzione 1 phie / phip equivalenti
      phip = psi; //azimuthal angle for the positron
      phie = phip + phi; //azimuthal angle for the electron
    }


  // Electron Kinematics 
  G4double dirX = sinTheta*cos(phie);
  G4double dirY = sinTheta*sin(phie);
  G4double dirZ = cosTheta;
  G4ThreeVector electronDirection(dirX,dirY,dirZ);

  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular
  // distribution with respect to the Z axis along the parent photon

  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ;

  SystemOfRefChange(gammaDirection0,electronDirection,
		    gammaPolarization0);

  G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
							electronDirection,
							electronKineEnergy);

  // The e+ is always created (even with kinetic energy = 0) for further annihilation
  Ene = positronTotEnergy/electron_mass_c2; // Normalized energy

  cosTheta = 0.;
  sinTheta = 0.;

  SetTheta(&cosTheta,&sinTheta,Ene);

  // Positron Kinematics
  dirX = sinTheta*cos(phip);
  dirY = sinTheta*sin(phip);
  dirZ = cosTheta;
  G4ThreeVector positronDirection(dirX,dirY,dirZ);

  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;
  SystemOfRefChange(gammaDirection0,positronDirection,
		    gammaPolarization0);

  // Create G4DynamicParticle object for the particle2
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
                                                       positronDirection, positronKineEnergy);
  fvect->push_back(particle1);
  fvect->push_back(particle2);

  // Kill the incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::ScreenFunction1(G4double screenVariable)
{
  // Compute the value of the screening function 3*phi1 - phi2
  G4double value;
  if (screenVariable > 1.)
    value = 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    value = 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);

  return value;
}



G4double G4LivermorePolarizedGammaConversionModel::ScreenFunction2(G4double screenVariable)
{
  // Compute the value of the screening function 1.5*phi1 - 0.5*phi2
  G4double value;

  if (screenVariable > 1.)
    value = 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    value = 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);

  return value;
}


void G4LivermorePolarizedGammaConversionModel::SetTheta(G4double* p_cosTheta, G4double* p_sinTheta, G4double Energy)
{
  // to avoid computational errors since Theta could be very small
  // Energy in Normalized Units (!)

  G4double Momentum = sqrt(Energy*Energy -1);
  G4double Rand = G4UniformRand();

  *p_cosTheta = (Energy*((2*Rand)- 1) + Momentum)/((Momentum*(2*Rand-1))+Energy);
  *p_sinTheta = (2*sqrt(Rand*(1-Rand)))/(Momentum*(2*Rand-1)+Energy);
}



G4double G4LivermorePolarizedGammaConversionModel::SetPhi(G4double Energy)
{
  G4double value = 0.;
  G4double Ene = Energy/MeV;

  G4double pl[4];
  G4double pt[2];
  G4double xi = 0;
  G4double xe = 0.;
  G4double n1=0.;
  G4double n2=0.;

  if (Ene>=50.)
    {
      const G4double ay0=5.6, by0=18.6, aa0=2.9, ba0 = 8.16E-3;
      const G4double aw = 0.0151, bw = 10.7, cw = -410.;

      const G4double axc = 3.1455, bxc = -1.11, cxc = 310.;

      pl[0] = Fln(ay0,by0,Ene);
      pl[1] = aa0 + ba0*(Ene);
      pl[2] = Poli(aw,bw,cw,Ene);
      pl[3] = Poli(axc,bxc,cxc,Ene);

      const G4double abf = 3.1216, bbf = 2.68;
      pt[0] = -1.4;
      pt[1] = abf + bbf/Ene;

      xi = 3.0;
      xe = Encu(pl,pt,xi);
      n1 = Fintlor(pl,pi) - Fintlor(pl,xe);
      n2 = Finttan(pt,xe) - Finttan(pt,0.);
    }
  else
    {
      const G4double ay0=0.144, by0=0.11;
      const G4double aa0=2.7, ba0 = 2.74;
      const G4double aw = 0.21, bw = 10.8, cw = -58.;
      const G4double axc = 3.17, bxc = -0.87, cxc = -6.;

      pl[0] = Fln(ay0, by0, Ene);
      pl[1] = Fln(aa0, ba0, Ene);
      pl[2] = Poli(aw,bw,cw,Ene);
      pl[3] = Poli(axc,bxc,cxc,Ene);

      n1 = Fintlor(pl,pi) - Fintlor(pl,xe);
    }


  G4double n=0.;
  n = n1+n2;

  G4double c1 = 0.;
  c1 = Glor(pl, xe);

  G4double r1,r2,r3;
  G4double xco=0.;

  if (Ene>=50.)
    {
      r1= G4UniformRand();
      if( r1>=n2/n)
        {
          do
	    {
              r2 = G4UniformRand();
              value = Finvlor(pl,xe,r2);
              xco = Glor(pl,value)/c1;
              r3 = G4UniformRand();
            } while(r3>=xco);
        }
      else
        {
          value = Finvtan(pt,n,r1);
        }
    }
  else
    {
      do
        {
          r2 = G4UniformRand();
          value = Finvlor(pl,xe,r2);
          xco = Glor(pl,value)/c1;
          r3 = G4UniformRand();
        } while(r3>=xco);
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::SetPsi(G4double Energy, G4double PhiLocal)
{
  G4double value = 0.;
  G4double Ene = Energy/MeV;

  G4double p0l[4];
  G4double ppml[4];
  G4double p0t[2];
  G4double ppmt[2];

  G4double xi = 0.;
  G4double xe0 = 0.;
  G4double xepm = 0.;

  if (Ene>=50.)
    {
      const G4double ay00 = 3.4, by00 = 9.8, aa00 = 1.34, ba00 = 5.3;
      const G4double aw0 = 0.014, bw0 = 9.7, cw0 = -2.E4;
      const G4double axc0 = 3.1423, bxc0 = -2.35, cxc0 = 0.;
      const G4double ay0p = 1.53, by0p = 3.2, aap = 0.67, bap = 8.5E-3;
      const G4double awp = 6.9E-3, bwp = 12.6, cwp = -3.8E4;
      const G4double axcp = 2.8E-3,bxcp = -3.133;
      const G4double abf0 = 3.1213, bbf0 = 2.61;
      const G4double abfpm = 3.1231, bbfpm = 2.84;

      p0l[0] = Fln(ay00, by00, Ene);
      p0l[1] = Fln(aa00, ba00, Ene);
      p0l[2] = Poli(aw0, bw0, cw0, Ene);
      p0l[3] = Poli(axc0, bxc0, cxc0, Ene);

      ppml[0] = Fln(ay0p, by0p, Ene);
      ppml[1] = aap + bap*(Ene);
      ppml[2] = Poli(awp, bwp, cwp, Ene);
      ppml[3] = Fln(axcp,bxcp,Ene);

      p0t[0] = -0.81;
      p0t[1] = abf0 + bbf0/Ene;
      ppmt[0] = -0.6;
      ppmt[1] = abfpm + bbfpm/Ene;

      xi = 3.0;
      xe0 = Encu(p0l, p0t, xi);
      xepm = Encu(ppml, ppmt, xi);
    }
  else
    {
      const G4double ay00 = 2.82, by00 = 6.35;
      const G4double aa00 = -1.75, ba00 = 0.25;

      const G4double aw0 = 0.028, bw0 = 5., cw0 = -50.;
      const G4double axc0 = 3.14213, bxc0 = -2.3, cxc0 = 5.7;
      const G4double ay0p = 1.56, by0p = 3.6;
      const G4double aap = 0.86, bap = 8.3E-3;
      const G4double awp = 0.022, bwp = 7.4, cwp = -51.;
      const G4double xcp = 3.1486;

      p0l[0] = Fln(ay00, by00, Ene);
      p0l[1] = aa00+pow(Ene, ba00);
      p0l[2] = Poli(aw0, bw0, cw0, Ene);
      p0l[3] = Poli(axc0, bxc0, cxc0, Ene);
      ppml[0] = Fln(ay0p, by0p, Ene);
      ppml[1] = aap + bap*(Ene);
      ppml[2] = Poli(awp, bwp, cwp, Ene);
      ppml[3] = xcp;
    }

  G4double a,b=0.;

  if (Ene>=50.)
    {
      if (PhiLocal>xepm)
	{
          b = (ppml[0]+2*ppml[1]*ppml[2]*Flor(ppml,PhiLocal));
        }
      else
        {
          b = Ftan(ppmt,PhiLocal);
        }
      if (PhiLocal>xe0)
        {
          a = (p0l[0]+2*p0l[1]*p0l[2]*Flor(p0l,PhiLocal));
        }
      else
        {
          a = Ftan(p0t,PhiLocal);
        }
    }
  else
    {
      b = (ppml[0]+2*ppml[1]*ppml[2]*Flor(ppml,PhiLocal));
      a = (p0l[0]+2*p0l[1]*p0l[2]*Flor(p0l,PhiLocal));
    }
  G4double nr =0.;

  if (b>a)
    {
      nr = 1./b;
    }
  else
    {
      nr = 1./a;
    }

  G4double r1,r2=0.;
  G4double r3 =-1.;
  do
    {
      r1 = G4UniformRand();
      r2 = G4UniformRand();
      //value = r2*pi;
      value = 2.*r2*pi;
      r3 = nr*(a*cos(value)*cos(value) + b*sin(value)*sin(value));
    }while(r1>r3);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Poli
(G4double a, G4double b, G4double c, G4double x)
{
  G4double value=0.;
  if(x>0.)
    {
      value =(a + b/x + c/(x*x*x));
    }
  else
    {
      //G4cout << "ERROR in Poli! " << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Fln
(G4double a, G4double b, G4double x)
{
  G4double value=0.;
  if(x>0.)
    {
      value =(a*log(x)-b);
    }
  else
    {
      //G4cout << "ERROR in Fln! " << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Encu
(G4double* p_p1, G4double* p_p2, G4double x0)
{
  G4int i=0;
  G4double fx = 1.;
  G4double x = x0;
  const G4double xmax = 3.0;

  for(i=0; i<100; ++i)
    {
      fx = (Flor(p_p1,x)*Glor(p_p1,x) - Ftan(p_p2, x))/
	(Fdlor(p_p1,x) - Fdtan(p_p2,x));
      x -= fx;
      if(x > xmax) { return xmax; }
      if(std::fabs(fx) <= x*1.0e-6) { break; }
    } 

  if(x < 0.0) { x = 0.0; }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Flor(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double w = p_p1[2];
  G4double xc = p_p1[3];

  value = 1./(pi*(w*w + 4.*(x-xc)*(x-xc)));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Glor(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double y0 = p_p1[0];
  G4double A = p_p1[1];
  G4double w = p_p1[2];
  G4double xc = p_p1[3];

  value = (y0 *pi*(w*w +  4.*(x-xc)*(x-xc)) + 2.*A*w);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Fdlor(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double A = p_p1[1];
  G4double w = p_p1[2];
  G4double xc = p_p1[3];

  value = (-16.*A*w*(x-xc))/
    (pi*(w*w+4.*(x-xc)*(x-xc))*(w*w+4.*(x-xc)*(x-xc)));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Fintlor(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double y0 = p_p1[0];
  G4double A = p_p1[1];
  G4double w = p_p1[2];
  G4double xc = p_p1[3];

  value = y0*x + A*atan( 2*(x-xc)/w) / pi;
  return value;
}


G4double G4LivermorePolarizedGammaConversionModel::Finvlor(G4double* p_p1, G4double x, G4double r)
{
  G4double value = 0.;
  G4double nor = 0.;
  G4double w = p_p1[2];
  G4double xc = p_p1[3];

  nor = atan(2.*(pi-xc)/w)/(2.*pi*w) - atan(2.*(x-xc)/w)/(2.*pi*w);
  value = xc - (w/2.)*tan(-2.*r*nor*pi*w+atan(2*(xc-x)/w));

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Ftan(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double a = p_p1[0];
  G4double b = p_p1[1];

  value = a /(x-b);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Fdtan(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double a = p_p1[0];
  G4double b = p_p1[1];

  value = -1.*a / ((x-b)*(x-b));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Finttan(G4double* p_p1, G4double x)
{
  G4double value =0.;
  G4double a = p_p1[0];
  G4double b = p_p1[1];

  value = a*log(b-x);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedGammaConversionModel::Finvtan(G4double* p_p1, G4double cnor, G4double r)
{
  G4double value =0.;
  G4double a = p_p1[0];
  G4double b = p_p1[1];

  value = b*(1-G4Exp(r*cnor/a));

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedGammaConversionModel::SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedGammaConversionModel::GetRandomPolarization(G4ThreeVector& direction0)
{
  G4ThreeVector d0 = direction0.unit();
  G4ThreeVector a1 = SetPerpendicularVector(d0); //different orthogonal
  G4ThreeVector a0 = a1.unit(); // unit vector

  G4double rand1 = G4UniformRand();
  
  G4double angle = twopi*rand1; // random polar angle
  G4ThreeVector b0 = d0.cross(a0); // cross product
  
  G4ThreeVector c;
  
  c.setX(std::cos(angle)*(a0.x())+std::sin(angle)*b0.x());
  c.setY(std::cos(angle)*(a0.y())+std::sin(angle)*b0.y());
  c.setZ(std::cos(angle)*(a0.z())+std::sin(angle)*b0.z());
  
  G4ThreeVector c0 = c.unit();

  return c0;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedGammaConversionModel::GetPerpendicularPolarization
(const G4ThreeVector& gammaDirection, const G4ThreeVector& gammaPolarization) const
{
  // 
  // The polarization of a photon is always perpendicular to its momentum direction.
  // Therefore this function removes those vector component of gammaPolarization, which
  // points in direction of gammaDirection
  //
  // Mathematically we search the projection of the vector a on the plane E, where n is the
  // plains normal vector.
  // The basic equation can be found in each geometry book (e.g. Bronstein):
  // p = a - (a o n)/(n o n)*n
  
  return gammaPolarization - gammaPolarization.dot(gammaDirection)/gammaDirection.dot(gammaDirection) * gammaDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedGammaConversionModel::SystemOfRefChange
    (G4ThreeVector& direction0,G4ThreeVector& direction1,
     G4ThreeVector& polarization0)
{
  // direction0 is the original photon direction ---> z
  // polarization0 is the original photon polarization ---> x
  // need to specify y axis in the real reference frame ---> y 
  G4ThreeVector Axis_Z0 = direction0.unit();
  G4ThreeVector Axis_X0 = polarization0.unit();
  G4ThreeVector Axis_Y0 = (Axis_Z0.cross(Axis_X0)).unit(); // to be confirmed;
  
  G4double direction_x = direction1.getX();
  G4double direction_y = direction1.getY();
  G4double direction_z = direction1.getZ();
  
  direction1 = (direction_x*Axis_X0 + direction_y*Axis_Y0 +  direction_z*Axis_Z0).unit();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LivermorePolarizedGammaConversionModel::InitialiseForElement(
								      const G4ParticleDefinition*, 
								      G4int Z)
{
  G4AutoLock l(&LivermorePolarizedGammaConversionModelMutex);
  if(!data[Z]) { ReadData(Z); }
  l.unlock();
}
