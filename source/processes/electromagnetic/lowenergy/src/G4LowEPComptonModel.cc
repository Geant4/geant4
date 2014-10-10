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
// ********************************************************************
// *********************************************************************
// |                                                                   |
// |             G4LowEPComptonModel-- Geant4 Monash University        |
// |                   low energy Compton scattering model.            |
// |             J. M. C. Brown, Monash University, Australia          |
// |                    ## Unpolarised photons only ##                 |
// |                                                                   |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | The following is a Geant4 class to simulate the process of        |
// | bound electron Compton scattering. General code structure is      |
// | based on G4LowEnergyCompton.cc and G4LivermoreComptonModel.cc.    |
// | Algorithms for photon energy, and ejected Compton electron        |
// | direction taken from:                                             |
// |                                                                   |
// | J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin,    |
// | "A low energy bound atomic electron Compton scattering model      |
// |  for Geant4", NIMA, submitted 2013.                               |
// |                                                                   |
// | The author acknowledges the work of the Geant4 collaboration      |
// | in developing the following algorithms that have been employed    |
// | or adapeted for the present software:                             |    
// |                                                                   |
// |  # sampling of photon scattering angle,                           |
// |  # target element selection in composite materials,               |
// |  # target shell selection in element,                             |
// |  # and sampling of bound electron momentum from Compton profiles. |
// |                                                                   |
// *********************************************************************
// |                                                                   |
// | History:                                                          |
// | --------                                                          |
// |                                                                   |
// | Nov. 2011 JMCB       - First version                              |
// | Feb. 2012 JMCB       - Migration to Geant4 9.5                    |
// | Sep. 2012 JMCB       - Final fixes for Geant4 9.6                 |
// | Feb. 2013 JMCB       - Geant4 9.6 FPE fix for bug 1426            |
// |                                                                   |
// *********************************************************************

#include <limits>
#include "G4LowEPComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4CrossSectionHandler.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Gamma.hh"

//****************************************************************************

using namespace std;

//****************************************************************************

G4LowEPComptonModel::G4LowEPComptonModel(const G4ParticleDefinition*,
						 const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false),
   scatterFunctionData(0),crossSectionHandler(0),fAtomDeexcitation(0)
{
  lowEnergyLimit = 250 * eV; 
  highEnergyLimit = 100 * GeV;

  verboseLevel=0 ;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(  verboseLevel>0 ) { 
    G4cout << "Low energy photon Compton model is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

}

//****************************************************************************

G4LowEPComptonModel::~G4LowEPComptonModel()
{  
  delete crossSectionHandler;
  delete scatterFunctionData;
}

//****************************************************************************

void G4LowEPComptonModel::Initialise(const G4ParticleDefinition* particle,
					 const G4DataVector& cuts)
{
  if (verboseLevel > 2) {
    G4cout << "Calling G4LowEPComptonModel::Initialise()" << G4endl;
  }

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }
  delete scatterFunctionData;

  // Reading of data files - all materials are read  
  crossSectionHandler = new G4CrossSectionHandler;
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(scatterFile);

  // For Doppler broadening
  shellData.SetOccupancyData();
  G4String file = "/doppler/shell-doppler";
  shellData.LoadData(file);

  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for low energy photon Compton model" << G4endl;
  }

  if(isInitialised) { return; }
  isInitialised = true;

  fParticleChange = GetParticleChangeForGamma();

  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

  if(  verboseLevel>0 ) { 
    G4cout << "Low energy photon Compton model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }  
}

//****************************************************************************

G4double G4LowEPComptonModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LowEPComptonModel" << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) { return 0.0; }
    
  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);  
  return cs;
}





//****************************************************************************


void G4LowEPComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						const G4MaterialCutsCouple* couple,
						const G4DynamicParticle* aDynamicGamma,
						G4double, G4double)
{

  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // then accepted or rejected depending on the Scattering Function multiplied
  // by factor from Klein - Nishina formula.
  // Expression of the angular distribution as Klein Nishina
  // angular and energy distribution and Scattering fuctions is taken from
  // D. E. Cullen "A simple model of photon transport" Nucl. Instr. Meth.
  // Phys. Res. B 101 (1995). Method of sampling with form factors is different
  // data are interpolated while in the article they are fitted.
  // Reference to the article is from J. Stepanek New Photon, Positron
  // and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10
  // TeV (draft).
  // The random number techniques of Butcher & Messel are used
  // (Nucl Phys 20(1960),15).

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy()/MeV;

  if (verboseLevel > 3) {
    G4cout << "G4LowEPComptonModel::SampleSecondaries() E(MeV)= " 
	   << photonEnergy0/MeV << " in " << couple->GetMaterial()->GetName() 
	   << G4endl;
  }
  
  // low-energy gamma is absorpted by this process
  if (photonEnergy0 <= lowEnergyLimit) 
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double LowEPCepsilon0 = 1. / (1. + 2. * e0m);
  G4double LowEPCepsilon0Sq = LowEPCepsilon0 * LowEPCepsilon0;
  G4double alpha1 = -std::log(LowEPCepsilon0);
  G4double alpha2 = 0.5 * (1. - LowEPCepsilon0Sq);

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  // Sample the energy of the scattered photon
  G4double LowEPCepsilon;
  G4double LowEPCepsilonSq;
  G4double oneCosT;
  G4double sinT2;
  G4double gReject;
 
  do
    {
      if ( alpha1/(alpha1+alpha2) > G4UniformRand())
	{
	  LowEPCepsilon = std::exp(-alpha1 * G4UniformRand());  
	  LowEPCepsilonSq = LowEPCepsilon * LowEPCepsilon;
	}
      else
	{
	  LowEPCepsilonSq = LowEPCepsilon0Sq + (1. - LowEPCepsilon0Sq) * G4UniformRand();
	  LowEPCepsilon = std::sqrt(LowEPCepsilonSq);
	}

      oneCosT = (1. - LowEPCepsilon) / ( LowEPCepsilon * e0m);
      sinT2 = oneCosT * (2. - oneCosT);
      G4double x = std::sqrt(oneCosT/2.) / (wlPhoton/cm);
      G4double scatteringFunction = scatterFunctionData->FindValue(x,Z-1);
      gReject = (1. - LowEPCepsilon * sinT2 / (1. + LowEPCepsilonSq)) * scatteringFunction;

    } while(gReject < G4UniformRand()*Z); 

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = std::sqrt(sinT2);
  G4double phi = twopi * G4UniformRand();
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  
  // Scatter photon energy and Compton electron direction - Method based on:
  // J. M. C. Brown, M. R. Dimmock, J. E. Gillam and D. M. Paganin'
  // "A low energy bound atomic electron Compton scattering model for Geant4"
  // NIMA ISSUE, PG, 2013
  
  // Set constants and initialize scattering parameters

  const G4double vel_c = c_light / (m/s);
  const G4double momentum_au_to_nat = halfpi* hbar_Planck / Bohr_radius / (kg*m/s);
  const G4double e_mass_kg =  electron_mass_c2 / c_squared / kg ;
    
  const G4int maxDopplerIterations = 1000;  
  G4double bindingE = 0.; 
  G4double pEIncident = photonEnergy0 ;
  G4double pERecoil =  -1.;
  G4double eERecoil = -1.;
  G4double e_alpha =0.;
  G4double e_beta = 0.;
 
  G4double CE_emission_flag = 0.;
  G4double ePAU = -1;
  G4int shellIdx = 0;
  G4double u_temp = 0;
  G4double cosPhiE =0;
  G4double sinThetaE =0;
  G4double cosThetaE =0;
  G4int iteration = 0; 
  do{
    
      
      // ******************************************
      // |     Determine scatter photon energy    |
      // ******************************************      
   
  do
    {
      iteration++;
      
      
      // ********************************************
      // |     Sample bound electron information    |
      // ********************************************
      
      // Select shell based on shell occupancy
      
      shellIdx = shellData.SelectRandomShell(Z);
      bindingE = shellData.BindingEnergy(Z,shellIdx)/MeV; 
      
      
      // Randomly sample bound electron momentum (memento: the data set is in Atomic Units)
      ePAU = profileData.RandomSelectMomentum(Z,shellIdx);  

      // Convert to SI units     
      G4double ePSI = ePAU * momentum_au_to_nat;
      
      //Calculate bound electron velocity and normalise to natural units
      u_temp = sqrt( ((ePSI*ePSI)*(vel_c*vel_c)) / ((e_mass_kg*e_mass_kg)*(vel_c*vel_c)+(ePSI*ePSI)) )/vel_c;  

      // Sample incident electron direction, amorphous material, to scattering photon scattering plane      

      e_alpha = pi*G4UniformRand();
      e_beta = twopi*G4UniformRand();   
      
      // Total energy of system  
      
      G4double eEIncident = electron_mass_c2 / sqrt( 1 - (u_temp*u_temp));
      G4double systemE = eEIncident + pEIncident;


      G4double gamma_temp = 1.0 / sqrt( 1 - (u_temp*u_temp));
      G4double numerator = gamma_temp*electron_mass_c2*(1 - u_temp * std::cos(e_alpha));
      G4double subdenom1 =  u_temp*cosTheta*std::cos(e_alpha);
      G4double subdenom2 = u_temp*sinTheta*std::sin(e_alpha)*std::cos(e_beta);
      G4double denominator = (1.0 - cosTheta) +  (gamma_temp*electron_mass_c2*(1 - subdenom1 - subdenom2) / pEIncident);
      pERecoil = (numerator/denominator);
      eERecoil = systemE - pERecoil; 
      CE_emission_flag = pEIncident - pERecoil;
    } while ( (iteration <= maxDopplerIterations) && (CE_emission_flag < bindingE));      
    
    
   
  // End of recalculation of photon energy with Doppler broadening



   // *******************************************************
   // |     Determine ejected Compton electron direction    |
   // *******************************************************      
      
      // Calculate velocity of ejected Compton electron   
      
      G4double a_temp = eERecoil / electron_mass_c2;
      G4double u_p_temp = sqrt(1 - (1 / (a_temp*a_temp)));

      // Coefficients and terms from simulatenous equations     
     
      G4double sinAlpha = std::sin(e_alpha);
      G4double cosAlpha = std::cos(e_alpha);
      G4double sinBeta = std::sin(e_beta);
      G4double cosBeta = std::cos(e_beta);
      
      G4double gamma = 1.0 / sqrt(1 - (u_temp*u_temp));
      G4double gamma_p = 1.0 / sqrt(1 - (u_p_temp*u_p_temp));
      
      G4double var_A = pERecoil*u_p_temp*sinTheta;
      G4double var_B = u_p_temp* (pERecoil*cosTheta-pEIncident);
      G4double var_C = (pERecoil-pEIncident) - ( (pERecoil*pEIncident) / (gamma_p*electron_mass_c2))*(1 - cosTheta);

      G4double var_D1 = gamma*electron_mass_c2*pERecoil;
      G4double var_D2 = (1 - (u_temp*cosTheta*cosAlpha) - (u_temp*sinTheta*cosBeta*sinAlpha));
      G4double var_D3 = ((electron_mass_c2*electron_mass_c2)*(gamma*gamma_p - 1)) - (gamma_p*electron_mass_c2*pERecoil);
      G4double var_D = var_D1*var_D2 + var_D3;     

      G4double var_E1 = ((gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*cosAlpha);
      G4double var_E2 = gamma_p*electron_mass_c2*pERecoil*u_p_temp*cosTheta;
      G4double var_E = var_E1 - var_E2;
      
      G4double var_F1 = ((gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*cosBeta*sinAlpha);
      G4double var_F2 = (gamma_p*electron_mass_c2*pERecoil*u_p_temp*sinTheta);
      G4double var_F = var_F1 - var_F2;
      
      G4double var_G = (gamma*gamma_p)*(electron_mass_c2*electron_mass_c2)*(u_temp*u_p_temp)*sinBeta*sinAlpha;
      
      // Two equations form a quadratic form of Wx^2 + Yx + Z = 0
      // Coefficents and solution to quadratic
      
      G4double var_W1 = (var_F*var_B - var_E*var_A)*(var_F*var_B - var_E*var_A);
      G4double var_W2 = (var_G*var_G)*(var_A*var_A) + (var_G*var_G)*(var_B*var_B);
      G4double var_W = var_W1 + var_W2;
      
      G4double var_Y = 2.0*(((var_A*var_D-var_F*var_C)*(var_F*var_B-var_E*var_A)) - ((var_G*var_G)*var_B*var_C));
      
      G4double var_Z1 = (var_A*var_D - var_F*var_C)*(var_A*var_D - var_F*var_C);
      G4double var_Z2 = (var_G*var_G)*(var_C*var_C) - (var_G*var_G)*(var_A*var_A);
      G4double var_Z = var_Z1 + var_Z2;
      G4double diff1 = var_Y*var_Y;
      G4double diff2 = 4*var_W*var_Z;      
      G4double diff = diff1 - diff2;
      
      
     // Check if diff is less than zero, if so ensure it is due to FPE
      
     //Determine number of digits (in decimal base) that G4double can accurately represent
     G4double g4d_order = G4double(numeric_limits<G4double>::digits10);      
     G4double g4d_limit = std::pow(10.,-g4d_order);
     //Confirm that diff less than zero is due FPE, i.e if abs of diff / diff1 and diff/ diff2 is less 
     //than 10^(-g4d_order), then set diff to zero
     
     if ((diff < 0.0) && (abs(diff / diff1) < g4d_limit) && (abs(diff / diff2) < g4d_limit) )    
     {
	   diff = 0.0;         
     }

  
      // Plus and minus of quadratic
      G4double X_p = (-var_Y + sqrt (diff))/(2*var_W);
      G4double X_m = (-var_Y - sqrt (diff))/(2*var_W);


      // Randomly sample one of the two possible solutions and determin theta angle of ejected Compton electron
      G4double ThetaE = 0.;
      G4double sol_select = G4UniformRand();
      
      if (sol_select < 0.5)
      {
          ThetaE = std::acos(X_p);
      }
      if (sol_select > 0.5)
      {
          ThetaE = std::acos(X_m);
      }
      
      cosThetaE = std::cos(ThetaE);
      sinThetaE = std::sin(ThetaE);
      G4double Theta = std::acos(cosTheta);
      
      //Calculate electron Phi
      G4double iSinThetaE = std::sqrt(1+std::tan((pi/2.0)-ThetaE)*std::tan((pi/2.0)-ThetaE));
      G4double iSinTheta = std::sqrt(1+std::tan((pi/2.0)-Theta)*std::tan((pi/2.0)-Theta)); 
      G4double ivar_A = iSinTheta/ (pERecoil*u_p_temp);
      // Trigs
      cosPhiE = (var_C - var_B*cosThetaE)*(ivar_A*iSinThetaE); 
      
     // End of calculation of ejection Compton electron direction
     
      //Fix for floating point errors

    } while ( (iteration <= maxDopplerIterations) && (abs(cosPhiE) > 1));        
      
   // Revert to original if maximum number of iterations threshold has been reached     
      
  
  if (iteration >= maxDopplerIterations)
    {
      pERecoil = photonEnergy0 ;
      bindingE = 0.;
      dirx=0.0;
      diry=0.0;
      dirz=1.0;
    }
  
  // Set "scattered" photon direction and energy
  
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;

  if (pERecoil > 0.)
    {
     fParticleChange->SetProposedKineticEnergy(pERecoil) ;

     // Set ejected Compton electron direction and energy
     G4double PhiE = std::acos(cosPhiE);
     G4double eDirX = sinThetaE * std::cos(phi+PhiE);
     G4double eDirY = sinThetaE * std::sin(phi+PhiE);
     G4double eDirZ = cosThetaE;
  
     G4double eKineticEnergy = pEIncident - pERecoil - bindingE;  
  
     G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
     eDirection.rotateUz(photonDirection0);
     G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
						   eDirection,eKineticEnergy) ;
     fvect->push_back(dp);

    }
  else
    {
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeTrackStatus(fStopAndKill);   
    }
      



  // sample deexcitation
  
  if(fAtomDeexcitation && iteration < maxDopplerIterations) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t i=nbefore; i<nafter; ++i) {
	  //Check if there is enough residual energy
          if (bindingE >= ((*fvect)[i])->GetKineticEnergy())
           {
             //Ok, this is a valid secondary: keep it
             bindingE -= ((*fvect)[i])->GetKineticEnergy();
           }
          else
           {
             //Invalid secondary: not enough energy to create it!
             //Keep its energy in the local deposit
             delete (*fvect)[i];
             (*fvect)[i]=0;
           }
	} 
      }
    }
  }

  //This should never happen
  if(bindingE < 0.0)
     G4Exception("G4LowEPComptonModel::SampleSecondaries()",
                 "em2051",FatalException,"Negative local energy deposit");

  fParticleChange->ProposeLocalEnergyDeposit(bindingE);
  
}

