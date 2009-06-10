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
// $Id: G4CrossSectionExcitationMillerGreenPartial.cc,v 1.4 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionExcitationMillerGreenPartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionExcitationMillerGreenPartial::G4CrossSectionExcitationMillerGreenPartial()
{
  nLevels = waterExcitation.NumberOfLevels();
  
  //PROTON
  kineticEnergyCorrection[0] = 1.;
  slaterEffectiveCharge[0][0] = 0.;
  slaterEffectiveCharge[1][0] = 0.;
  slaterEffectiveCharge[2][0] = 0.;
  sCoefficient[0][0] = 0.;
  sCoefficient[1][0] = 0.;
  sCoefficient[2][0] = 0.;

  //ALPHA++
  kineticEnergyCorrection[1] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][1]=0.;
  slaterEffectiveCharge[1][1]=0.;
  slaterEffectiveCharge[2][1]=0.;
  sCoefficient[0][1]=0.;
  sCoefficient[1][1]=0.;
  sCoefficient[2][1]=0.;

  // ALPHA+
  kineticEnergyCorrection[2] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][2]=2.0;
  slaterEffectiveCharge[1][2]=1.15;
  slaterEffectiveCharge[2][2]=1.15;
  sCoefficient[0][2]=0.7;
  sCoefficient[1][2]=0.15;
  sCoefficient[2][2]=0.15;

  // HELIUM
  kineticEnergyCorrection[3] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][3]=1.7;
  slaterEffectiveCharge[1][3]=1.15;
  slaterEffectiveCharge[2][3]=1.15;
  sCoefficient[0][3]=0.5;
  sCoefficient[1][3]=0.25;
  sCoefficient[2][3]=0.25;  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionExcitationMillerGreenPartial::~G4CrossSectionExcitationMillerGreenPartial()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::CrossSection(G4double k, G4int excitationLevel, 
								  const G4ParticleDefinition* particleDefinition)
{
  //                               ( ( z * aj ) ^ omegaj ) * ( t - ej ) ^ nu
  // sigma(t) = zEff^2 * sigma0 * --------------------------------------------
  //                               jj ^ ( omegaj + nu ) + t ^ ( omegaj + nu )
  //
  // where t is the kinetic energy corrected by Helium mass over proton mass for Helium ions
  //
  // zEff is:
  //  1 for protons
  //  2 for alpha++
  //  and  2 - c1 S_1s - c2 S_2s - c3 S_2p for alpha+ and He
  //
  // Dingfelder et al., RPC 59, 255-275, 2000 from Miller and Green (1973)
  // Formula (34) and Table 2
  
  const G4double sigma0(1.E+8 * barn);
  const G4double nu(1.);
  const G4double aj[]={876.*eV, 2084.* eV, 1373.*eV, 692.*eV, 900.*eV};
  const G4double jj[]={19820.*eV, 23490.*eV, 27770.*eV, 30830.*eV, 33080.*eV};
  const G4double omegaj[]={0.85, 0.88, 0.88, 0.78, 0.78};
  
  G4int particleTypeIndex = 0;
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition()) particleTypeIndex=0;
  if (particleDefinition == instance->GetIon("alpha++")) particleTypeIndex=1;
  if (particleDefinition == instance->GetIon("alpha+")) particleTypeIndex=2;
  if (particleDefinition == instance->GetIon("helium")) particleTypeIndex=3;
 
  G4double tCorrected;
  tCorrected = k * kineticEnergyCorrection[particleTypeIndex];

  // SI - added protection 
  if (tCorrected < waterExcitation.ExcitationEnergy(excitationLevel)) return 0;
  //
  
  G4int z = 10;

  G4double numerator;
  numerator = std::pow(z * aj[excitationLevel], omegaj[excitationLevel]) * 
    std::pow(tCorrected - waterExcitation.ExcitationEnergy(excitationLevel), nu);

  G4double power;
  power = omegaj[excitationLevel] + nu;

  G4double denominator;
  denominator = std::pow(jj[excitationLevel], power) + std::pow(tCorrected, power);

  G4double zEff = particleDefinition->GetPDGCharge() / eplus + particleDefinition->GetLeptonNumber();

  zEff -= ( sCoefficient[0][particleTypeIndex] * S_1s(k, waterExcitation.ExcitationEnergy(excitationLevel), slaterEffectiveCharge[0][particleTypeIndex], 1.) +
	    sCoefficient[1][particleTypeIndex] * S_2s(k, waterExcitation.ExcitationEnergy(excitationLevel), slaterEffectiveCharge[1][particleTypeIndex], 2.) +
	    sCoefficient[2][particleTypeIndex] * S_2p(k, waterExcitation.ExcitationEnergy(excitationLevel), slaterEffectiveCharge[2][particleTypeIndex], 2.) );

  G4double cross = sigma0 * zEff * zEff * numerator / denominator;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4CrossSectionExcitationMillerGreenPartial::RandomSelect(G4double k, 
							       const G4ParticleDefinition* particle)
{
  G4int i = nLevels;
  G4double value = 0.;
  std::deque<double> values;
  
  // ---- MGP ---- The following algorithm is wrong: it works if the cross section 
  // is a monotone increasing function.
  // The algorithm should be corrected by building the cumulative function 
  // of the cross section and comparing a random number in the range 0-1 against
  // the cumulative value at each bin 

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  // ELECTRON CORRECTION
 
  if ( particle == instance->GetIon("alpha++")||
       particle == G4Proton::ProtonDefinition())
       
  {  while (i > 0)
     {
	i--;
	G4double partial = CrossSection(k,i,particle);
	values.push_front(partial);
	value += partial;
     }

     value *= G4UniformRand();
    
     i = nLevels;

     while (i > 0)
     {
	i--;
	if (values[i] > value) return i;
	value -= values[i];
     }
  } 

  // add ONE or TWO electron-water excitation for alpha+ and helium
   
  if ( particle == instance->GetIon("alpha+") 
       ||
       particle == instance->GetIon("helium")
     ) 
  {
    while (i>0)
    {
	  i--;
         
	  G4CrossSectionExcitationEmfietzoglouPartial* excitationXS = 
	    new G4CrossSectionExcitationEmfietzoglouPartial();
         
	  G4double sigmaExcitation=0;
	  if (k*0.511/3728 > 7.4*eV && k*0.511/3728 < 10*keV) sigmaExcitation = excitationXS->CrossSection(k*0.511/3728,i);
 
	  G4double partial = CrossSection(k,i,particle);
	  if (particle == instance->GetIon("alpha+")) partial = CrossSection(k,i,particle) + sigmaExcitation;
	  if (particle == instance->GetIon("helium")) partial = CrossSection(k,i,particle) + 2*sigmaExcitation;
	  values.push_front(partial);
	  value += partial;
	  delete excitationXS;
    }
  
    value*=G4UniformRand();
  
    i=5;
    while (i>0)
    {
	  i--;
   
	  if (values[i]>value) return i;
  
	  value-=values[i];
    }
  }    
  //	

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::Sum(G4double k, const G4ParticleDefinition* particle)
{
  G4double totalCrossSection = 0.;

  for (G4int i=0; i<nLevels; i++)
  {
    totalCrossSection += CrossSection(k,i,particle);
  }
  return totalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::S_1s(G4double t,
							  G4double energyTransferred,
							  G4double slaterEffectiveCharge,
							  G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (7)
 
  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value = 1. - std::exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );
  
  return value;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::S_2s(G4double t,
							  G4double energyTransferred,
							  G4double slaterEffectiveCharge,
							  G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::S_2p(G4double t,
							  G4double energyTransferred,
							  G4double slaterEffectiveCharge,
							  G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionExcitationMillerGreenPartial::R(G4double t,
						       G4double energyTransferred,
						       G4double slaterEffectiveCharge,
						       G4double shellNumber) 
{
  // tElectron = m_electron / m_alpha * t
  // Hardcoded in Riccardo's implementation; to be corrected
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double tElectron = 0.511/3728. * t;
  G4double value = 2. * tElectron * slaterEffectiveCharge / (energyTransferred * shellNumber);
  
  return value;
}

