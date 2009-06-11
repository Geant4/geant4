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
// $Id: G4FinalStateIonisationRudd.cc,v 1.11 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4FinalStateIonisationRudd.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateIonisationRudd::G4FinalStateIonisationRudd()
{
  lowEnergyLimitDefault = 100 * eV;
  highEnergyLimitDefault = 100 * MeV;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");

  G4String proton;
  G4String hydrogen;
  G4String alphaPlusPlus;
  G4String alphaPlus;
  G4String helium;

  proton = protonDef->GetParticleName();
  lowEnergyLimit[proton] = 100. * eV;
  highEnergyLimit[proton] = 500. * keV;

  hydrogen = hydrogenDef->GetParticleName();
  lowEnergyLimit[hydrogen] = 100. * eV;
  highEnergyLimit[hydrogen] = 100. * MeV;

  alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlusPlus] = 1. * keV;
  highEnergyLimit[alphaPlusPlus] = 10. * MeV;

  alphaPlus = alphaPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlus] = 1. * keV;
  highEnergyLimit[alphaPlus] = 10. * MeV;

  helium = heliumDef->GetParticleName();
  lowEnergyLimit[helium] = 1. * keV;
  highEnergyLimit[helium] = 10. * MeV;

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4FinalStateIonisationRudd is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateIonisationRudd::~G4FinalStateIonisationRudd()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4FinalStateProduct& G4FinalStateIonisationRudd::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  product.Clear();

  const G4DynamicParticle* particle = track.GetDynamicParticle();

  G4double lowLim = lowEnergyLimitDefault;
  G4double highLim = highEnergyLimitDefault;

  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  if (k >= lowLim && k <= highLim)
  {
      G4ParticleDefinition* definition = particle->GetDefinition();
      G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
      G4double particleMass = definition->GetPDGMass();
      G4double totalEnergy = k + particleMass;
      G4double pSquare = k*(totalEnergy+particleMass);
      G4double totalMomentum = std::sqrt(pSquare);

      const G4String& particleName = definition->GetParticleName();
  
      G4int ionizationShell = cross.RandomSelect(k,particleName);
  
      G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);
  
      G4double bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

      G4double cosTheta = 0.;
      G4double phi = 0.; 
      RandomizeEjectedElectronDirection(definition, k,secondaryKinetic, cosTheta, phi);

      G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
      G4double dirX = sinTheta*std::cos(phi);
      G4double dirY = sinTheta*std::sin(phi);
      G4double dirZ = cosTheta;
      G4ThreeVector deltaDirection(dirX,dirY,dirZ);
      deltaDirection.rotateUz(primaryDirection);

      G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));

      G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
      G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
      G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;

      product.ModifyPrimaryParticle(finalPx,finalPy,finalPz,k-bindingEnergy-secondaryKinetic);
      product.AddEnergyDeposit(bindingEnergy);

      G4DynamicParticle* aElectron = new G4DynamicParticle(G4Electron::Electron(),deltaDirection,secondaryKinetic);
      product.AddSecondary(aElectron);
  }

  if (k < lowLim) 
  {  
    product.KillPrimaryParticle();
  }
 
  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
								    G4double k, 
								    G4int shell)
{
  G4double maximumKineticEnergyTransfer = 0.;
 
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition() 
      || particleDefinition == instance->GetIon("hydrogen"))
  { 
      maximumKineticEnergyTransfer= 4.* (electron_mass_c2 / proton_mass_c2) * k;
  }

  if (particleDefinition == instance->GetIon("helium") 
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("alpha++"))
  { 
      maximumKineticEnergyTransfer= 4.* (0.511 / 3728) * k;
  }

  G4double crossSectionMaximum = 0.;
  
  for(G4double value=waterStructure.IonisationEnergy(shell); value<=4.*waterStructure.IonisationEnergy(shell) ; value+=0.1*eV)
  {
      G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k, value, shell);
      if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
  }
  
  G4double secElecKinetic = 0.;
  
  do
  {
    secElecKinetic = G4UniformRand() * maximumKineticEnergyTransfer;
  } while(G4UniformRand()*crossSectionMaximum > DifferentialCrossSection(particleDefinition, 
									 k,
									 secElecKinetic+waterStructure.IonisationEnergy(shell),
									 shell));

  return(secElecKinetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4FinalStateIonisationRudd::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
								   G4double k, 
								   G4double secKinetic, 
								   G4double & cosTheta, 
								   G4double & phi )
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double maxSecKinetic = 0.;
 
  if (particleDefinition == G4Proton::ProtonDefinition() 
      || particleDefinition == instance->GetIon("hydrogen")) 
  { 
      maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
  }
  
  if (particleDefinition == instance->GetIon("helium") 
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("alpha++"))
  { 
      maxSecKinetic = 4.* (0.511 / 3728) * k;
  }
  
  phi = twopi * G4UniformRand();
  cosTheta = std::sqrt(secKinetic / maxSecKinetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4FinalStateIonisationRudd::DifferentialCrossSection(G4ParticleDefinition* particleDefinition, 
							      G4double k, 
							      G4double energyTransfer, 
							      G4int ionizationLevelIndex)
{
  // Shells ids are 0 1 2 3 4 (4 is k shell)
  // !!Attention, "energyTransfer" here is the energy transfered to the electron which means
  //             that the secondary kinetic energy is w = energyTransfer - bindingEnergy
  //
  //   ds            S                F1(nu) + w * F2(nu)
  //  ---- = G(k) * ----     -------------------------------------------
  //   dw            Bj       (1+w)^3 * [1 + exp{alpha * (w - wc) / nu}]
  //
  // w is the secondary electron kinetic Energy in eV
  //
  // All the other parameters can be found in Rudd's Papers
  //
  // M.Eugene Rudd, 1988, User-Friendly model for the energy distribution of
  // electrons from protons or electron collisions. Nucl. Tracks Rad. Meas.Vol 16 N0 2/3 pp 219-218
  //

  const G4int j=ionizationLevelIndex;

  G4double A1 ; 
  G4double B1 ; 
  G4double C1 ; 
  G4double D1 ; 
  G4double E1 ;
  G4double A2 ; 
  G4double B2 ; 
  G4double C2 ; 
  G4double D2 ; 
  G4double alphaConst ;

  if (j == 4) 
  {
      //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
      A1 = 1.25; 
      B1 = 0.5; 
      C1 = 1.00; 
      D1 = 1.00; 
      E1 = 3.00; 
      A2 = 1.10; 
      B2 = 1.30;
      C2 = 1.00; 
      D2 = 0.00; 
      alphaConst = 0.66;
  }
  else 
  {
      //Data For Liquid Water from Dingfelder (Protons in Water)
      A1 = 1.02; 
      B1 = 82.0; 
      C1 = 0.45; 
      D1 = -0.80; 
      E1 = 0.38; 
      A2 = 1.07; 
      B2 = 14.6;
      C2 = 0.60; 
      D2 = 0.04; 
      alphaConst = 0.64;
  }
  
  const G4double n = 2.;
  const G4double Gj[5] = {0.99, 1.11, 1.11, 0.52, 1.};

  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double wBig = (energyTransfer - waterStructure.IonisationEnergy(ionizationLevelIndex));
  G4double w = wBig / waterStructure.IonisationEnergy(ionizationLevelIndex);
  G4double Ry = 13.6*eV;

  G4double tau = 0.;

  if (particleDefinition == G4Proton::ProtonDefinition() 
      || particleDefinition == instance->GetIon("hydrogen")) 
  {
      tau = (electron_mass_c2/proton_mass_c2) * k ;
  }
   
  if ( particleDefinition == instance->GetIon("helium") 
       || particleDefinition == instance->GetIon("alpha+")
       || particleDefinition == instance->GetIon("alpha++")) 
  {
      tau = (0.511/3728.) * k ;
  }
 
  G4double S = 4.*pi * Bohr_radius*Bohr_radius * n * std::pow((Ry/waterStructure.IonisationEnergy(ionizationLevelIndex)),2);
  G4double v2 = tau / waterStructure.IonisationEnergy(ionizationLevelIndex);
  G4double v = std::sqrt(v2);
  G4double wc = 4.*v2 - 2.*v - (Ry/(4.*waterStructure.IonisationEnergy(ionizationLevelIndex)));

  G4double L1 = (C1* std::pow(v,(D1))) / (1.+ E1*std::pow(v, (D1+4.)));
  G4double L2 = C2*std::pow(v,(D2));
  G4double H1 = (A1*std::log(1.+v2)) / (v2+(B1/v2));
  G4double H2 = (A2/v2) + (B2/(v2*v2));

  G4double F1 = L1+H1;
  G4double F2 = (L2*H2)/(L2+H2);

  G4double sigma = CorrectionFactor(particleDefinition, k/eV) 
    * Gj[j] * (S/waterStructure.IonisationEnergy(ionizationLevelIndex)) 
    * ( (F1+w*F2) / ( std::pow((1.+w),3) * ( 1.+std::exp(alphaConst*(w-wc)/v))) );
    
  if (    particleDefinition == G4Proton::ProtonDefinition() 
	  || particleDefinition == instance->GetIon("hydrogen")
	  ) 
  {
      return(sigma);
  }

  if (particleDefinition == instance->GetIon("alpha++") ) 
  {
      slaterEffectiveCharge[0]=0.;
      slaterEffectiveCharge[1]=0.;
      slaterEffectiveCharge[2]=0.;
      sCoefficient[0]=0.;
      sCoefficient[1]=0.;
      sCoefficient[2]=0.;
  }

  if (particleDefinition == instance->GetIon("alpha+") ) 
  {
      slaterEffectiveCharge[0]=2.0;
      slaterEffectiveCharge[1]=1.15;
      slaterEffectiveCharge[2]=1.15;
      sCoefficient[0]=0.7;
      sCoefficient[1]=0.15;
      sCoefficient[2]=0.15;
  }

  if (particleDefinition == instance->GetIon("helium") ) 
  {
      slaterEffectiveCharge[0]=1.7;
      slaterEffectiveCharge[1]=1.15;
      slaterEffectiveCharge[2]=1.15;
      sCoefficient[0]=0.5;
      sCoefficient[1]=0.25;
      sCoefficient[2]=0.25;
  }
  
  if (    particleDefinition == instance->GetIon("helium") 
	  || particleDefinition == instance->GetIon("alpha+")
	  || particleDefinition == instance->GetIon("alpha++")
	  ) 
  {
      sigma = Gj[j] * (S/waterStructure.IonisationEnergy(ionizationLevelIndex)) * ( (F1+w*F2) / ( std::pow((1.+w),3) * ( 1.+std::exp(alphaConst*(w-wc)/v))) );
    
      G4double zEff = particleDefinition->GetPDGCharge() / eplus + particleDefinition->GetLeptonNumber();
  
      zEff -= ( sCoefficient[0] * S_1s(k, energyTransfer, slaterEffectiveCharge[0], 1.) +
		sCoefficient[1] * S_2s(k, energyTransfer, slaterEffectiveCharge[1], 2.) +
		sCoefficient[2] * S_2p(k, energyTransfer, slaterEffectiveCharge[2], 2.) );
	   
      return zEff * zEff * sigma ;
  }  
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::S_1s(G4double t, 
					  G4double energyTransferred, 
					  G4double slaterEffectiveChg, 
					  G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (7)
 
  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value = 1. - std::exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::S_2s(G4double t,
					  G4double energyTransferred, 
					  G4double slaterEffectiveChg, 
					  G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::S_2p(G4double t, 
					  G4double energyTransferred,
					  G4double slaterEffectiveChg, 
					  G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::R(G4double t,
				       G4double energyTransferred,
				       G4double slaterEffectiveChg,
				       G4double shellNumber) 
{
  // tElectron = m_electron / m_alpha * t
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double tElectron = 0.511/3728. * t;
  G4double value = 2. * tElectron * slaterEffectiveChg / (energyTransferred * shellNumber);
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateIonisationRudd::CorrectionFactor(G4ParticleDefinition* particleDefinition, G4double k) 
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::Proton()) 
  {
      return(1.);
  }
  else 
    if (particleDefinition == instance->GetIon("hydrogen")) 
    { 
	G4double value = (std::log(k/eV)-4.2)/0.5;
	return((0.8/(1+std::exp(value))) + 0.9);
    }
    else 
    {    
	return(1.);
    }
}
