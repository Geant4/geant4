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
// $Id: G4DNARuddIonisationExtendedModel.cc,v 1.3 2010-11-04 14:52:17 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


// Modified by Z. Francis to handle HZE && inverse rudd function sampling 26-10-2010

#include "G4DNARuddIonisationExtendedModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::G4DNARuddIonisationExtendedModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

  lowEnergyLimitForA[1] = 0 * eV; 
  lowEnergyLimitForA[2] = 0 * eV;
  lowEnergyLimitForA[3] = 0 * eV; 
  lowEnergyLimitOfModelForA[1] = 100 * eV; 
  lowEnergyLimitOfModelForA[4] = 1 * keV; 
  lowEnergyLimitOfModelForA[5] = 0.5 * MeV; // For A = 3 or above, limit is MeV/uma
  killBelowEnergyForA[1] = lowEnergyLimitOfModelForA[1]; 
  killBelowEnergyForA[4] = lowEnergyLimitOfModelForA[4]; 
  killBelowEnergyForA[5] = lowEnergyLimitOfModelForA[5];


  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if( verboseLevel>0 ) 
  { 
    G4cout << "Rudd ionisation model is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationExtendedModel::~G4DNARuddIonisationExtendedModel()
{  
  // Cross section
  
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4DNARuddIonisationExtendedModel::Initialise()" << G4endl;

  // Energy limits

  G4String fileProton("dna/sigma_ionisation_p_rudd");
  G4String fileHydrogen("dna/sigma_ionisation_h_rudd");
  G4String fileAlphaPlusPlus("dna/sigma_ionisation_alphaplusplus_rudd");
  G4String fileAlphaPlus("dna/sigma_ionisation_alphaplus_rudd");
  G4String fileHelium("dna/sigma_ionisation_he_rudd");
  G4String fileCarbon("dna/sigma_ionisation_c_rudd");
  G4String fileNitrogen("dna/sigma_ionisation_n_rudd");
  G4String fileOxygen("dna/sigma_ionisation_o_rudd");
  G4String fileIron("dna/sigma_ionisation_fe_rudd");

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");
  G4ParticleDefinition* carbonDef = instance->GetIon("carbon");
  G4ParticleDefinition* nitrogenDef = instance->GetIon("nitrogen");
  G4ParticleDefinition* oxygenDef = instance->GetIon("oxygen");
  G4ParticleDefinition* ironDef = instance->GetIon("iron");

  G4String proton;
  G4String hydrogen;
  G4String alphaPlusPlus;
  G4String alphaPlus;
  G4String helium;
  G4String carbon;
  G4String nitrogen;
  G4String oxygen;
  G4String iron;

  G4double scaleFactor = 1 * m*m;

  if (protonDef != 0)
  {
    proton = protonDef->GetParticleName();
    tableFile[proton] = fileProton;
    lowEnergyLimit[proton] = lowEnergyLimitForA[1];
    highEnergyLimit[proton] = 500. * keV;

    // Cross section

    G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
    tableProton->LoadData(fileProton);
    tableData[proton] = tableProton;
  }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: proton is not defined");
  }

  if (hydrogenDef != 0)
  {
    hydrogen = hydrogenDef->GetParticleName();
    tableFile[hydrogen] = fileHydrogen;

    lowEnergyLimit[hydrogen] = lowEnergyLimitForA[1];
    highEnergyLimit[hydrogen] = 100. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
									     eV,
									     scaleFactor );
    tableHydrogen->LoadData(fileHydrogen);
      
    tableData[hydrogen] = tableHydrogen;
  }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: hydrogen is not defined");
  }

  if (alphaPlusPlusDef != 0)
  {
    alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
    tableFile[alphaPlusPlus] = fileAlphaPlusPlus;

    lowEnergyLimit[alphaPlusPlus] = lowEnergyLimitForA[4];
    highEnergyLimit[alphaPlusPlus] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
										  eV,
										  scaleFactor );
    tableAlphaPlusPlus->LoadData(fileAlphaPlusPlus);
      
    tableData[alphaPlusPlus] = tableAlphaPlusPlus;
  }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: alphaPlusPlus is not defined");
  }

  if (alphaPlusDef != 0)
  {
    alphaPlus = alphaPlusDef->GetParticleName();
    tableFile[alphaPlus] = fileAlphaPlus;

    lowEnergyLimit[alphaPlus] = lowEnergyLimitForA[4];
    highEnergyLimit[alphaPlus] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									      eV,
									      scaleFactor );
    tableAlphaPlus->LoadData(fileAlphaPlus);
    tableData[alphaPlus] = tableAlphaPlus;
  }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: alphaPlus is not defined");
  }

  if (heliumDef != 0)
  {
    helium = heliumDef->GetParticleName();
    tableFile[helium] = fileHelium;

    lowEnergyLimit[helium] = lowEnergyLimitForA[4];
    highEnergyLimit[helium] = 400. * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
    tableHelium->LoadData(fileHelium);
    tableData[helium] = tableHelium;
    }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: helium is not defined");
  }

 if (carbonDef != 0)
  {
    carbon = carbonDef->GetParticleName();
    tableFile[carbon] = fileCarbon;

    lowEnergyLimit[carbon] = lowEnergyLimitForA[5] * particle->GetAtomicMass();
    highEnergyLimit[carbon] = 1e6* particle->GetAtomicMass() * MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableCarbon = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
    tableCarbon->LoadData(fileCarbon);
    tableData[carbon] = tableCarbon;
    }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: Carbon is not defined");
  }

 if (oxygenDef != 0)
  {
    oxygen = oxygenDef->GetParticleName();
    tableFile[oxygen] = fileOxygen;

    lowEnergyLimit[oxygen] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    highEnergyLimit[oxygen] = 1e6* particle->GetAtomicMass()* MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableOxygen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
      tableOxygen->LoadData(fileOxygen);
      tableData[oxygen] = tableOxygen;
    }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: Oxygen is not defined");
  }

  if (nitrogenDef != 0)
  {
    nitrogen = nitrogenDef->GetParticleName();
    tableFile[nitrogen] = fileNitrogen;

    lowEnergyLimit[nitrogen] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    highEnergyLimit[nitrogen] = 1e6* particle->GetAtomicMass()* MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableNitrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
      tableNitrogen->LoadData(fileNitrogen);
      tableData[nitrogen] = tableNitrogen;
    }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: Nitrogen is not defined");
  }

  if (ironDef != 0)
  {
    iron = ironDef->GetParticleName();
    tableFile[iron] = fileIron;

    lowEnergyLimit[iron] = lowEnergyLimitForA[5]* particle->GetAtomicMass();
    highEnergyLimit[iron] = 1e6* particle->GetAtomicMass()* MeV;

    // Cross section

    G4DNACrossSectionDataSet* tableIron = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
									   eV,
									   scaleFactor );
      tableIron->LoadData(fileIron);
      tableData[iron] = tableIron;
    }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::Initialise: Iron is not defined");
  }

// ZF Following lines can be replaced by:
    SetLowEnergyLimit(lowEnergyLimit[particle->GetParticleName()]);
    SetHighEnergyLimit(highEnergyLimit[particle->GetParticleName()]);
// at least for HZE

  /*
  if (particle==protonDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if (particle==hydrogenDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[hydrogen]);
    SetHighEnergyLimit(highEnergyLimit[hydrogen]);
  }

  if (particle==heliumDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[helium]);
    SetHighEnergyLimit(highEnergyLimit[helium]);
  }

  if (particle==alphaPlusDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  if (particle==alphaPlusPlusDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlusPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlusPlus]);
  }

  if (particle==carbonDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[carbon]);
    SetHighEnergyLimit(highEnergyLimit[carbon]);
  }

  if (particle==oxygenDef) 
  {
    SetLowEnergyLimit(lowEnergyLimit[oxygen]);
    SetHighEnergyLimit(highEnergyLimit[oxygen]);
  }*/
  
//----------------------------------------------------------------------

  if( verboseLevel>0 ) 
  { 
    G4cout << "Rudd ionisation model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
           << particle->GetParticleName()
           << G4endl;
  }

  //
  
  if(!isInitialised) 
  {
    isInitialised = true;
  
    if(pParticleChange)
      fParticleChangeForGamma = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChangeForGamma = new G4ParticleChangeForGamma();
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARuddIonisationExtendedModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* particleDefinition,
					   G4double k,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNARuddIonisationExtendedModel" << G4endl;

 // Calculate total cross section for model

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != instance->GetIon("hydrogen")
      &&
      particleDefinition != instance->GetIon("alpha++")
      &&
      particleDefinition != instance->GetIon("alpha+")
      &&
      particleDefinition != instance->GetIon("helium")
      &&
      particleDefinition != instance->GetIon("carbon")
      &&
      particleDefinition != instance->GetIon("nitrogen")
      &&
      particleDefinition != instance->GetIon("oxygen")
      &&
      particleDefinition != instance->GetIon("iron")
     )
   	    
    return 0;

  G4double lowLim = 0;
  
  if (     particleDefinition == G4Proton::ProtonDefinition()
       ||  particleDefinition == instance->GetIon("hydrogen")
     )
       
       lowLim = lowEnergyLimitOfModelForA[1];
       
  else if (     particleDefinition == instance->GetIon("alpha++")
       ||       particleDefinition == instance->GetIon("alpha+")
       ||       particleDefinition == instance->GetIon("helium")
     )
       
       lowLim = lowEnergyLimitOfModelForA[4];

  else lowLim = lowEnergyLimitOfModelForA[5];
  
  G4double highLim = 0;
  G4double sigma=0;

  if (material->GetName() == "G4_WATER")
  {
    const G4String& particleName = particleDefinition->GetParticleName();
   
    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);

    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

    if (k <= highLim)
    {
      
      //SI : XS must not be zero otherwise sampling of secondaries method ignored

      if (k < lowLim) k = lowLim;

      //      
      
      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);
	
      if (pos != tableData.end())
      {
         G4DNACrossSectionDataSet* table = pos->second;
         if (table != 0)
         {
	      sigma = table->FindValue(k);
         }
      }
      else
      {
        G4Exception
        ("G4DNARuddIonisationExtendedModel::CrossSectionPerVolume: attempting to calculate cross section for wrong particle");
      }
  
    } // if (k >= lowLim && k < highLim)
    
    if (verboseLevel > 3)
    {
      G4cout << "---> Kinetic energy(eV)=" << k/eV << G4endl;
      G4cout << " - Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*densityWater/(1./cm) << G4endl;
    } 
 
  } // if (waterMaterial)
 
 return sigma*material->GetAtomicNumDensityVector()[1];	   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationExtendedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* /*couple*/,
					      const G4DynamicParticle* particle,
					      G4double,
					      G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNARuddIonisationExtendedModel" << G4endl;

  G4double lowLim = 0;
  G4double highLim = 0;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  // ZF: the following line summarizes the commented part 
  if(particle->GetDefinition()->GetAtomicMass() <= 4) lowLim = killBelowEnergyForA[particle->GetDefinition()->GetAtomicMass()];
 else lowLim = killBelowEnergyForA[5]*particle->GetDefinition()->GetAtomicMass();

/*if(particle->GetDefinition()->GetAtomicMass() >= 5) lowLim = killBelowEnergyForA[5]*particle->GetDefinition()->GetAtomicMass();


  if (     particle->GetDefinition() == G4Proton::ProtonDefinition()
       ||  particle->GetDefinition() == instance->GetIon("hydrogen")
     )
       
       lowLim = killBelowEnergyForA[1];
       
  if (     particle->GetDefinition() == instance->GetIon("alpha++")
       ||  particle->GetDefinition() == instance->GetIon("alpha+")
       ||  particle->GetDefinition() == instance->GetIon("helium")
     )
       
       lowLim = killBelowEnergyForA[4];*/


    
  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  // SI - the following is useless since lowLim is already defined 
  /* 
  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }
  */

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())highLim = pos2->second;

  if (k >= lowLim && k < highLim)
  {
      G4ParticleDefinition* definition = particle->GetDefinition();
      G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
      /*
      G4double particleMass = definition->GetPDGMass();
      G4double totalEnergy = k + particleMass;
      G4double pSquare = k*(totalEnergy+particleMass);
      G4double totalMomentum = std::sqrt(pSquare);
      */
      
      G4int ionizationShell = RandomSelect(k,particleName);

      G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);
  
      G4double bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

      G4double cosTheta = 0.;
      G4double phi = 0.; 
      RandomizeEjectedElectronDirection(definition, k,secondaryKinetic, cosTheta, phi, ionizationShell);
    
      G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
      G4double dirX = sinTheta*std::cos(phi);
      G4double dirY = sinTheta*std::sin(phi);
      G4double dirZ = cosTheta;
      G4ThreeVector deltaDirection(dirX,dirY,dirZ);
      deltaDirection.rotateUz(primaryDirection);
    
      // Ignored for ions on electrons
      /*
      G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));

      G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
      G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
      G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;

      G4ThreeVector direction;
      direction.set(finalPx,finalPy,finalPz);

      fParticleChangeForGamma->ProposeMomentumDirection(direction.unit()) ;
      */
      fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);

      fParticleChangeForGamma->SetProposedKineticEnergy(k-bindingEnergy-secondaryKinetic);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy);

      G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
      fvect->push_back(dp);

  }

  // SI - not useful since low energy of model is 0 eV
  
  if (k < lowLim) 
  {  
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
								    G4double k, 
								    G4int shell)
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  //-- Fast sampling method -----
  G4double proposed_energy;
  G4double random1;
  G4double value_sampling;
  G4double max1;

  do 
  {
       proposed_energy = ProposedSampledEnergy(particleDefinition, k, shell); // Proposed energy by inverse function sampling
       
       max1=0.;

       for(G4double en=0.; en<20.; en+=1.) if(RejectionFunction(particleDefinition, k, en, shell) > max1)
         max1=RejectionFunction(particleDefinition, k, en, shell);
       
       random1 = G4UniformRand()*max1;
       
       value_sampling = RejectionFunction(particleDefinition, k, proposed_energy, shell);
  
  } while(random1 > value_sampling);

  return(proposed_energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4DNARuddIonisationExtendedModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
								   G4double k, 
								   G4double secKinetic, 
								   G4double & cosTheta, 
								   G4double & phi,
									G4int shell )
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double maxSecKinetic = 0.;
  G4double maximumEnergyTransfer = 0.;
 
 /* if (particleDefinition == G4Proton::ProtonDefinition() 
      || particleDefinition == instance->GetIon("hydrogen")) 
  { 
      if(k/MeV < 100.)maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
      else {	
		G4double beta2 = 1.-(1.+k*);
		maxSecKinetic = 	
		}
  }
  
  if (particleDefinition == instance->GetIon("helium") 
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("alpha++"))
  { 
      maxSecKinetic = 4.* (0.511 / 3728) * k;
  }
  
  if (particleDefinition == G4Carbon::Carbon()) 
  { 
      maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k / 12.;
  }*/

// ZF. generalized & relativistic version 

  if( (k/MeV)/(particleDefinition->GetPDGMass()/MeV)  <= 0.1 )
  {
    maximumEnergyTransfer= 4.* (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
    maximumEnergyTransfer+=waterStructure.IonisationEnergy(shell);
  }
  else
  {
    G4double approx_nuc_number = particleDefinition->GetPDGMass() / proton_mass_c2;
    G4double en_per_nucleon = k/approx_nuc_number;
    G4double beta2 = 1. - 1./pow( (1.+(en_per_nucleon/electron_mass_c2)*(electron_mass_c2/proton_mass_c2)), 2.);
    G4double gamma = 1./sqrt(1.-beta2);
    maximumEnergyTransfer = 2.*electron_mass_c2*(gamma*gamma-1.)/(1.+2.*gamma*(electron_mass_c2/particleDefinition->GetPDGMass())+pow(electron_mass_c2/particleDefinition->GetPDGMass(), 2.) );
    maximumEnergyTransfer+=waterStructure.IonisationEnergy(shell);
  }

  maxSecKinetic = maximumEnergyTransfer-waterStructure.IonisationEnergy(shell);

  phi = twopi * G4UniformRand();
  
  if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
  else cosTheta = (2.*G4UniformRand())-1.;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::RejectionFunction(G4ParticleDefinition* particleDefinition, 
							      G4double k, 
							      G4double proposed_ws, 
							      G4int ionizationLevelIndex)
{
  const G4int j=ionizationLevelIndex;
  G4double Bj_energy, alphaConst;
  G4double Ry = 13.6*eV;
  const G4double Gj[5] = {0.99, 1.11, 1.11, 0.52, 1.};
  
  // const G4double Bj[5] = {12.61*eV, 14.73*eV, 18.55*eV, 32.20*eV, 539.7*eV}; //Ding Paper
  
  // Following values provided by M. Dingfelder (priv. comm)
  const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

  if (j == 4) 
  {
      alphaConst = 0.66;
      //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
      Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
      //---
  }
  else 
  {
      alphaConst = 0.64;
      Bj_energy = Bj[ionizationLevelIndex];
  }
 
  G4double energyTransfer = proposed_ws + Bj_energy;
  proposed_ws/=Bj_energy;
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4double tau = 0.;
  G4double A_ion = 0.;
  tau = (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
  A_ion = particleDefinition->GetAtomicMass();

  G4double v2;
  G4double beta2;

  if((tau/MeV)<5.447761194e-2)
  {
    v2 = tau / Bj_energy; 
    beta2 = 2.*tau / electron_mass_c2; 
  }
  // Relativistic
  else 		
  {	
    v2 = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
    beta2 =1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion)); 		
  }
  
  G4double v = std::sqrt(v2);
  G4double wc = 4.*v2 - 2.*v - (Ry/(4.*Bj_energy)); 
  G4double rejection_term = 1.+std::exp(alphaConst*(proposed_ws - wc) / v);
  rejection_term = (1./rejection_term)*CorrectionFactor(particleDefinition,k,ionizationLevelIndex) * Gj[j]; 
    //* (S/Bj_energy) ; Not needed anymore


  if (    particleDefinition == G4Proton::ProtonDefinition() 
	  || particleDefinition == instance->GetIon("hydrogen")
      ) 
  {
      return(rejection_term);
  }
  
  if(particleDefinition->GetAtomicMass() > 4) // anything above Helium
  {
      G4double Z = particleDefinition->GetAtomicNumber();
      G4double x = 100.*std::sqrt(beta2)/std::pow(Z,(2./3.));
      G4double Zeffion = Z*(1.-std::exp(-1.316*x+0.112*x*x-0.0650*x*x*x));
      rejection_term*=Zeffion*Zeffion;
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
// The following values are provided by M. Dingfelder (priv. comm)      
      slaterEffectiveCharge[1]=2.0;
      slaterEffectiveCharge[2]=2.0;
//
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

      G4double zEff = particleDefinition->GetPDGCharge() / eplus + particleDefinition->GetLeptonNumber();
  
      zEff -= ( sCoefficient[0] * S_1s(k, energyTransfer, slaterEffectiveCharge[0], 1.) +
		sCoefficient[1] * S_2s(k, energyTransfer, slaterEffectiveCharge[1], 2.) +
		sCoefficient[2] * S_2p(k, energyTransfer, slaterEffectiveCharge[2], 2.) );
	   
      rejection_term*= zEff * zEff;
  }  

  return (rejection_term);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4DNARuddIonisationExtendedModel::ProposedSampledEnergy(G4ParticleDefinition* particle, 
							      G4double k,  
							      G4int ionizationLevelIndex)
{

  const G4int j=ionizationLevelIndex;

  G4double A1, B1, C1, D1, E1, A2, B2, C2, D2, alphaConst ; 
  G4double Bj_energy;
  
  // const G4double Bj[5] = {12.61*eV, 14.73*eV, 18.55*eV, 32.20*eV, 539.7*eV}; //Ding Paper
  // Following values provided by M. Dingfelder (priv. comm)
  
  const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

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
      //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
      Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
      //---
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
      //B2 = 14.6; From Ding Paper
      // Value provided by M. Dingfelder (priv. comm)
      B2 = 11.6;
      //
      C2 = 0.60; 
      D2 = 0.04; 
      alphaConst = 0.64;

      Bj_energy = Bj[ionizationLevelIndex];
  }
  
  G4double tau = 0.;
  G4double A_ion = 0.;
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();
  tau = (electron_mass_c2 / particle->GetPDGMass()) * k;

  A_ion = particle->GetAtomicMass();

  G4double v2;
  G4double beta2;
  if((tau/MeV)<5.447761194e-2)
  {
    v2 = tau / Bj_energy; 
    beta2 = 2.*tau / electron_mass_c2; 
  }
  // Relativistic
  else
  {	
    v2 = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
    beta2 =1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion)); 		
  }
  
  G4double v = std::sqrt(v2);
  //G4double wc = 4.*v2 - 2.*v - (Ry/(4.*Bj_energy)); 
  G4double L1 = (C1* std::pow(v,(D1))) / (1.+ E1*std::pow(v, (D1+4.)));
  G4double L2 = C2*std::pow(v,(D2));
  G4double H1 = (A1*std::log(1.+v2)) / (v2+(B1/v2));
  G4double H2 = (A2/v2) + (B2/(v2*v2));
  G4double F1 = L1+H1;
  G4double F2 = (L2*H2)/(L2+H2);

// ZF. generalized & relativistic version 
  G4double maximumEnergy;

  //---- maximum kinetic energy , non relativistic ------
  if( (k/MeV)/(particle->GetPDGMass()/MeV)  <= 0.1 )
  {
    maximumEnergy = 4.* (electron_mass_c2 / particle->GetPDGMass()) * k; 
  }
  //---- relativistic -----------------------------------
  else
  {
    G4double gamma = 1./sqrt(1.-beta2);
    maximumEnergy = 2.*electron_mass_c2*(gamma*gamma-1.)/
       (1.+2.*gamma*(electron_mass_c2/particle->GetPDGMass())+pow(electron_mass_c2/particle->GetPDGMass(), 2.) );
  }
  
  //either it is transfered energy or secondary electron energy ...	
  //maximumEnergy-=Bj_energy;

  //-----------------------------------------------------
  G4double wmax = maximumEnergy/Bj_energy;
  G4double c = wmax*(F2*wmax+F1*(2.+wmax))/(2.*(1.+wmax)*(1.+wmax));
  c=1./c; //!!!!!!!!!!! manual calculus leads to  c=1/c
  G4double randVal = G4UniformRand();
  G4double proposed_ws = F1*F1*c*c + 2.*F2*c*randVal - 2.*F1*c*randVal;
  proposed_ws = -F1*c+2.*randVal+std::sqrt(proposed_ws);
//  proposed_ws = -F1*c+2.*randVal-std::sqrt(proposed_ws);
  proposed_ws/= ( F1*c + F2*c - 2.*randVal );
  proposed_ws*=Bj_energy;

  return(proposed_ws); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::S_1s(G4double t, 
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

G4double G4DNARuddIonisationExtendedModel::S_2s(G4double t,
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

G4double G4DNARuddIonisationExtendedModel::S_2p(G4double t, 
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

G4double G4DNARuddIonisationExtendedModel::R(G4double t,
				       G4double energyTransferred,
				       G4double slaterEffectiveChg,
				       G4double shellNumber) 
{
  // tElectron = m_electron / m_alpha * t
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double tElectron = 0.511/3728. * t;
  // The following values are provided by M. Dingfelder (priv. comm)    
  G4double H = 2.*13.60569172 * eV;
  G4double value = std::sqrt ( 2. * tElectron / H ) / ( energyTransferred / H ) *  (slaterEffectiveChg/shellNumber);
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::CorrectionFactor(G4ParticleDefinition* particleDefinition, G4double k, G4int shell) 
{
// ZF Shortened 
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

    if (particleDefinition == instance->GetIon("hydrogen") && shell < 4) 
    { 
	G4double value = (std::log10(k/eV)-4.2)/0.5;
	// The following values are provided by M. Dingfelder (priv. comm)    
        return((0.6/(1+std::exp(value))) + 0.9);
    }
    else 
    {    
	return(1.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationExtendedModel::RandomSelect(G4double k, const G4String& particle )
{   
  
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
 
  G4int level = 0;

  // Retrieve data table corresponding to the current particle type  

  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  pos = tableData.find(particle);

  if (pos != tableData.end())
  {
      G4DNACrossSectionDataSet* table = pos->second;

      if (table != 0)
      {
	  G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
	    
	  const size_t n(table->NumberOfComponents());
	  size_t i(n);
	  G4double value = 0.;
	    
	  while (i>0)
	  { 
	      i--;
	      valuesBuffer[i] = table->GetComponent(i)->FindValue(k);

	      value += valuesBuffer[i];
	  }
	    
	  value *= G4UniformRand();
	    
	  i = n;
	    
	  while (i > 0)
	  {
	      i--;
		
	      if (valuesBuffer[i] > value)
	      {
		  delete[] valuesBuffer;
		  return i;
	      }
	      value -= valuesBuffer[i];
	  }
	    
	  if (valuesBuffer) delete[] valuesBuffer;
	    
      }
  }
  else
  {
    G4Exception("G4DNARuddIonisationExtendedModel::RandomSelect: attempting to calculate cross section for wrong particle");
  }
     
  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::PartialCrossSection(const G4Track& track )
{
  G4double sigma = 0.;

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  
  G4double lowLim = 0;
  G4double highLim = 0;

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
      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);
	
      if (pos != tableData.end())
      {
	  G4DNACrossSectionDataSet* table = pos->second;
	  if (table != 0)
          {
	      sigma = table->FindValue(k);
          }
      }
      else
      {
	  G4Exception("G4DNARuddIonisationExtendedModel::PartialCrossSection: attempting to calculate cross section for wrong particle");
      }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationExtendedModel::Sum(G4double /* energy */, const G4String& /* particle */)
{
  return 0;
}

