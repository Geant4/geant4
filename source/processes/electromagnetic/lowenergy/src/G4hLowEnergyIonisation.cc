// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisation physics process -------
//                by Vladimir Ivanchenko, 14 July 1999 
//                was made on the base of G4hIonisation class
//                developed by Laszlo Urban
// ************************************************************
// It is the extention of the ionisation process for the slow 
// charged hadrons.
// ************************************************************
// 28 July   1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli added ICRU parametrisations for protons  
// 20 August 1999 G.Mancinelli added ICRU tables for alpha 
// 31 August 1999 V.Ivanchenko update and cleen up 
// 30 Sept.  1999 V.Ivanchenko minor upgrade 
// 12 Dec.   1999 S. Chauvie added Barkas correction 
// 19 Jan.   2000 V.Ivanchenko minor changing in Barkas corrections
// 02 April  2000 S. Chauvie linearization of barkas effect
// 03 April  2000 V.Ivanchenko Nuclear Stopping power for antiprotons
// 23 May    2000 MG Pia  Clean up for QAO model 
// 24 May    2000 MG Pia  Code properly indented to improve legibility
// 17 July   2000 V.Ivanchenko Bug in scaling AlongStepDoIt method
// 25 July   2000 V.Ivanchenko New design iteration
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include "G4hLowEnergyIonisation.hh"
#include "globals.hh"
#include "G4Poisson.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::G4hLowEnergyIonisation(const G4String& processName)
  : G4hLowEnergyLoss(processName),
    theMeanFreePathTable(NULL),
    protonLowEnergy(1.*keV),
    protonHighEnergy(2.*MeV),
    antiProtonLowEnergy(1.*keV),
    antiProtonHighEnergy(2.*MeV),
    theProtonTable("ICRU_R49p"),
    theAntiProtonTable("ICRU_R49p"),
    theNuclearTable("ICRU_R49Mollere"),
    theBetheBlochModel(NULL),
    theProtonModel(NULL),
    theAntiProtonModel(NULL),
    theNuclearStoppingModel(NULL),
    theIonEffChargeModel(NULL),
    nStopping(true),
    theBarkas(true),
    theProton (G4Proton::Proton()),
    theAntiProton (G4AntiProton::AntiProton()),
    theElectron ( G4Electron::Electron() ),
    factor(twopi_mc2_rcl2),
    protonMass(proton_mass_c2)
{ 
  SetPhysicsTableBining(10.0*eV, 100.0*TeV, 200) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::~G4hLowEnergyIonisation() 
{
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  if(theBetheBlochModel)delete theBetheBlochModel;
  if(theProtonModel)delete theProtonModel;
  if(theAntiProtonModel)delete theAntiProtonModel;
  if(theNuclearStoppingModel)delete theNuclearStoppingModel;
  if(theIonEffChargeModel)delete theIonEffChargeModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetPhysicsTableBining(G4double lowE, 
                                                   G4double highE,
						   G4int nBins)
{
  LowestKineticEnergy = lowE ;
  HighestKineticEnergy = highE ;
  TotBin = nBins ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetElectronicStoppingPowerModel(
                             const G4ParticleDefinition* aParticle,
                             const G4String& dedxTable)
  // This method defines the ionisation parametrisation method via its name 
{
  if(0 < aParticle->GetPDGCharge()) {
    SetProtonElectronicStoppingPowerModel(dedxTable) ;
  } else {
    SetAntiProtonElectronicStoppingPowerModel(dedxTable) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::InicialiseParametrisation() 

{
  // cuts for  electron 
  deltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  // Define models for parametrisation of electronic energy losses
  theBetheBlochModel = new G4hBetheBlochModel("Bethe-Bloch") ;
  theProtonModel = new G4hParametrisedLossModel(theProtonTable) ;
  theAntiProtonModel = new G4QAOLowEnergyLoss(theAntiProtonTable) ;
  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  theIonEffChargeModel = new G4hIonEffChargeSquare("Ziegler1988") ;

  // Energy limits for parametrisation of electronic energy losses
  protonLowEnergy = G4std::max(protonLowEnergy,
                               theProtonModel->LowEnergyLimit(theProton)) ;   

  G4double x1 = theBetheBlochModel->LowEnergyLimit(theProton) ;   
  G4double x2 = theProtonModel->HighEnergyLimit(theProton) ;   
  if(protonHighEnergy < x1) protonHighEnergy = x1 ;
  if(protonHighEnergy > x2) protonHighEnergy = x2 ;
  
  antiProtonLowEnergy = G4std::max(antiProtonLowEnergy,
                      theAntiProtonModel->LowEnergyLimit(theAntiProton)) ;   

  x2 = theAntiProtonModel->HighEnergyLimit(theAntiProton) ;   
  if(antiProtonHighEnergy < x1) antiProtonHighEnergy = x1 ;
  if(antiProtonHighEnergy > x2) antiProtonHighEnergy = x2 ;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildPhysicsTable(
                       const G4ParticleDefinition& aParticleType) 
  
  //  just call BuildLossTable+BuildLambdaTable
{
  InicialiseParametrisation() ;

  G4double charge = aParticleType.GetPDGCharge()/eplus ;

  G4double electronCutInRange = theElectron->GetCuts(); 

  if(0.0 < charge) 
    {

      if( (ptableElectronCutInRange != electronCutInRange)  
	  || (theDEDXpTable == NULL))
	{
	  BuildLossTable(*theProton) ;
	  RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
	  CounterOfpProcess++;
	}

    } else{

      if( (pbartableElectronCutInRange != electronCutInRange)  
	  || (theDEDXpbarTable == NULL))
	{
	  BuildLossTable(*theAntiProton) ;
	  RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
	  CounterOfpbarProcess++;
	}
    }
  
  BuildLambdaTable(aParticleType) ;

  BuildDEDXTable(aParticleType) ;

  if((&aParticleType == theProton) ) PrintInfoDefinition() ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLossTable(
                             const G4ParticleDefinition& aParticleType) 
{
  // Tables for different hadrons will be different because of
  // small difference in Tmax connected with RateMass !?
  
  // Inicialisation
  G4double lowEdgeEnergy , ionloss, ionlossBB, paramB ;
  G4double lowEnergy, highEnergy, charge;

  if(aParticleType == *theProton) {
    lowEnergy = protonLowEnergy ;
    highEnergy = protonHighEnergy ;
    charge = 1.0 ;
  } else {
    lowEnergy = antiProtonLowEnergy ;
    highEnergy = antiProtonHighEnergy ;
    charge = -1.0 ;
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  
  //  create table
  G4int numOfMaterials = theMaterialTable->length();
  
  if ( theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }

  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  //  loop for materials  
  for (G4int j=0; j<numOfMaterials; j++) {

    // create physics vector and fill it  
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                         HighestKineticEnergy, 
                                                         TotBin);
      
    // get material parameters needed for the energy loss calculation  
    const G4Material* material= (*theMaterialTable)[j];

    // low energy of Bethe-Bloch formula for this material
    G4double highE = G4std::max(highEnergy,theBetheBlochModel->
                                LowEnergyLimit(&aParticleType,material)) ;

    if ( charge > 0.0 ) {
      ionloss = ProtonParametrisedDEDX(material,highE) ;
    } else {
      ionloss = AntiProtonParametrisedDEDX(material,highE) ;
    }

    ionlossBB = theBetheBlochModel->TheValue(&aParticleType,material,highE) ;
    ionlossBB -= DeltaRaysEnergy(material,highE,protonMass) ;

    if(theBarkas) {
      ionlossBB += BarkasTerm(material,highE)*charge ;
      ionlossBB += BlochTerm(material,highE,protonMass,1.0) ;
    }
    paramB =  ionloss/ionlossBB - 1.0 ; 

    // now comes the loop for the kinetic energy values      
    for (G4int i = 0 ; i < TotBin ; i++) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;

      // low energy part for this material, parametrised energy loss formulae
      if ( lowEdgeEnergy < highE ) {
	 
        if ( charge > 0.0 ) {
          ionloss = ProtonParametrisedDEDX(material,lowEdgeEnergy) ;
	} else {
          ionloss = AntiProtonParametrisedDEDX(material,lowEdgeEnergy) ;
	}

      } else {
	    
        // high energy part for this material, Bethe-Bloch formula
        ionloss = theBetheBlochModel->TheValue(theProton,material, 
                                               lowEdgeEnergy) ; 

        ionloss -= DeltaRaysEnergy(material,lowEdgeEnergy,protonMass) ;

        if(theBarkas) {
          ionloss += BarkasTerm(material,lowEdgeEnergy)*charge ;
          ionloss += BlochTerm(material,lowEdgeEnergy,protonMass,1.0) ;
        }
	ionloss *= (1.0 + paramB*highEnergy/lowEdgeEnergy) ;
      }
	  
      // now put the loss into the vector
      aVector->PutValue(i,ionloss) ;
    }
    // Insert vector for this material into the table
    theLossTable->insert(aVector) ;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLambdaTable(
                       const G4ParticleDefinition& aParticleType) 
  
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 
  
  G4double lowEdgeEnergy , value ,sigma ;
  G4bool isOutRange ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    
  //create table
  G4int numOfMaterials = theMaterialTable->length();
  
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  
  theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);
  
  // get electron and particle cuts in kinetic energy
  
  G4double* deltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
  
  // loop for materials 
  
  for (G4int J=0 ; J < numOfMaterials; J++) { 

    //create physics vector then fill it ....      
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                         HighestKineticEnergy, 
                                                         TotBin);
      
    // compute the (macroscopic) cross section first
      
    const G4Material* material = (*theMaterialTable)[J];   
    const G4ElementVector* theElementVector = 
                           material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector = 
                           material->GetAtomicNumDensityVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;
    G4double excEnergy = material->GetIonisation()->GetMeanExcitationEnergy();
      
      // get the electron kinetic energy cut for the actual material,
      //  it will be used in ComputeMicroscopicCrossSection
      // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
      //   ------------------------------------------------------
      
    // Cut in Delta energy is limited by exitation energy
    G4double deltaCut = G4std::max(excEnergy,deltaCutInKineticEnergy[J]) ;
      
    for ( G4int i = 0 ; i < TotBin ; i++ ) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      G4double chargeSquare = theIonEffChargeModel->TheValue(&aParticleType,
                              material,lowEdgeEnergy) ;
	  
      G4double sigma = 0.0 ;  
      for (G4int iel=0; iel<NumberOfElements; iel++ ) {
	sigma += theAtomicNumDensityVector[iel]*chargeSquare*
		 ComputeMicroscopicCrossSection(aParticleType,
			lowEdgeEnergy,(*theElementVector)(iel)->GetZ(),
                        excEnergy,deltaCut ) ;
      }
	  
      // mean free path = 1./macroscopic cross section
	  
      value = sigma<=0 ? DBL_MAX : 1./sigma ;     
	  
      aVector->PutValue(i, value) ;
    }
      
    theMeanFreePathTable->insert(aVector);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ComputeMicroscopicCrossSection(
                           const G4ParticleDefinition& aParticleType,
			         G4double kineticEnergy,
				 G4double atomicNumber,
				 G4double excEnergy,
				 G4double deltaCutInEnergy) const
{
  //******************************************************************
  // cross section formula is OK for spin=0, 1/2, 1 only !
  // *****************************************************************
  
  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = z )
  
  G4double energy, beta2, tmax, var ;
  G4double totalCrossSection = 0.0 ;
  
  G4double eexc2 = excEnergy*excEnergy ;
  G4double particleMass = aParticleType.GetPDGMass() ;
  
  // get particle data ...................................
  
  energy = kineticEnergy + particleMass;
  
  // some kinematics......................
  
  beta2  = kineticEnergy*(energy+particleMass) / (energy*energy);
  var    = particleMass+electron_mass_c2;
  tmax   = 2.0*electron_mass_c2*kineticEnergy * (energy+particleMass)
         / (var*var + 2.0*electron_mass_c2*kineticEnergy) ;
  
  // now you can calculate the total cross section 
  
  if( tmax > deltaCutInEnergy ) {
    
    var=deltaCutInEnergy/tmax;
    totalCrossSection = (1.0 - var*(1.0 - beta2*log(var))) / deltaCutInEnergy ;
    G4double spin = aParticleType.GetPDGSpin() ;
    
    // +term for spin=1/2 particle
    if( 0.5 == spin )
      totalCrossSection +=  0.5 * (tmax - deltaCutInEnergy) / (energy*energy);
    
    // +term for spin=1 particle
    else if( 0.9 < spin )
      totalCrossSection += -log(var)/(3.0*deltaCutInEnergy) +
	(tmax - deltaCutInEnergy) * ( (5.0+ 1.0/var)*0.25 / (energy*energy) - 
	 beta2 / (tmax * deltaCutInEnergy) ) / 3.0 ;
    
    totalCrossSection *= twopi_mc2_rcl2 * atomicNumber / beta2 ;
  }
  
  return totalCrossSection ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetConstraints(
                                 const G4DynamicParticle* particle,
				 const G4Material* material) 
{
  // returns the Step limit
  // dEdx is calculated as well as the range  
  // based on Effective Charge Approach
  
  G4double stepLimit ;
  G4bool isOut ;
  
  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double kineticEnergy = particle->GetKineticEnergy() ;
  G4double charge = particle->GetCharge() ;
  
  // Scale the kinetic energy
  
  G4double tscaled= kineticEnergy*massRatio ; 
  G4double chargeSquare = theIonEffChargeModel->TheValue(particle,material) ;
  G4double dx, s, highEnergy;
  
  if(charge > 0.0) {
    
    fdEdx = G4EnergyLossTables::GetDEDX(theProton, tscaled, material) 
          * chargeSquare ;
    
    fRangeNow = G4EnergyLossTables::GetRange(theProton, tscaled, material) ;
    highEnergy = protonHighEnergy ;
    dx = G4EnergyLossTables::GetRange(theProton, highEnergy, material) ;
    
    if(tscaled < highEnergy) {
      // For Bragg's peak the limit in range is estimated 
      // in order to be inside linLossLimit on each step
      fdEdx = ProtonParametrisedDEDX(material, tscaled) * chargeSquare ;
    }
    
    // Antiprotons and negative hadrons 
  } else {
    
    fdEdx = G4EnergyLossTables::GetDEDX(theAntiProton, tscaled, material) 
          * chargeSquare ;
    
    fRangeNow = G4EnergyLossTables::GetRange(theAntiProton, tscaled, material);
    highEnergy = antiProtonHighEnergy ;
    dx = G4EnergyLossTables::GetRange(theAntiProton, highEnergy, material) ;
    
    if(tscaled < highEnergy) {
      // For Bragg's peak the limit in range is estimated 
      // in order to be inside linLossLimit on each step
      fdEdx = AntiProtonParametrisedDEDX(material, tscaled) * chargeSquare ;
    }
  }
  
  //
  fRangeNow /= (chargeSquare*massRatio) ;
  dx        /= (chargeSquare*massRatio) ;
  stepLimit  = fRangeNow ;
  
  // compute the (random) Step limit in standard energy range
  if(tscaled > highEnergy ) {
    if(fRangeNow > finalRange) {
      stepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;
      
      //  randomise this value
      if(rndmStepFlag) stepLimit = finalRange + 
                                  (stepLimit-finalRange)*G4UniformRand() ;
    }
    // cut step in front of Bragg's peak
    if(stepLimit > fRangeNow - dx*0.9) stepLimit = fRangeNow - dx*0.9 ;

  // Step limit in low energy range
  } else {
    stepLimit = G4std::min (fRangeNow, dx * linLossLimit) ;    
  }
  
  return stepLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::AlongStepDoIt(
                                           const G4Track& trackData, 
                                           const G4Step& stepData) 
{
  // compute the energy loss after a step 
  
  G4double finalT = 0.0 ;
  
  aParticleChange.Initialize(trackData) ;
  G4Material* material = trackData.GetMaterial() ;
  
  // get the actual (true) Step length from stepData 
  const G4double step = stepData.GetStepLength() ;
  
  const G4DynamicParticle* particle = trackData.GetDynamicParticle() ;
  
  G4int index = material->GetIndex() ;
  G4double kineticEnergy = particle->GetKineticEnergy() ;
  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double charge    = (particle->GetCharge())/eplus ;
  G4double tscaled= kineticEnergy*massRatio ; 
  G4double chargeSquare = charge*charge ;
  G4double eloss = 0.0 ;
  G4double nloss = 0.0 ;
  
    // very small particle energy
  if(kineticEnergy < MinKineticEnergy) {
    eloss = kineticEnergy ;  

    // particle energy outside tabulated energy range
  } else if(( kineticEnergy > HighestKineticEnergy) || 
            ( kineticEnergy <= LowestKineticEnergy)) {
    eloss = step*fdEdx ; 
 
    // proton parametrisation model  
  } else if(tscaled < protonHighEnergy && charge > 0.0) {
    
    if(nStopping) {
      nloss = theNuclearStoppingModel->TheValue(particle, material) ;
    }
    
    G4double eFinal = kineticEnergy - step*(fdEdx + nloss) ;
    
    if(0.0 < eFinal) {
      eloss = (fdEdx + 
               ProtonParametrisedDEDX(material,eFinal*massRatio)*chargeSquare)
            *  step * 0.5 ;
      if(nStopping) {
	nloss = ( nloss + 
                 (theNuclearStoppingModel->TheValue(particle, material)))
              *   step*0.5 ;
      }
    } else {
      eloss = kineticEnergy ;
    }

    // antiproton parametrisation model
  } else if(tscaled < antiProtonHighEnergy && charge < 0.0) {
    
    if(nStopping) {
      nloss = theNuclearStoppingModel->TheValue(particle, material) ;
    }
    
    G4double eFinal = kineticEnergy - step*(fdEdx + nloss) ;
    
    if(0.0 < eFinal) {
      eloss = (fdEdx + 
            AntiProtonParametrisedDEDX(material,eFinal*massRatio)*chargeSquare)
            *  step * 0.5 ;
      if(nStopping) {
	nloss = ( nloss + 
                 (theNuclearStoppingModel->TheValue(particle, material)))
              *   step*0.5 ;
      }
    } else {
      eloss = kineticEnergy ;
    }

    // big step
  } else if(step >= fRangeNow ) {
    eloss = kineticEnergy ;

    // tabulated range 
  } else {
    
    // step longer than linear step limit 
    if(step > linLossLimit*fRangeNow) {
      
      G4double rscaled= fRangeNow*massRatio*chargeSquare ;
      G4double sscaled=   step   *massRatio*chargeSquare ;
      
      if(charge>0.) {
        eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theProton,rscaled, material) -
	        G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theProton,rscaled-sscaled,material) ;
       
      } else {
        eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theAntiProton,rscaled,material) -
	        G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theAntiProton,rscaled-sscaled,material) ;
      }
      eloss /= massRatio ;

    // step shorter than linear step limit       
    } else {
      eloss = step*fdEdx ;
    }
  }
  
  finalT = kineticEnergy - eloss - nloss;
  
  if(finalT > MinKineticEnergy) {
    
    //  now the electron loss with fluctuation
    if((EnlossFlucFlag) && (finalT < kineticEnergy)  
                        && (kineticEnergy > LowestKineticEnergy)) {
      
      eloss = GetLossWithFluct(particle, material, eloss/chargeSquare)
	    * chargeSquare ;
      if(nStopping) {
        nloss = NuclearLossFluctuation(particle, material, nloss) ;
      }
      finalT = kineticEnergy - eloss - nloss ;
    }    
  }
  
  //  kill the particle if the kinetic energy <= 0  
  if (finalT <= 0.0 )
    {
      finalT = 0.0 ;
      if( "proton" == (particle->GetDefinition()->GetParticleName()) )
	aParticleChange.SetStatusChange(fStopAndKill);
      else  
	aParticleChange.SetStatusChange(fStopButAlive); 
    } 
  
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(kineticEnergy-finalT) ;
  
  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ProtonParametrisedDEDX(
                                 const G4Material* material, 
				       G4double kineticEnergy) const
{
  G4double eloss = 0.0;

    // Free Electron Gas Model  
  if(kineticEnergy < protonLowEnergy) {
    eloss = (theProtonModel->TheValue(theProton, material, protonLowEnergy))
          * sqrt(kineticEnergy/protonLowEnergy) ;
    
    // Parametrisation
  } else {
    eloss = theProtonModel->TheValue(theProton, material, kineticEnergy) ;
  }
  
  // Delta rays energy
  eloss -= DeltaRaysEnergy(material,kineticEnergy,protonMass) ;

  if(eloss < 0.0) eloss = 0.0 ;  
  
  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::AntiProtonParametrisedDEDX(
                                 const G4Material* material, 
				       G4double kineticEnergy) const
{
  G4double eloss = 0.0 ;

  // Choose the model
  G4VLowEnergyModel * theModel ;
  if(theAntiProtonModel->IsInCharge(theAntiProton,material)) {
    theModel = theAntiProtonModel ;
  } else { 
    theModel = theProtonModel ;
  }

    // Free Electron Gas Model  
  if(kineticEnergy < antiProtonLowEnergy) {
    eloss = theModel->TheValue(theAntiProton, material, antiProtonLowEnergy);
    
    // Parametrisation
  } else {
    eloss = theModel->TheValue(theAntiProton, material, kineticEnergy) ;
  }

    // Proton model is used
  if(theBarkas && (theModel == theProtonModel)) 
                  eloss -= 2.0*BarkasTerm(material, kineticEnergy);
  
  // Delta rays energy
  eloss -= DeltaRaysEnergy(material,kineticEnergy,protonMass) ;

  if(eloss < 0.0) eloss = 0.0 ;  
  
  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::DeltaRaysEnergy(
                                 const G4Material* material, 
                                       G4double kineticEnergy,
                                       G4double particleMass) const
{
  G4double dloss = 0.0 ;

  G4double deltaCutNow = deltaCutInKineticEnergy[(material->GetIndex())] ;   
  G4double electronDensity = material->GetElectronDensity();
  G4double eexc = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc ;

  G4double tau = kineticEnergy/particleMass ;    
  G4double rateMass = electron_mass_c2/particleMass ;
  
  // some local variables 
  
  G4double gamma,bg2,beta2,tmax,x ;
  
  gamma = tau + 1.0 ;
  bg2 = tau*(tau+2.0) ;
  beta2 = bg2/(gamma*gamma) ;
  tmax = 2.*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

  // Validity range for delta electron cross section
  G4double deltaCut = G4std::max(deltaCutNow, eexc);
    
  if ( deltaCut < tmax) {
    x = deltaCut / tmax ;
    dloss = ( beta2 * (x - 1.0) - log(x) ) * factor * electronDensity / beta2 ;
  }
  return dloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::PostStepDoIt(
                                           const G4Track& trackData,   
					   const G4Step& stepData) 
{
  // Units are expressed in GEANT4 internal units.
  
  G4double KineticEnergy,TotalEnergy,TotalMomentum,
           betasquare,MaxKineticEnergyTransfer,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,te2,grej,Psquare,Esquare,summass,rate,grejc,finalMomentum ;
  
  aParticleChange.Initialize(trackData) ;
  G4Material* aMaterial = trackData.GetMaterial() ;
  G4double Eexc = aMaterial->GetIonisation()->GetMeanExcitationEnergy();
  G4double Eexc2 = Eexc*Eexc ;
  
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;
  
  ParticleMass=aParticle->GetDefinition()->GetPDGMass();
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy ;
  summass = ParticleMass + electron_mass_c2 ;    
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;
  
  //  get kinetic energy cut for the electron....
  G4double DeltaCutInKineticEnergyNow = 
           deltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics......................
  
  betasquare=Psquare/Esquare ;
  MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
    /(summass*summass+2.*electron_mass_c2*KineticEnergy);

  // Validity range for delta electron cross section
  G4double DeltaCut = G4std::max(DeltaCutInKineticEnergyNow,Eexc);
    
  // sampling kinetic energy of the delta ray 
  
  if( MaxKineticEnergyTransfer <= DeltaCut )
    {
      // pathological case (it should not happen ,
      // there is no change at all).....
      
      // return &aParticleChange;
      return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
    }
  else
    {
      // normal case ......................................
      xc   = DeltaCut / MaxKineticEnergyTransfer ;
      rate = MaxKineticEnergyTransfer / TotalEnergy ;
      rate = rate*rate ;
      G4double spin = aParticle->GetDefinition()->GetPDGSpin() ;   
      
      // sampling follows ...      
      do {
	x=xc/(1.-(1.-xc)*G4UniformRand());

        if(0.0 == spin) {
          grej = 1.0 - betasquare * x ;

        } else if (0.5 == spin) {
          grej = (1.0 - betasquare * x + 0.5*x*x*rate) / (1.0 + 0.5 * rate) ;

        } else {
          grej = (1.0 - betasquare * x ) * (1.0 + x/ (3.0*xc)) + 
	    x * x * rate * (1.0 + 0.5 * x / xc) / 3.0 / 
	    (1.0 + 1.0/(3.0*xc) + rate *(1.0+ 0.5/xc) /3.0) ;
        }

      } while( G4UniformRand() > grej );
    }
  
  DeltaKineticEnergy = x * MaxKineticEnergyTransfer ;
  
  if(DeltaKineticEnergy <= 0.)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  
  DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
						  2. * electron_mass_c2 )) ;
  TotalMomentum = sqrt(Psquare) ;
  costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
    /(DeltaTotalMomentum * TotalMomentum) ;
  
  //  protection against costheta > 1 or < -1   ---------------
  if ( costheta < -1. ) 
    costheta = -1. ;
  if ( costheta > +1. ) 
    costheta = +1. ;
  
  //  direction of the delta electron  ........
  phi = twopi * G4UniformRand() ; 
  sintheta = sqrt((1.+costheta)*(1.-costheta));
  dirx = sintheta * cos(phi) ;
  diry = sintheta * sin(phi) ;
  dirz = costheta ;
  
  G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
  DeltaDirection.rotateUz(ParticleDirection) ;
  
  // create G4DynamicParticle object for delta ray
  G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
  theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
  theDeltaRay->SetMomentumDirection(DeltaDirection.x(),
                                    DeltaDirection.y(),
                                    DeltaDirection.z()); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
  
  // fill aParticleChange 
  finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ;
  G4double Edep = 0 ;
  
  if (finalKineticEnergy > MinKineticEnergy)
    {
      finalPx = TotalMomentum*ParticleDirection.x()
	- DeltaTotalMomentum*DeltaDirection.x();
      finalPy = TotalMomentum*ParticleDirection.y()
	- DeltaTotalMomentum*DeltaDirection.y();
      finalPz = TotalMomentum*ParticleDirection.z()
	- DeltaTotalMomentum*DeltaDirection.z();
      finalMomentum =
	sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
      finalPx /= finalMomentum ;
      finalPy /= finalMomentum ;
      finalPz /= finalMomentum ;
      
      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    }
  else
    {
      finalKineticEnergy = 0. ;
      Edep = finalKineticEnergy ;
      if (aParticle->GetDefinition()->GetParticleName() == "proton")
	aParticleChange.SetStatusChange(fStopAndKill);
      else  aParticleChange.SetStatusChange(fStopButAlive);
    }
  
  aParticleChange.SetEnergyChange( finalKineticEnergy );
  aParticleChange.SetNumberOfSecondaries(1);   
  aParticleChange.AddSecondary( theDeltaRay );
  aParticleChange.SetLocalEnergyDeposit (Edep);
  
  //ResetNumberOfInteractionLengthLeft();
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ComputeDEDX(
                                 const G4ParticleDefinition* aParticle,
                                 const G4Material* material, 
				       G4double kineticEnergy) 
{  
  G4double dedx = 0.0 ;
    
  G4double tscaled = kineticEnergy*protonMass/(aParticle->GetPDGMass()) ; 
  G4double charge  = aParticle->GetPDGCharge() ;
  
  if(charge>0.0) {
    if(tscaled > protonHighEnergy) {
      dedx=G4EnergyLossTables::GetPreciseDEDX(theProton,tscaled,material) ;

    } else {
      dedx=ProtonParametrisedDEDX(material,tscaled) ;
    }

  } else {
    if(tscaled > antiProtonHighEnergy) {
      dedx=G4EnergyLossTables::GetPreciseDEDX(theAntiProton,tscaled,material); 
      
    } else {
      dedx=AntiProtonParametrisedDEDX(material,tscaled) ;
    }
  }    
  dedx *= theIonEffChargeModel->TheValue(aParticle, material, kineticEnergy) ;
  
  return dedx ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::BarkasTerm(const G4Material* material,
  				                  G4double kineticEnergy) const
//Function to compute the Barkas term from:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397 
//
{
  static double FTable[47][2] = { 
    0.02, 21.5,	0.03, 20.0,	0.04, 18.0,	0.05, 15.6,
    0.06, 15.0,	0.07, 14.0,   	0.08, 13.5,   	0.09, 13,
    0.1,  12.2,  	0.2,   9.25,  	0.3,   7.0,     0.4,   6.0,   
    0.5,  4.5, 	0.6,   3.5,   	0.7,   3.0,     0.8,   2.5,   
    0.9,  2.0,  	1.0,   1.7,  	1.2,   1.2,   	1.3,   1.0,     
    1.4,  0.86,  	1.5,   0.7,	1.6,   0.61,  	1.7,   0.52,  
    1.8,  0.5,   	1.9,   0.43,	2.0,   0.42,  	2.1,   0.3,   
    2.4,  0.2,	3.0,   0.13,  	3.08,  0.1,   	3.1,   0.09, 
    3.3,  0.08,     3.5,   0.07,  	3.8,   0.06,	4.0,   0.051, 
    4.1,  0.04,  	4.8,   0.03,    5.0,   0.024,   5.1,   0.02,
    6.0,  0.013,    6.5,   0.01,	7.0,   0.009, 	7.1,   0.008,
    8.0,  0.006, 	9.0,   0.0032, 10.0,   0.0025 };

  // Information on particle and material
  G4double beta2 = 1.0 + kineticEnergy / protonMass ;
  beta2 = 1.0 - beta2*beta2 ;
  if(0.0 >= beta2) return 0.0;
  
  G4double BarkasTerm = 0.0;
  G4double AMaterial = 0.0;
  G4double ZMaterial = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int numberOfElements = material->GetNumberOfElements();
  
  for (G4int i = 0; i<numberOfElements; ++i) {

    AMaterial = (*theElementVector)(i)->GetA()*mole/g;
    ZMaterial = (*theElementVector)(i)->GetZ();
    
    G4double X = 137.0 * 137.0 * beta2 / ZMaterial;
  
    // Variables to compute L_1
    G4double Eta0Chi = 0.8;
    G4double EtaChi = Eta0Chi * ( 1.0 + 6.02*pow( ZMaterial,-1.19 ) );
    G4double W = ( EtaChi * pow( ZMaterial,1.0/6.0 ) ) / sqrt(X); 
    G4double FunctionOfW = FTable[46][1]*FTable[46][0]/W ;
    
    for(int IndexOfFTable=0; IndexOfFTable<47; IndexOfFTable++) {
    
      if(W<FTable[IndexOfFTable][0] ) {
    
        if(0 == IndexOfFTable) {
     	  FunctionOfW = FTable[0][1] ;
        }
    
        else {
          FunctionOfW = ( W - FTable[IndexOfFTable-1][0]) 
	    / (FTable[IndexOfFTable][0] - FTable[IndexOfFTable-1][0]) ;    	
	  FunctionOfW *=(FTable[IndexOfFTable][1]-FTable[IndexOfFTable-1][1]) ;
	  FunctionOfW += FTable[IndexOfFTable-1][1] ;
	}
    
        break;
      }
    
    }
    
    BarkasTerm += FunctionOfW /( sqrt(ZMaterial * X) * X);
  }

  BarkasTerm *= factor * (material->GetElectronDensity()) / beta2 ;

  return BarkasTerm;
}
        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::BlochTerm(const G4Material* material,
                                                 G4double kineticEnergy,
                                                 G4double particleMass,
  		        		         G4double chargeSquare) const
//Function to compute the Bloch term from:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397 
//
{
  G4double eloss = 0.0 ;
  G4double beta2 = 1.0 + kineticEnergy / particleMass ;
  beta2 = 1.0 - beta2*beta2 ;
  G4double y = chargeSquare / (137.0*137.0*beta2) ;

  if(y < 0.05) {
    eloss = -1.202 ;
     
  } else {
    eloss = 1.0 / (1.0 + y) ;
    G4double de = eloss ;

    for(G4int i=2; de>eloss*0.01; i++) {
      de = 1.0/( i * (i*i + y)) ;
      eloss += de ;
    }
  }
  eloss *= y * factor * (material->GetElectronDensity()) / beta2 ;
 
  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetLossWithFluct(
                                 const G4DynamicParticle* particle,
                                 const G4Material* material,
                                       G4double    MeanLoss) const
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same as in Glandz in Geant3.
{
   static const G4double minLoss = 1.*eV ;
   static const G4double probLim = 0.01 ;
   static const G4double sumaLim = -log(probLim) ;
   static const G4double alim = 10.;
   static const G4double kappa = 10. ;


// check if the material has changed ( cache mechanism) - is not used 

//  if (material != lastMaterial)
//    {
//      lastMaterial = material;

  G4int    imat        = material->GetIndex(); 
  G4double excEnergy   = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double f1Fluct     = material->GetIonisation()->GetF1fluct();
  G4double f2Fluct     = material->GetIonisation()->GetF2fluct();
  G4double e1Fluct     = material->GetIonisation()->GetEnergy1fluct();
  G4double e2Fluct     = material->GetIonisation()->GetEnergy2fluct();
  G4double e1LogFluct  = material->GetIonisation()->GetLogEnergy1fluct();
  G4double e2LogFluct  = material->GetIonisation()->GetLogEnergy2fluct();
  G4double rateFluct   = material->GetIonisation()->GetRateionexcfluct();
  G4double ipotFluct   = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double ipotLogFluct= material->GetIonisation()->GetLogMeanExcEnergy();
//    }
  G4double threshold,w1,w2,C,
           beta2,suma,e0,loss,lossc ,w;
  G4double a1,a2,a3;
  G4int p1,p2,p3;
  G4int nb;
  G4double Corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp1,dp3;
  G4double siga ;

  // shortcut for very very small loss 
  if(MeanLoss < minLoss) return MeanLoss ;

  // get particle data
  G4double Tkin   = particle->GetKineticEnergy();
  ParticleMass = particle->GetMass() ;
  G4double DeltaCutInKineticEnergyNow = 
           deltaCutInKineticEnergy[material->GetIndex()];

  // Validity range for delta electron cross section
  threshold = G4std::max(DeltaCutInKineticEnergyNow,excEnergy);

  G4double rmass = electron_mass_c2/ParticleMass;
  G4double tau   = Tkin/ParticleMass, tau1 = tau+1., tau2 = tau*(tau+2.);
  G4double Tm    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);

  if (Tm <= ipotFluct) Tm = ipotFluct ;
  
  if(Tm > threshold) Tm = threshold;
  beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation ?
  if(MeanLoss >= kappa*Tm)
  {
    siga = sqrt(MeanLoss*Tm*(1.-0.5*beta2)) ;
    loss = RandGauss::shoot(MeanLoss,siga) ;
    if(loss < 0.) loss = 0. ;
    return loss ;
  }

  w1 = Tm/ipotFluct;
  w2 = log(2.*electron_mass_c2*tau2);

  C = MeanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  if(Tm > ipotFluct)
    a3 = rateFluct*MeanLoss*(Tm-ipotFluct)/(ipotFluct*Tm*log(w1));
  else
  {
     a1 /= 1.-rateFluct ;
     a2 /= 1.-rateFluct ;
     a3  = 0. ;
  } 

  suma = a1+a2+a3;

  loss = 0. ;

  if(suma < sumaLim)             // very small Step
    {
      e0 = material->GetIonisation()->GetEnergy0fluct();

      if(Tm == ipotFluct)
      {
        a3 = MeanLoss/e0;

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,int(RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        loss = p3*e0 ;

        if(p3 > 0)
          loss += (1.-2.*G4UniformRand())*e0 ;

      }
      else
      {
        Tm = Tm-ipotFluct+e0 ;
        a3 = MeanLoss*(Tm-e0)/(Tm*e0*log(Tm/e0));

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,int(RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        if(p3 > 0)
        {
          w = (Tm-e0)/Tm ;
          if(p3 > nmaxCont2)
          {
            dp3 = G4float(p3) ;
            Corrfac = dp3/G4float(nmaxCont2) ;
            p3 = nmaxCont2 ;
          }
          else
            Corrfac = 1. ;

          for(G4int i=0; i<p3; i++) loss += 1./(1.-w*G4UniformRand()) ;
          loss *= e0*Corrfac ;  
        }        
      }
    }
    
  else                              // not so small Step
    {
      // excitation type 1
      if(a1>alim)
      {
        siga=sqrt(a1) ;
        p1 = G4std::max(0,int(RandGauss::shoot(a1,siga)+0.5));
      }
      else
       p1 = G4Poisson(a1);

      // excitation type 2
      if(a2>alim)
      {
        siga=sqrt(a2) ;
        p2 = G4std::max(0,int(RandGauss::shoot(a2,siga)+0.5));
      }
      else
        p2 = G4Poisson(a2);

      loss = p1*e1Fluct+p2*e2Fluct;
 
      // smearing to avoid unphysical peaks
      if(p2 > 0)
        loss += (1.-2.*G4UniformRand())*e2Fluct;   
      else if (loss>0.)
        loss += (1.-2.*G4UniformRand())*e1Fluct;   

      // ionisation .......................................
     if(a3 > 0.)
     {
      if(a3>alim)
      {
        siga=sqrt(a3) ;
        p3 = G4std::max(0,int(RandGauss::shoot(a3,siga)+0.5));
      }
      else
        p3 = G4Poisson(a3);

      lossc = 0.;
      if(p3 > 0)
      {
        na = 0.; 
        alfa = 1.;
        if (p3 > nmaxCont2)
        {
          dp3        = G4float(p3);
          rfac       = dp3/(G4float(nmaxCont2)+dp3);
          namean     = G4float(p3)*rfac;
          sa         = G4float(nmaxCont1)*rfac;
          na         = RandGauss::shoot(namean,sa);
          if (na > 0.)
          {
            alfa   = w1*G4float(nmaxCont2+p3)/(w1*G4float(nmaxCont2)+G4float(p3));
            alfa1  = alfa*log(alfa)/(alfa-1.);
            ea     = na*ipotFluct*alfa1;
            sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
            lossc += RandGauss::shoot(ea,sea);
          }
        }

        nb = G4int(G4float(p3)-na);
        if (nb > 0)
        {
          w2 = alfa*ipotFluct;
          w  = (Tm-w2)/Tm;      
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }
      }        
      loss += lossc;  
     }
    } 

  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4hLowEnergyIonisation::NuclearLossFluctuation(
                const G4DynamicParticle* particle,
                const G4Material* material,
                      G4double nuclearEnergyLoss) const
{
  G4double x = nuclearEnergyLoss ;
  return x ;
}
        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::PrintInfoDefinition() const
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n         Good description above the mean excitation energy.\n";
  comments += "         delta ray energy sampled from  differential Xsection.";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << LowestKineticEnergy / eV << " eV " 
         << " to " << HighestKineticEnergy / TeV << " TeV "
         << " in " << TotBin << " bins."
 << "\n        Parametrisation model for protons is  " 
 << theProtonTable
         << "\n        from " << protonLowEnergy / keV << " keV "
         << " to " << protonHighEnergy / MeV << " MeV " << "." << G4endl ;
  G4cout << "\n        Parametrisation model for antiprotons is  "
         << theAntiProtonTable
         << "\n        from " << antiProtonLowEnergy / keV << " keV "
         << " to " << antiProtonHighEnergy / MeV << " MeV " << "." << G4endl ;
  if(theBarkas){
  G4cout << "        Parametrization of the Barkas effect is switched on." 
         << G4endl ;
  }
  if(nStopping) {
  G4cout << "        Nuclear stopping power model is " << theNuclearTable 
         << G4endl ; 
  }

  G4bool printHead = true;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  
  // loop for materials 
  
  for (G4int j=0 ; j < numOfMaterials; j++) { 

    const G4Material* material= (*theMaterialTable)[j];
    G4double deltaCutNow = deltaCutInKineticEnergy[(material->GetIndex())] ;   
    G4double eexc = material->GetIonisation()->GetMeanExcitationEnergy();

    if(eexc > deltaCutNow) {
      if(printHead) {
        printHead = false ;

        G4cout << "           material       min.delta energy(keV) " << G4endl;
        G4cout << G4endl;
      }

      G4cout << G4std::setw(20) << material->GetName()
	     << G4std::setw(15) << eexc/keV << G4endl;
    }
  }    
}
