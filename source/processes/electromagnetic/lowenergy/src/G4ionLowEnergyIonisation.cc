// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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
//      ---------- G4ionLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 6 September 1999 
//                was made on the base of G4hLowEnergyIonisation class
// ************************************************************
// It is the extention of the ionisation process for the slow 
// charged ions.
// ************************************************************
//  6 September 1999 V.Ivanchenko create
// 30 September 1999 V.Ivanchenko minor upgrade
// ------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionLowEnergyIonisation.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionLowEnergyIonisation::G4ionLowEnergyIonisation(const G4String& processName)
  : G4hLowEnergyIonisation(processName),
  theIon (G4Proton::Proton())
{ 
  LowestKineticEnergy = 10.*eV ;
  HighestKineticEnergy = 100.*TeV ;
  TotBin = 200 ;
  MassRatio = 1.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionLowEnergyIonisation::~G4ionLowEnergyIonisation() 
{
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionLowEnergyIonisation::SetIonDefinition(G4ParticleDefinition* theIonType)
{
  theIon    = theIonType ;
  MassRatio = proton_mass_c2/(theIonType->GetPDGMass()) ;
  Charge    = (theIonType->GetPDGCharge())/eplus ;
  cout << "New ion with Q = " << Charge << "; MassR = " << MassRatio << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetLowEnergyForParametrisation(const G4Material* material) 

{
  // The low limit of paramerisation of ionisation energy from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Below this limit the free electron gas model is used

  // hadrons or ions with charge +-1
  if(Charge < 1.5) return ParamLowEnergy/MassRatio ;

  // helium or ions with charge = +2
  if(Charge < 2.5) return ParamLowEnergy ;

  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = material->GetNumberOfElements() ;
  G4double Z = 0.0, Norm = 0.0 ; 
  
  // only 1 element in the material
  if( 1 == NumberOfElements ) {
    Z = material->GetZ() ;

  //  loop for the elements in the material
  //  to find out average value of Z
  } else {
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        const G4Element* element = (*theElementVector)(iel) ;
        G4double Z2 = element->GetZ() ;
        const G4double W2 = theAtomicNumDensityVector[iel] ;
        Norm += W2 ;
        Z    += Z2 * W2 ;
      }
    Z  /= Norm ;
  }
  G4double E1 = 3.25 * keV ;
  G4double E2 = 25.0 * keV / pow(Z, 0.667) ;
  E1 = max (E1, E2) ;
  return max(ParamLowEnergy, E1) / MassRatio ; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetConstraints(const G4DynamicParticle *aParticle,
                                              G4Material *aMaterial)
{
  // returns the Step limit
  // dRoverRange is the max. allowed relative range loss in one step
  // it calculates dEdx and the range as well....

  G4double KineticEnergy,StepLimit;
  G4bool isOut ;

  Charge = aParticle->GetDefinition()->GetPDGCharge()/eplus ;

  KineticEnergy = aParticle->GetKineticEnergy();

  G4double massratio=proton_mass_c2/
           aParticle->GetDefinition()->GetPDGMass() ;

  G4double Tscaled= KineticEnergy*massratio ; 
  G4double ChargeSquare = GetIonEffChargeSquare(aMaterial,KineticEnergy,Charge) ;

     if(Charge>0.)
     {
       fRangeNow = G4EnergyLossTables::GetRange( theProton,
                                            Tscaled,aMaterial) ;
        fdEdx     = G4EnergyLossTables::GetDEDX( theProton,
                                            Tscaled,aMaterial) ;
     }
     else
     {
       fRangeNow = G4EnergyLossTables::GetRange( theAntiProton,
                                             Tscaled,aMaterial) ;
       fdEdx     = G4EnergyLossTables::GetDEDX( theAntiProton,
                                             Tscaled,aMaterial) ;
     }
     fdEdx     *= ChargeSquare ;
     fRangeNow /= (ChargeSquare*massratio) ;

  // compute the (random) Step limit ..............
  if(fRangeNow > finalRange)
  {
    StepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;

    //  randomise this value
    if(rndmStepFlag) StepLimit = 
                finalRange+(StepLimit-finalRange)*G4UniformRand() ;
    if(StepLimit > fRangeNow) StepLimit = fRangeNow ;
  }
  else StepLimit = fRangeNow ;


  return StepLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4ionLowEnergyIonisation::AlongStepDoIt( 
                              const G4Track& trackData,const G4Step& stepData) 
 // compute the energy loss after a step 
{
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4double finalT,Step,MeanLoss ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;
  
  // get the actual (true) Step length from stepData 
  Step = stepData.GetStepLength() ;

  aParticle = trackData.GetDynamicParticle() ;

  G4int index = aMaterial->GetIndex() ;
  G4double E = aParticle->GetKineticEnergy() ;
  G4double ParticleCharge = aParticle->GetDefinition()->GetPDGCharge() ;
  G4double ChargeSquare = GetIonEffChargeSquare(aMaterial, E, ParticleCharge) ;

  if(E < MinKineticEnergy) MeanLoss = E ;
  else
  {
    if(Step >= fRangeNow ) MeanLoss = E ;

    else if(( E > HighestKineticEnergy)||( E <= LowestKineticEnergy))
              MeanLoss = Step*fdEdx ; 
     
    else
    {
      if(Step>linLossLimit*fRangeNow)
      {
        G4double massratio=proton_mass_c2/
                 aParticle->GetDefinition()->GetPDGMass() ;

        G4double rscaled= fRangeNow*massratio*ChargeSquare ;
        G4double sscaled=   Step   *massratio*ChargeSquare ;

        if(Charge>0.)
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        else
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        MeanLoss /= (massratio*ChargeSquare) ;
      }
      else MeanLoss = Step*fdEdx ;
    }
  } 
  finalT = E - MeanLoss ;

  if(finalT < MinKineticEnergy) finalT = 0. ;

  //  now the loss with fluctuation
  if((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowestKineticEnergy))
  {
    MeanLoss /= ChargeSquare ;
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss)*ChargeSquare ;
    if (finalT < 0.) finalT = E-MeanLoss ;
  }

  //  kill the particle if the kinetic energy <= 0  
  if (finalT <= 0. )
  {
    finalT = 0.;
    if(aParticle->GetDefinition()->GetParticleName() == "proton")
      aParticleChange.SetStatusChange(fStopAndKill);
    else  
      aParticleChange.SetStatusChange(fStopButAlive); 
  } 

  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetIonParametrisedLoss(const G4Material* material, 
                                                          const G4double KinEnergy, 
                                                          const G4double DeltaRayCutNow)
{
  // Inicialisation
  G4double               Se = 0.0 ;
  G4double               Sn = 0.0 ;
  G4double          ionloss = 0.0 ; 
  G4double           ion125 = 0.0 ; 
  G4double  ExpStopPower125 = 0.0 ;
  G4double ReducedKinEnergy = KinEnergy * MassRatio ;
  G4double     ChargeSquare = GetIonEffChargeSquare(material, KinEnergy, Charge) ;
  G4double               Z1 = Charge ;
  G4double               A1 = ProtonMassAMU / MassRatio ;

  // First of all check tables for specific materials for ICRU_49 parametrisation
  G4int molecIndex = (MolecIsInICRU_R49p(material))+1; 
 
  if ((molecIndex > 0) && (DEDXtable == "ICRU_R49p")) {

    G4double NbOfAtomsPerVolume  = material->GetTotNbOfAtomsPerVolume();
    ionloss = GetStoppingPowerICRU_R49p(molecIndex, ReducedKinEnergy, "Mol")  
            * NbOfAtomsPerVolume * ZieglerFactor * ChargeSquare ;

    // Second - check the table for chemical factors
  } else {
    G4double ExpStopPower125 = MolecIsInZiegler1988(material); 
  }


  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = material->GetNumberOfElements() ;
  
  //  loop for the elements in the material
  //  calculation based on Bragg's rule 
  for (G4int iel=0; iel<NumberOfElements; iel++)
    {
      const G4Element* element = (*theElementVector)(iel) ;
      G4double Z2 = element->GetZ() ;
      G4double A2 = element->GetA()*mole/g ;
      G4int iz = int(Z2) ;
      if( iz <= 0 ) iz = 1 ;
      if( iz > 92 ) iz = 92 ; 
      
  // Electronic Stopping Power 
  // Choose the parametrisation using the table name
      
  // The "Ziegler1977H" table
      if(DEDXtable == "Ziegler1977H") { 
        Se = GetStoppingPower1977H(iz, ReducedKinEnergy)  
           * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	
        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125 += GetStoppingPower1977H(iz, 125.0*keV)  
                  * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}

	// Nuclear Stopping Power
        if(nStopping) {
          Sn += GetStoppingPower1977n(Z1, Z2, A1, A2, KinEnergy) 
              * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
  // The "Ziegler1977He" table
      } else if(DEDXtable == "Ziegler1977He") {
	G4double HeKinEnergy = ReducedKinEnergy*HeMassAMU/ProtonMassAMU ;
        Se = GetStoppingPower1977He(iz, HeKinEnergy) 
           * theAtomicNumDensityVector[iel]*ZieglerFactor 
           / GetHeEffChargeSquare(iz, HeKinEnergy) ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125 += GetStoppingPower1977H(iz, 125.0*keV) 
                  * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}

	// Nuclear Stopping Power
        if(nStopping) {
          Sn += GetStoppingPower1977n(Z1, Z2, A1, A2, KinEnergy) 
              * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
      
  // The "ICRU_R49p" table
      } else if(DEDXtable == "ICRU_R49p") {

        // The material is not in the list of materials
        if(molecIndex < 0) { 
          Se = GetStoppingPowerICRU_R49p(iz, ReducedKinEnergy, "Ele")  
             * theAtomicNumDensityVector[iel]*ZieglerFactor ;
        } 

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125 += GetStoppingPowerICRU_R49p(iz, 125.0*keV, "Ele") 
                  * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
	
	// Nuclear Stopping Power 
        if(nStopping) {
          Sn += GetStoppingPowerMoliere(Z1, Z2, A1, A2, KinEnergy) 
              * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}

  // The "ICRU_R49He" table
      } else if(DEDXtable == "ICRU_R49He") {
	G4double HeKinEnergy = ReducedKinEnergy*HeMassAMU/ProtonMassAMU ;
        Se = GetStoppingPowerICRU_R49He(iz, HeKinEnergy) 
           * theAtomicNumDensityVector[iel]*ZieglerFactor 
           / GetHeEffChargeSquare(iz, HeKinEnergy) ; 

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125 += GetStoppingPowerICRU_R49p(iz, 125.0*keV, "Ele")  
                  * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}

	// Nuclear Stopping Power
        if(nStopping) {
          Sn += GetStoppingPower1985n(Z1, Z2, A1, A2, KinEnergy) 
              * theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
		
  // The G4 beta version model
      } else if(DEDXtable == "UrbanModel") {
        Se = theAtomicNumDensityVector[iel]*GetUrbanModel(element, ReducedKinEnergy) ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125 += theAtomicNumDensityVector[iel]*GetUrbanModel(element, 125.0*keV) ;
	}
      }
      
      ionloss     += Se * ChargeSquare ;
    }

  // Chemical factor is taken into account
  if(ExpStopPower125 > 0.0) {
    ionloss *= GetChemicalFactor(ExpStopPower125, ReducedKinEnergy, ion125) ;
  }

  // Correction due to delta-electrons energy loss. 
  // Bethe-Bloch formulae was used. 
  if(DEDXtable != "UrbanModel") {
    ionloss -= GetDeltaRaysEnergy(material, ReducedKinEnergy, DeltaRayCutNow) 
             * ChargeSquare ;
  }

  // Nuclear Stopping Power
  if(nStopping) ionloss += Sn ;
    
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetIonBetheBlochLoss(const G4Material* material, 
                                                        const G4double KinEnergy,
                                                        const G4double DeltaRayCutNow)
{
  G4double ionloss ;
  G4double taul = material->GetIonisation()->GetTaul() ;
  G4double tau  = MassRatio*KinEnergy/proton_mass_c2 ;    // tau is relative energy
  G4double ChargeSquare = GetIonEffChargeSquare(material, KinEnergy, Charge) ;
  
  if ( tau < taul ) {
    
    //  low energy part , parametrised L.Urban energy loss formulae
    
    const G4ElementVector* theElementVector=
      material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
      material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
      material->GetNumberOfElements() ;
    
    ionloss = 0. ;
    
    //  loop for the elements in the material
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
	const G4Element* element = (*theElementVector)(iel) ;
	ionloss += GetUrbanModel(element, KinEnergy*MassRatio) * theAtomicNumDensityVector[iel] ;
      }
    
  } else {
    // Standard Bethe-Bloch formulae
    
    // some local variables 
    
    G4double gamma,bg2,beta2,Tmax,rcut,x,delta,sh ;
    G4double ElectronDensity = material->GetElectronDensity();
    G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
    G4double Eexc2 = Eexc*Eexc ;
    G4double Cden = material->GetIonisation()->GetCdensity();
    G4double Mden = material->GetIonisation()->GetMdensity();
    G4double Aden = material->GetIonisation()->GetAdensity();
    G4double X0den = material->GetIonisation()->GetX0density();
    G4double X1den = material->GetIonisation()->GetX1density();
    G4double* ShellCorrectionVector;
    ShellCorrectionVector = material->GetIonisation()->
      GetShellCorrectionVector();
    
    gamma = tau + 1.0 ;
    bg2 = tau*(tau+2.0) ;
    beta2 = bg2/(gamma*gamma) ;
    Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass) ;
    
    if ( DeltaRayCutNow < Tmax)
      rcut = DeltaRayCutNow/Tmax ;
    else
      rcut = 1.;
    
    ionloss = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)+log(rcut)-(1.+rcut)*beta2 ;
    
    // density correction 
    
    x = log(bg2)/twoln10 ;
    if ( x < X0den )
      delta = 0. ;
    else 
      {
	delta = twoln10*x - Cden ;
	if ( x < X1den )
	  delta += Aden*pow((X1den-x),Mden) ;
      } 
    
    // shell correction 
    
    if ( bg2 > bg2lim ) {
      sh = 0. ;      
      x = 1. ;
      for (G4int k=0; k<=2; k++) {
	x *= bg2 ;
	sh += ShellCorrectionVector[k]/x;
      }
    } else {
      sh = 0. ;      
      x = 1. ;
      for (G4int k=0; k<=2; k++) {
	x *= bg2lim ;
	sh += ShellCorrectionVector[k]/x;
      }
      sh *= log(tau/taul)/log(taulim/taul) ;     
    }
    
    // now you can compute the total ionisation loss
    
    ionloss -= delta + sh ;
    ionloss *= Factor*ElectronDensity*ChargeSquare/beta2 ;
  }
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetIonLossWithFluct(const G4DynamicParticle* aParticle,
                                                             G4Material* aMaterial,
                                                             G4double    MeanLoss)

//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is the same as in Glandz in Geant3.
{
  G4double ChargeSquare = Charge*Charge ;

  G4double loss = GetLossWithFluct(aParticle, aMaterial, MeanLoss/ChargeSquare) * ChargeSquare ;

  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionLowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n         Good description above the mean excitation energy.\n";
  comments += "         delta ray energy sampled from  differential Xsection.";
  
  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
							  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins."
         << "\n        Low energy losses approximation is taken from  " << DEDXtable
         << "\n        from " << G4BestUnit(ParamLowEnergy,"Energy")
         << " to " << G4BestUnit(ParamHighEnergy,"Energy") << "." << endl ;
  if(nStopping) {
    G4cout << "        Simulation of nuclear stopping is switched on.  \n" << endl ; 
  }
}











































