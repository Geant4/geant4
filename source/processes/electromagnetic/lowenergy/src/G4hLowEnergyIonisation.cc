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
// 23 May 2000    MG Pia  Clean up for QAO model 
// 28 July   1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli added ICRU parametrisations for protons  
// 20 August 1999 G.Mancinelli added ICRU tables for alpha 
// 31 August 1999 V.Ivanchenko update and cleen up 
// 30 Sept.  1999 V.Ivanchenko minor upgrade 
// 12 Dec.   1999 S. Chauvie added Barkas correction 
// 19 Jan.   2000 V.Ivanchenko minor changing in Barkas corrections
// 02 April  2000 S. Chauvie linearization of barkas effect
// 03 April  2000 V.Ivanchenko Nuclear Stopping power for antiprotons
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hLowEnergyIonisation.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::G4hLowEnergyIonisation(const G4String& processName)
  : G4hLowEnergyLoss(processName),
    theMeanFreePathTable(NULL),
    ParamLowEnergy(1.*keV),
    ParamHighEnergy(2.*MeV),
    DEDXtable("Ziegler1977H"),
    nStopping(false),
    pbarStop(false),
    theProton (G4Proton::Proton()),
    theAntiProton (G4AntiProton::AntiProton()),
    theElectron ( G4Electron::Electron() ),
    twoln10(2.*log(10.)),
    Factor(twopi_mc2_rcl2),
    bg2lim(0.0169), 
    taulim(8.4146e-3),
    RateMass(electron_mass_c2/proton_mass_c2),
    ProtonMassAMU(1.007276),
    HeMassAMU(4.0026),
    ZieglerFactor(eV*cm2*1.0e-15) 
{ 
    LowestKineticEnergy = 10.*eV ;
    HighestKineticEnergy = 100.*TeV ;
    TotBin = 200 ;
    MassRatio = 1.0 ;
    DeltaCutInKineticEnergy = 0; 

    qaoLoss = new G4QAOLowEnergyLoss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::~G4hLowEnergyIonisation() 
{
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }

  delete qaoLoss;
  qaoLoss = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetStoppingPowerTableName(const G4String& dedxTable)
{
  if(dedxTable == "Ziegler1977H") { 
    DEDXtable = "Ziegler1977H";
    ParamHighEnergy = 2.0*MeV;
    
  } else if(dedxTable == "Ziegler1977He") {
    DEDXtable = "Ziegler1977He";
    ParamHighEnergy = 2.0*MeV;
    
  } else if(dedxTable == "ICRU_R49p") {
    DEDXtable = "ICRU_R49p";
    ParamHighEnergy = 2.0*MeV;
    
    // set at 2 MeV. The ICRU report affirm their parametrisations are
    // valid up to 1 MeV for protons. They have used Ziegler-like
    // parametrisations up to an energy value depending on the material
    // (typical values from 0.1 to 0.8 MeV), used Bethe-Bloch for
    // energies over another value (from 0.5 to 3 MeV), and interpolated
    // in between. Unfortunately they give NO values for the paramteres
    // on the interpolation, so Ziegler-like parametrisation is here used
    // up to 2 MeV (better boundary conditions there wrt 1 MeV) and
    // Bethe-Bloch for higher values (applying continuity constraint) 
    
    ParamHighEnergy = 2.0*MeV;
    
  } else if(dedxTable == "ICRU_R49He") {
    DEDXtable = "ICRU_R49He";
    ParamHighEnergy = 2.0*MeV;

  } else if(dedxTable == "ICRU_R49PowersHe") {
    DEDXtable = "ICRU_R49PowersHe";
    ParamHighEnergy = 2.0*MeV;
    
  } else if(dedxTable == "UrbanModel") {
    DEDXtable = "UrbanModel";
    ParamHighEnergy = 2.0*MeV;
    
  } else {
  G4cout << "G4hLowEnergyIonisation Warning: There is no table with the name ="
         << dedxTable << G4endl; 
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetNuclearStoppingOn()
{
  nStopping = true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetNuclearStoppingOff()
{
  nStopping = false ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetAntiProtonStoppingOn()
{
  pbarStop = true ;
  ParamLowEnergy = 0.05*MeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetAntiProtonStoppingOff()
{
  pbarStop = false ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // Tables for different hadrons will be different because of
  // small difference in Tmax connected with RateMass
  // RateMass = electron_mass_c2 / (aParticleType.GetPDGMass()) ;

  // cuts for  electron ....................
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
  
  G4double LowEdgeEnergy , ionloss, ionlossBB;
  G4double paramA, paramB;
  static const G4MaterialTable* theMaterialTable=
    G4Material::GetMaterialTable();
  
  //  create table
  
  G4int numOfMaterials = theMaterialTable->length();
  
  if ( theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  //  loop for materials
  
  for (G4int J=0; J<numOfMaterials; J++)
    {
      
      // create physics vector and fill it
      
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                           HighestKineticEnergy, 
                                                           TotBin);
      
      // get material parameters needed for the energy loss calculation
      
      G4Material* material= (*theMaterialTable)[J];
        
      // get  electron cut in kin. energy for the material
      
      DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;
            
      // define constants A and B for this material  
      
      paramA  = GetParametrisedLoss(material, ParamLowEnergy, 
                                   DeltaCutInKineticEnergyNow)
				   /sqrt(ParamLowEnergy) ; 

      ionloss = GetParametrisedLoss(material, ParamHighEnergy, 
                                    DeltaCutInKineticEnergyNow) ;

      ionlossBB = GetBetheBlochLoss(material, ParamHighEnergy, 
                                    DeltaCutInKineticEnergyNow) ; 

      paramB =  ionloss/ionlossBB - 1.0 ; 
      
      // now comes the loop for the kinetic energy values
      
      for (G4int i = 0 ; i < TotBin ; i++)
	{
	  LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
	  
	  if ( LowEdgeEnergy < ParamHighEnergy ) {
	    //  low energy part , parametrised energy loss formulae
	    
	    if ( LowEdgeEnergy < ParamLowEnergy ) {
	      
	      // The model of free electron gas
	      ionloss = GetFreeElectronGasLoss(paramA, LowEdgeEnergy) ;
	      
	    } else {
	      // Parametrisation for intermediate energy range
	      ionloss = GetParametrisedLoss(material, LowEdgeEnergy, 
                                            DeltaCutInKineticEnergyNow) ;
	    }
	  } else {
	    
	    // high energy part , Bethe-Bloch formula
	    ionloss = GetBetheBlochLoss(material, LowEdgeEnergy, 
                                        DeltaCutInKineticEnergyNow) ; 
	    ionloss *= (1.0 + paramB*ParamHighEnergy/LowEdgeEnergy) ;
	  }

	  // now put the loss into the vector
	  aVector->PutValue(i,ionloss) ;
	}
      // insert vector for this material into the table
      theLossTable->insert(aVector) ;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetPreciseDEDX  (G4Material* aMaterial, 
                                 const G4double KinEnergy,
                                 const G4ParticleDefinition* aParticleType)
{
  // Calculation for different hadrons will be different because of
  // small difference in Tmax connected with RateMass
  //  RateMass = electron_mass_c2 / (aParticleType.GetPDGMass()) ;

  G4double ionloss, ionlossBB ;
  G4double paramA, paramB, dedx ;

  ParticleMass = aParticleType->GetPDGMass() ;
  Charge       = aParticleType->GetPDGCharge()/eplus ;
  MassRatio    = proton_mass_c2/ParticleMass ;

  G4double Tscaled = KinEnergy*MassRatio ; 
  G4double ChargeSquare = GetIonEffChargeSquare(aMaterial,KinEnergy,Charge) ;

  if(Tscaled > ParamHighEnergy) {
    if(Charge>0.) {

      dedx = G4EnergyLossTables::GetPreciseDEDX( theProton,Tscaled,aMaterial) 
           * ChargeSquare ;

    } else {

      dedx = G4EnergyLossTables::GetPreciseDEDX( theAntiProton,Tscaled,aMaterial) 
           * ChargeSquare ;
    }
  } else {

              
    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[(aMaterial->GetIndex())] ;   
    if ( Tscaled < ParamLowEnergy ) {
      
      // define constants A for this material  
      paramA  = GetParametrisedLoss(aMaterial, ParamLowEnergy, 
                                    DeltaCutInKineticEnergyNow)/sqrt(ParamLowEnergy) ; 
	      
      // The model of free electron gas
      ionloss = GetFreeElectronGasLoss(paramA, Tscaled) ;
	      
    } else {

      // Parametrisation for intermediate energy range
      ionloss = GetParametrisedLoss(aMaterial, Tscaled, 
                                    DeltaCutInKineticEnergyNow) ;
    }
    ionloss *= ChargeSquare ;

  }

  return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetPhysicsTableBining(G4double lowE, G4double highE,
						   G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE;
  TotBin = nBins ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)

  //  just call BuildLossTable+BuildLambdaTable
{
  ParticleMass = aParticleType.GetPDGMass() ;
  
  Charge = aParticleType.GetPDGCharge()/eplus ;
  
  G4double ElectronCutInRange = G4Electron::Electron()->GetCuts(); 
  
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

    if(Charge>0.) 
      {
	if( (ptableElectronCutInRange != ElectronCutInRange)  
	    || (theDEDXpTable == NULL))
	  {
	    BuildLossTable(aParticleType) ;
	    RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
	    CounterOfpProcess++;
	  }
      }
    else
      {
	if( (pbartableElectronCutInRange != ElectronCutInRange)  
	    || (theDEDXpbarTable == NULL))
	  {
	    BuildLossTable(aParticleType) ;
	    RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
	    CounterOfpbarProcess++;
	  }
      }

  BuildLambdaTable(aParticleType) ;
  
  BuildDEDXTable(aParticleType) ;
  
  if((&aParticleType == G4Proton::Proton()) )
        PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)

{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 
  
  G4double LowEdgeEnergy , Value ,sigma ;
  G4bool isOutRange ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  //Particle properties

  //ParticleMass = aParticleType.GetPDGMass() ;
  //G4double Charge = aParticle.GetPDGCharge()/eplus ;
  
  //create table
  
  G4int numOfMaterials = theMaterialTable->length();
  
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  
  theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);
  
  // get electron and particle cuts in kinetic energy
  
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
  
  // loop for materials 
  
  for (G4int J=0 ; J < numOfMaterials; J++)
    { 
      //create physics vector then fill it ....
      
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                           HighestKineticEnergy, 
                                                           TotBin);
      
      // compute the (macroscopic) cross section first
      
      const G4Material* material= (*theMaterialTable)[J];
      
      const G4ElementVector* theElementVector = material->GetElementVector() ;
      const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
      const G4int NumberOfElements = material->GetNumberOfElements() ;
      G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
      
      // get the electron kinetic energy cut for the actual material,
      //  it will be used in ComputeMicroscopicCrossSection
      // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
      //   ------------------------------------------------------
      
      DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;
      
      for ( G4int i = 0 ; i < TotBin ; i++ )
        {
	  LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
          G4double ChargeSquare = GetIonEffChargeSquare(material,LowEdgeEnergy,Charge) ;
	  
	  sigma = 0.0 ;
	  
	  for (G4int iel=0; iel<NumberOfElements; iel++ )
	    {
	      sigma +=  theAtomicNumDensityVector[iel]*ChargeSquare*
		ComputeMicroscopicCrossSection(aParticleType,
					       LowEdgeEnergy,
					       (*theElementVector)(iel)->GetZ(),
                                               Eexc ) ;
	    }
	  
	  // mean free path = 1./macroscopic cross section
	  
	  Value = sigma<=0 ? DBL_MAX : 1./sigma ;     
	  
	  aVector->PutValue(i, Value) ;
        }
      
      theMeanFreePathTable->insert(aVector);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber,
                                 G4double ExcEnergy)
{
  //******************************************************************
  // cross section formula is OK for spin=0, 1/2, 1 only !
  // *****************************************************************
  
  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = Z )
  
  G4double TotalEnergy, Beta2, Tmax, var ;
  G4double TotalCrossSection = 0.0 ;

  G4double Eexc2 = ExcEnergy*ExcEnergy ;
  
  // get particle data ...................................
  
  TotalEnergy=KineticEnergy + ParticleMass;
  
  // some kinematics......................
  
  Beta2  = KineticEnergy*(TotalEnergy+ParticleMass) / (TotalEnergy*TotalEnergy);
  var    = ParticleMass+electron_mass_c2;
  Tmax   = 2.*electron_mass_c2*KineticEnergy * (TotalEnergy+ParticleMass)
         / (var*var+2.*electron_mass_c2*KineticEnergy);

  // Validity range for delta electron cross section
  G4double DeltaCut = G4std::max(DeltaCutInKineticEnergyNow, ExcEnergy);
  
  // now you can calculate the total cross section ------------------
  
  if( Tmax > DeltaCut ) {

      var=DeltaCut/Tmax;
      TotalCrossSection = (1.-var*(1.-Beta2*log(var))) / DeltaCut ;
      G4double spin = aParticleType.GetPDGSpin() ;
      
      // +term for spin=1/2 particle
      if( 0.5 == spin )
	  TotalCrossSection +=  0.5 * (Tmax - DeltaCut) / (TotalEnergy*TotalEnergy);

      // +term for spin=1 particle
      else if( 0.9 < spin )
	  TotalCrossSection += -log(var)/(3.0*DeltaCut) +
                               (Tmax - DeltaCut) *
         ( (5.0+ 1.0/var)*0.25 / (TotalEnergy*TotalEnergy) - 
           Beta2 / (Tmax * DeltaCut) ) / 3.0 ;
	
      TotalCrossSection *= twopi_mc2_rcl2 * AtomicNumber / Beta2 ;
    }
  
  return TotalCrossSection ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetConstraints(const G4DynamicParticle *aParticle,
                                                  G4Material *aMaterial)
{
  // returns the Step limit
  // it calculates dEdx and the range as well 
  // based on Effective charge approach

  G4double KineticEnergy,StepLimit ;
  G4bool isOut ;

  theParticle = aParticle->GetDefinition() ;
  MassRatio = proton_mass_c2/(theParticle->GetPDGMass()) ;
  Charge    = (theParticle->GetPDGCharge())/eplus ;

  KineticEnergy = aParticle->GetKineticEnergy() ;

  // Scale the kinetic energy

  G4double Tscaled= KineticEnergy*MassRatio ; 
  G4double ChargeSquare = GetIonEffChargeSquare(aMaterial,KineticEnergy,Charge) ;
  G4double dx, s ;

  if(Charge>0.) {

    fdEdx = G4EnergyLossTables::GetDEDX( theProton, Tscaled, aMaterial) 
          * ChargeSquare ;

    fRangeNow = G4EnergyLossTables::GetRange( theProton, Tscaled, aMaterial) ;
    s = fRangeNow ;

    if(Tscaled < ParamHighEnergy) {
      // For Bragg's peak the limit in range is estimated 
      // in order to be inside linLossLimit on each step
      fdEdx = GetPreciseDEDX (aMaterial, KineticEnergy, theParticle) ;

      dx = G4EnergyLossTables::GetRange( theProton,
                                 ParamHighEnergy, aMaterial) * linLossLimit ;
      fRangeNow = G4std::min (fRangeNow, dx) ;
    }

    // Antiprotons and negative hadrons 
  } else {

    fdEdx = G4EnergyLossTables::GetDEDX( theAntiProton, Tscaled, aMaterial) 
          * ChargeSquare ;

    fRangeNow = G4EnergyLossTables::GetRange( theAntiProton, Tscaled, aMaterial) ;

    if(Tscaled < ParamHighEnergy) {
      // For Bragg's peak the limit in range is estimated 
      // in order to be inside linLossLimit on each step
      fdEdx = GetPreciseDEDX (aMaterial, KineticEnergy, theParticle) ;
      dx = G4EnergyLossTables::GetRange( theAntiProton,
                                 ParamHighEnergy, aMaterial) * linLossLimit ;
      fRangeNow = G4std::min (fRangeNow, dx) ;
    }
  }

  //
  fRangeNow /= (ChargeSquare*MassRatio) ;
  StepLimit  = fRangeNow ;

  // compute the (random) Step limit ..............
  if(fRangeNow > finalRange) {
    if(Tscaled > ParamHighEnergy ) {
      StepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;

      //  randomise this value
      if(rndmStepFlag) StepLimit = 
                finalRange+(StepLimit-finalRange)*G4UniformRand() ;
      if(StepLimit > fRangeNow) StepLimit = fRangeNow ;
    }
  }

  return StepLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::AlongStepDoIt( 
                              const G4Track& trackData, const G4Step& stepData) 

{
  // compute the energy loss after a step 

  G4double finalT = 0.0 ;

  aParticleChange.Initialize(trackData) ;
  G4Material* aMaterial = trackData.GetMaterial() ;
  
  // get the actual (true) Step length from stepData 
  const G4double Step = stepData.GetStepLength() ;

  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;

  G4int index = aMaterial->GetIndex() ;
  G4double E = aParticle->GetKineticEnergy() ;

  if( (aParticle->GetDefinition()) != theParticle ) {

    theParticle = aParticle->GetDefinition() ;
    MassRatio = proton_mass_c2/(theParticle->GetPDGMass()) ;
    Charge    = (theParticle->GetPDGCharge())/eplus ;
  }

  G4double Tscaled= E*MassRatio ; 
  G4double ChargeSquare = Charge*Charge ;
  G4double Eloss = 0.0 ;
  G4double Nloss = 0.0 ;

    if(E < MinKineticEnergy) Eloss = E ;
  
    else if(( E > HighestKineticEnergy)||( E <= LowestKineticEnergy))
              Eloss = Step*fdEdx ; 

    else if(Tscaled < ParamHighEnergy) {

      // Nuclear Stopping Power
      if(nStopping) {
        Nloss = GetNuclearDEDX(aMaterial, E, theParticle) ;
      }

      G4double E1 = E - Step*(fdEdx + Nloss) ;

      if(0.0 < E1) {
        Eloss = (fdEdx + GetPreciseDEDX (aMaterial, E1, theParticle))*Step*0.5 ;
        if(nStopping) {
          Nloss = (Nloss + GetNuclearDEDX (aMaterial, E1, theParticle))*Step*0.5 ;
	}
      } else Eloss = E ;

    } else if(Step >= fRangeNow ) Eloss = E ;
    
    else {

      if(Step>linLossLimit*fRangeNow) {

        G4double rscaled= fRangeNow*MassRatio*ChargeSquare ;
        G4double sscaled=   Step   *MassRatio*ChargeSquare ;

        if(Charge>0.)
        {
          Eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        else
        {
          Eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        Eloss /= (MassRatio*ChargeSquare) ;

      } else Eloss = Step*fdEdx ;
    }
 
  finalT = E - Eloss - Nloss;

  if(finalT > MinKineticEnergy) {

    //  now the electron loss with fluctuation
    if((EnlossFlucFlag) && (finalT < E) && (E > LowestKineticEnergy)) {
 
      Eloss = GetLossWithFluct(aParticle,aMaterial,Eloss/ChargeSquare)
               * ChargeSquare ;
      //        if(nStopping) {
      //      Nloss = GetNuclearLossWithFluct(theParticle,aMaterial,Nloss) ;
      // }
      finalT = E - Eloss - Nloss ;
    }    
  }

  //  kill the particle if the kinetic energy <= 0  
  if (finalT <= 0.0 )
  {
    finalT = 0.0 ;
    if(theParticle->GetParticleName() == "proton")
      aParticleChange.SetStatusChange(fStopAndKill);
    else  
      aParticleChange.SetStatusChange(fStopButAlive); 
  } 

  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::PostStepDoIt(const G4Track& trackData,   
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
  DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

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

G4double G4hLowEnergyIonisation::GetParametrisedLoss(G4Material* material, 
                                                     const G4double KinEnergy,
						     const G4double DeltaRayCutNow)

{

  G4double ionloss, ion, ionloss125, ion125;
  G4int molecIndex = 0 ;
  ionloss = 0.0 ;
   
  // start with aniproton management with quantum harmonic
  // oscillator loss model

  G4ParticleDefinition* pbarDef = G4AntiProton::AntiProtonDefinition();
  if ( pbarStop && (qaoLoss->IsInCharge(KinEnergy,pbarDef,material) ))
    {
      G4DynamicParticle dynamicPart();
      dynamicPart.SetDefinition(pbarDef);
      dynamicPart.SetKineticEnergy(KinEnergy);
      // does it require also the following?
      // dynamicPart.SetCharge(-1);

      ionloss = qaoloss->EnergyLoss(dynamicPart,material);
      ionloss -= GetDeltaRaysEnergy(material,KinEnergy,DeltaRayCutNow);
    }
  //  if(-0.5>Charge && pbarStop && qaoloss.IsMaterial(material->GetName())){
  //  ionloss = qaoloss.EnergyLoss(Charge, material, KinEnergy);
  //  ionloss-=GetDeltaRaysEnergy(material,KinEnergy,DeltaRayCutNow);
  return ionloss;
  }
  
  
    
  // First of all check tables for specific materials for ICRU_49 parametrisation
 
  // Ziegler parametrisation in ICRU49
  if ( DEDXtable == "ICRU_R49p" ) {

    molecIndex = (MolecIsInICRU_R49p(material))+1; 

    if ((molecIndex > 0)) {
      ionloss = GetMolecICRU_R49Loss(material, KinEnergy, 
                                     DeltaRayCutNow, molecIndex);

      // Barkas corrections for negative particles
      if( (-0.5 > Charge) && pbarStop) {
        ionloss += ComputeBarkasTerm( material, KinEnergy ) * (Charge - 1.0) ;
      }  
 
      if ( ionloss <= 0.0) ionloss = 0.0 ;
      return ionloss;
    }
  }

  // Powers parametrisation in ICRU49 
  if ( DEDXtable == "ICRU_R49PowersHe" ) {

    molecIndex = (MolecIsInICRU_R49PowersHe(material))+1;  
    if ( molecIndex > 0 ) {
      ionloss = GetMolecICRU_R49Loss(material, KinEnergy, 
                                     DeltaRayCutNow, molecIndex) ; 

      // Barkas corrections for negative particles
      if( (-0.5 > Charge) && pbarStop) {
        ionloss += ComputeBarkasTerm( material, KinEnergy ) * (Charge - 1.0) ;
      }  
 
      if ( ionloss <= 0.0) ionloss = 0.0 ;
      return ionloss;
    }
    DEDXtable = "ICRU_R49He" ;
  }

  G4double ExpStopPower125 = MolecIsInZiegler1988(material); 

  // Now cycle over elements - calculation based on Bragg's rule 

  // get elements in the actual material,
  const G4ElementVector* theElementVector  =material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector=material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements             = material->GetNumberOfElements() ;
  
  ionloss    = 0.0 ;
  ionloss125 = 0.0 ;
  
  //  loop for the elements in the material
  for (G4int iel=0; iel<NumberOfElements; iel++)
    {
      const G4Element* element = (*theElementVector)(iel) ;
      G4double Z2 = element->GetZ() ;
      G4int iz = int(Z2) ;
      if( iz <= 0 ) iz = 1 ;
      if( iz > 92 ) iz = 92 ; 
      
  // Electronic Stopping Power 
  // Choose the parametrisation using the table name
      
  // The "Ziegler1977H" table
      if(DEDXtable == "Ziegler1977H") { 
        ion = GetStoppingPower1977H(iz, KinEnergy) ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
	
        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125  = GetStoppingPower1977H(iz, 125.0*keV) ; 
          ion125 *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}

  // The "Ziegler1977He" table
      } else if(DEDXtable == "Ziegler1977He") {
	G4double HeKinEnergy = KinEnergy*HeMassAMU/ProtonMassAMU ;
        ion = GetStoppingPower1977He(iz, HeKinEnergy) /
              GetHeEffChargeSquare(iz, HeKinEnergy) ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125  = GetStoppingPower1977H(iz, 125.0*keV) ; 
          ion125 *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
      
  // The "ICRU_R49p" table
      } else if(DEDXtable == "ICRU_R49p") { 
        ion = GetStoppingPowerICRU_R49p(iz, KinEnergy, "Ele") ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125  = GetStoppingPowerICRU_R49p(iz, 125.0*keV, "Ele") ; 
          ion125 *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
	
  // The "ICRU_R49He" table
      } else if(DEDXtable == "ICRU_R49He") {
	G4double HeKinEnergy = KinEnergy*HeMassAMU/ProtonMassAMU ;
        ion = GetStoppingPowerICRU_R49He(iz, HeKinEnergy) /
              GetHeEffChargeSquare(iz, HeKinEnergy) ; 
        ion *= theAtomicNumDensityVector[iel]*ZieglerFactor ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125  = GetStoppingPowerICRU_R49p(iz, 125.0*keV, "Ele") ; 
          ion125 *= theAtomicNumDensityVector[iel]*ZieglerFactor ;
	}
		
  // The G4 beta version model
      } else if(DEDXtable == "UrbanModel") {
        ion = theAtomicNumDensityVector[iel]*GetUrbanModel(element, KinEnergy) ;

        // Chemical factor calculation
        if(ExpStopPower125 > 0.0){
          ion125  = theAtomicNumDensityVector[iel]*GetUrbanModel(element, 125.0*keV) ;
	}
      }
            
      ionloss    += ion ;
      ionloss125 += ion125 ;

    }

  // Chemical factor is taken into account
  if(ExpStopPower125 > 0.0){
    G4double x = GetChemicalFactor(ExpStopPower125, KinEnergy, ionloss125) ;
    ionloss *= x ;
  }
  
  // Correction due to delta-electrons energy loss. 
  // Bethe-Bloch formulae is used. 
  // Exception for the G4 beta version model, 
  // where the effect is already taken into account.
  if(DEDXtable != "UrbanModel") {
    ionloss    -= GetDeltaRaysEnergy(material, KinEnergy, DeltaRayCutNow) ;
  }
  
  // Correction term for the Barkas effect applied if pbarStop = true
  // and only for negative charged particles
  // Barkas term is taken into account in Ziegler/ICRU tables,
  // so for antiprotons a correction term must be multiplied by factor (-2) 
  if( (-0.5 > Charge) && pbarStop) {

    ionloss += ComputeBarkasTerm( material, KinEnergy ) * (Charge - 1.0) ;

  }  
 
  if ( ionloss <= 0.0) ionloss = 0.0 ;
  
  return ionloss;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetNuclearDEDX(G4Material* material, 
                                                const G4double KinEnergy,
                        		        const G4ParticleDefinition* aParticleType)

{

  G4double ionloss = 0.0 ;
 
  // Now cycle over elements - calculation based on Bragg's rule 

  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = material->GetNumberOfElements() ;

  MassRatio = proton_mass_c2/(aParticleType->GetPDGMass()) ;
  Charge    = (aParticleType->GetPDGCharge())/eplus ;

  G4double A1 = ProtonMassAMU/MassRatio ;
  
  //  loop for the elements in the material
  for (G4int iel=0; iel<NumberOfElements; iel++) {
      const G4Element* element = (*theElementVector)(iel) ;
      G4double Z1 = abs(Charge) ;
      G4double Z2 = element->GetZ() ;
      G4double A2 = element->GetA()*mole/g ;
      G4int iz = G4int(Z2) ;
      if( iz <= 0 ) iz = 1 ;
      if( iz > 92 ) iz = 92 ; 
  // Choose the parametrisation using the table name
      
  // The "Ziegler1977H" table
      if(DEDXtable == "Ziegler1977H") { 
	ionloss = GetStoppingPower1977n(Z1, Z2, A1, A2, KinEnergy) 
	        * theAtomicNumDensityVector[iel]*ZieglerFactor ;
    
  // The "ICRU_R49p" table
      //      } else if(DEDXtable == "ICRU_R49p") { 
      } else { 
	ionloss = GetStoppingPowerMoliere(Z1, Z2, A1, A2, KinEnergy)  
	        * theAtomicNumDensityVector[iel]*ZieglerFactor ;
      }
    }
  return ionloss;
}

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//Function to compute the Barkas term from:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397 
//
G4double G4hLowEnergyIonisation::ComputeBarkasTerm(const G4Material* material,
  				                   const G4double KinEnergy)

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

  // Internal variable for Kinetic Energy 
  // in order to keep Barkas correction to be constant below 500 keV 

  G4double KineticEnergy = KinEnergy;
  if( 500*keV > KineticEnergy ) KineticEnergy = 500*keV;

  // Information on particle and material
  
  G4double BarkasTerm = 0.0;
  G4double AMaterial = 0.0;
  G4double ZMaterial = 0.0;
  G4double RoMaterial = material->GetDensity()/6.2415063631e18;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int i=0;
  for (i = 0; i<material->GetNumberOfElements(); ++i) {

    AMaterial = (*theElementVector)(i)->GetA()*mole/g;
    ZMaterial = (*theElementVector)(i)->GetZ();
    
    G4double Beta = sqrt( 2.0*KineticEnergy / proton_mass_c2 );
    G4double X = ( (137.0*Beta) * (137.0*Beta) ) / ZMaterial;
  
    // Variables to compute L_1
    G4double Eta0Chi = 0.8;
    G4double EtaChi = Eta0Chi * ( 1.0 + 6.02*pow( ZMaterial,-1.19 ) );
    G4double W = ( EtaChi * pow( ZMaterial,1.0/6.0 ) ) / sqrt(X); 
    G4double FunctionOfW = 0.0;
    
    for(int IndexOfFTable=0; IndexOfFTable<47; IndexOfFTable++) {
    
      if(W<FTable[IndexOfFTable][0]) {
    
        if(0 == IndexOfFTable) {
     	  FunctionOfW = FTable[0][1] ;
        }
    
     	else if(46 == IndexOfFTable) {
          FunctionOfW = FTable[46][1] ;
        }
    
        else {
          FunctionOfW = ( W - FTable[IndexOfFTable-1][0]) 
	   		/ (FTable[IndexOfFTable][0] - FTable[IndexOfFTable-1][0]) ;    	
	  FunctionOfW *= (FTable[IndexOfFTable][1] - FTable[IndexOfFTable-1][1]) ;
	  FunctionOfW += FTable[IndexOfFTable-1][1] ;
	}
    
        break;
      }
    
    }
    
    G4double BarkasCoeffLbyARB = FunctionOfW / ( sqrt(ZMaterial) * pow(X,1.5) );
    BarkasTerm += BarkasCoeffLbyARB * ( 0.030708 * ZMaterial * RoMaterial )
               	          / ( AMaterial*Beta*Beta );
  }

  return BarkasTerm;
}
        
       
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetMolecICRU_R49Loss(const G4Material* material, 
					              const G4double KinEnergy, 
					              const G4double DeltaRayCutNow,
                                                      const G4int molecIndex)
{
    G4double NbOfAtomsPerVolume  = material->GetTotNbOfAtomsPerVolume();
    const G4int NumberOfElements = material->GetNumberOfElements() ;
  
    G4double ionloss = 0.0 ;
  
    if(DEDXtable == "ICRU_R49p") {
      ionloss = GetStoppingPowerICRU_R49p(molecIndex, KinEnergy, "Mol") ; 

    } else  if(DEDXtable == "ICRU_R49PowersHe") {
      G4double HeKinEnergy = KinEnergy*HeMassAMU/ProtonMassAMU ;
      ionloss = GetStoppingPowerICRU_R49PowersHe(molecIndex, HeKinEnergy) 
              / GetIonEffChargeSquare(material, HeKinEnergy, 2.0);
    }

    ionloss *= NbOfAtomsPerVolume*ZieglerFactor ;
    
    G4int nAtoms = 0;
    
    const G4int* theAtomsVector = material->GetAtomsVector() ;
    const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;

    //    const G4ElementVector* theElementVector=
    //      material->GetElementVector() ;
    
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
	nAtoms += theAtomsVector[iel];
      }
    ionloss /= nAtoms;
    
    // Correction due to delta-electrons energy loss Bethe-Bloch 
    // formula is used.
    ionloss -= GetDeltaRaysEnergy(material, KinEnergy, DeltaRayCutNow) ;
     
    if ( ionloss <= 0.) ionloss = 0. ;
    return ionloss;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetBetheBlochLoss(const G4Material* material, 
                                                   const G4double KinEnergy,
                                                   const G4double DeltaRayCutNow)
{
  G4double ionloss ;
  G4double taul = material->GetIonisation()->GetTaul() ;
  G4double tau  = KinEnergy/proton_mass_c2 ;         // tau is relative energy
  
  if ( tau < taul ) {
    
    //  low energy part , parametrised L.Urban energy loss formulae
    
    const G4ElementVector* theElementVector = material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements = material->GetNumberOfElements() ;
    
    ionloss = 0. ;
    
    //  loop for the elements in the material
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
	const G4Element* element = (*theElementVector)(iel) ;
	ionloss += GetUrbanModel(element, KinEnergy) * theAtomicNumDensityVector[iel] ;
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
    ShellCorrectionVector = material->GetIonisation()->GetShellCorrectionVector();
    
    gamma = tau + 1.0 ;
    bg2 = tau*(tau+2.0) ;
    beta2 = bg2/(gamma*gamma) ;
    Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass) ;
    
    // Validity range for delta electron cross section
    G4double DeltaCut = G4std::max(DeltaRayCutNow, Eexc);

    if ( DeltaCut < Tmax)
      rcut = DeltaCut/Tmax ;
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
    
    // now you can compute the total ionization loss
    
    ionloss -= delta + sh ;
    ionloss *= Factor*ElectronDensity/beta2 ;
  }

  // Barkas correction term is switch on 
  if( pbarStop) {
    ionloss += ComputeBarkasTerm( material, KinEnergy ) * Charge ;
  }  

  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetDeltaRaysEnergy(const G4Material* material, 
                                                    const G4double KinEnergy,
                                                    const G4double DeltaRayCutNow)
{
  G4double ionloss = 0.0 ;
  G4double tau  = KinEnergy/proton_mass_c2 ;         // tau is relative energy
  
  // some local variables 
  
  G4double gamma,bg2,beta2,Tmax,x ;
  G4double ElectronDensity = material->GetElectronDensity();
  
  gamma = tau + 1.0 ;
  bg2 = tau*(tau+2.0) ;
  beta2 = bg2/(gamma*gamma) ;
  Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass) ;

  G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double Eexc2 = Eexc*Eexc ;

  // Validity range for delta electron cross section
  G4double DeltaCut = G4std::max(DeltaRayCutNow, Eexc);

    
  if ( DeltaCut < Tmax) {
      x = DeltaCut / Tmax ;
      ionloss = ( beta2 * (x - 1.0) - log(x) ) * Factor * ElectronDensity / beta2 ;
  }
  return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetUrbanModel(const G4Element* element, 
                                                     G4double KinEnergy)
{
  // Parametrisation of low energy protons energy loss by L.Urban 
  // for the class G4hIonisation for the G4 beta version
  
  G4double tau  = KinEnergy/proton_mass_c2 ;         // tau is relative energy
  G4double ionloss = 0. ;
  
  if ( tau < element->GetIonisation()->GetTau0()) {  
    ionloss = ( element->GetIonisation()->GetAlow()*sqrt(tau)
		+ element->GetIonisation()->GetBlow()*tau) ;
  }
  else {
    ionloss = element->GetIonisation()->GetClow()/sqrt(tau) ;
  }
  
  return ionloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4hLowEnergyIonisation::MolecIsInICRU_R49p(const G4Material*
						 material)
{  
  const size_t NumberOfMolecula = 11 ;
  
  // Note the Graphite implementation is temporary. Right now it assume
  // that a chemical formula = Graphite is given. I'm not sure of the best
  // way to deal with this...
  
  static G4String Name[NumberOfMolecula] = {
    "AlO",                     "C_2O",                      "CH_4",  
    "(C_2H_4)_N-Polyethylene", "(C_2H_4)_N-Polypropylene",  "(C_8H_8)_N",  
    "C_3H_8",                  "SiO_2",                     "H_2O",
    "H_2O-Gas",                "Graphite"
  } ;
  
  G4String chFormula = material->GetChemicalFormula() ;
  G4int iMol;
  
  // Special treatment for water in gas state
  
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && "H_2O" == chFormula) {
    chFormula = "H_2O-Gas";
  }
  
  if (" " !=  chFormula ) {
    
    // Search for the material in the table
    for (iMol=0; iMol<NumberOfMolecula; iMol++) {
      if (chFormula == Name[iMol]) {
	return iMol;
      }
    }
  }
  iMol = -1;
  return iMol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4hLowEnergyIonisation::MolecIsInICRU_R49PowersHe(const G4Material*
						        material)
{  
  const size_t NumberOfMolecula = 30 ;
    
  static G4String Name[NumberOfMolecula] = {
   "H_2", "Be", "C", "Graphite", "N_2",
   "O_2", "Al", "Si", "Ar", "Cu",
   "Ge", "W", "Au", "Pb", "C_2H_2",
   "CO_2", "Cellulose-Nitrat", "C_2H_4", "LiF",
   "CH_4", "Nylon", "Polycarbonate", "(CH_2)_N-Polyetilene", "PMMA",
   "C_8H_8)_N", "SiO_2", "CsI", "H_2O", "H_2O-Gas"} ;      
  
  G4String chFormula = material->GetChemicalFormula() ;
  G4int iMol;
  
  // Special treatment for water in gas state
  
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && "H_2O" == chFormula) {
    chFormula = "H_2O-Gas";
  }
  
  if (" " !=  chFormula ) {
    
    // Search for the material in the table
    for (iMol=0; iMol<NumberOfMolecula; iMol++) {
      if (chFormula == Name[iMol]) {
	return iMol;
      }
    }
  }
  iMol = -1;
  return iMol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977H(G4int iz, G4double KinEnergy)
{
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom
  
  // The data and the fit from: 
  // H.H.Andersen & J.F.Ziegler Hydrogen Stopping Powers and
  // Ranges in All Elements, Vol.3, Pergamon Press, 1977
  
  G4double E = KinEnergy/(keV*ProtonMassAMU) ; // energy in keV/amu
  
  static G4double A[92][12] = {
    1.262,1.440,  242.6,12000.,0.115900,0.0005099,54360.0,-5.0520,2.0490,-0.30440,0.019660,-0.0004659,
    1.229,1.397,  484.5,5873.0,0.052250,0.0010200,24510.0,-2.1580,0.8278,-0.11720,0.007259,-0.000166,
    1.411,1.600,  725.6,3013.0,0.045780,0.0015300,21470.0,-0.5831,0.5620,-0.11830,0.009298,-0.0002498,
    2.248,2.590,  966.0, 153.8,0.034750,0.0020390,16300.0, 0.2779,0.1745,-0.05684,0.005155,-0.0001488,
    2.474,2.815, 1206.0,1060.0,0.028550,0.0025490,13450.0,-2.4450,1.2830,-0.22050,0.015600,-0.0003930,
    2.631,2.989, 1445.0, 957.2,0.028190,0.0030590,13220.0,-4.3800,2.0440,-0.32830,0.022210,-0.0005417,
    2.954,3.350, 1683.0,1900.0,0.025130,0.0035690,11790.0,-5.0540,2.3250,-0.37130,0.025060,-0.0006109,
    2.652,3.000, 1920.0,2000.0,0.022300,0.0040790,10460.0,-6.7340,3.0190,-0.47480,0.031710,-0.0007669,
    2.085,2.352, 2157.0,2634.0,0.018160,0.0045890,8517.0, -5.5710,2.4490,-0.37810,0.024830,-0.0005919,
    1.951,2.199, 2393.0,2699.0,0.015680,0.0050990,7353.0, -4.4080,1.8790,-0.28140,0.017960,-0.0004168,
    2.542,2.869, 2628.0,1854.0,0.014720,0.0056090,6905.0, -4.9590,2.0730,-0.30540,0.019210,-0.0004403,
    3.792,4.293, 2862.0,1009.0,0.013970,0.0061180,6551.0, -5.5100,2.2660,-0.32950,0.020470,-0.0004637,
    4.154,4.739, 2766.0, 164.5,0.020230,0.0066280,6309.0, -6.0610,2.4600,-0.35350,0.021730,-0.0004871,
    4.150,4.700, 3329.0, 550.0,0.013210,0.0071380,6194.0, -6.2940,2.5380,-0.36280,0.022200,-0.0004956,
    3.232,3.647, 3561.0,1560.0,0.012670,0.0076480,5942.0, -6.5270,2.6160,-0.37210,0.022670,-0.0005040,
    3.447,3.891, 3792.0,1219.0,0.012110,0.0081580,5678.0, -6.7610,2.6940,-0.38140,0.023140,-0.0005125,
    5.047,5.714, 4023.0, 878.6,0.011780,0.0086680,5524.0, -6.9940,2.7730,-0.39070,0.023610,-0.0005209,
    5.731,6.500, 4253.0, 530.0,0.011230,0.0091780,5268.0, -7.2270,2.8510,-0.40000,0.024070,-0.0005294,
    5.151,5.833, 4482.0, 545.7,0.011290,0.0096870,5295.0, -7.4400,2.9230,-0.40940,0.024620,-0.0005411,
    5.521,6.252, 4710.0, 553.3,0.011120,0.0102000,5214.0, -7.6530,2.9950,-0.41870,0.025160,-0.0005529,
    5.201,5.884, 4938.0, 560.9,0.009995,0.0107100,4688.0, -8.0120,3.1230,-0.43500,0.026050,-0.0005707,
    4.862,5.496, 5165.0, 568.5,0.009474,0.0112200,4443.0, -8.3710,3.2510,-0.45130,0.026940,-0.0005886,
    4.480,5.055, 5391.0, 952.3,0.009117,0.0117300,4276.0, -8.7310,3.3790,-0.46760,0.027830,-0.0006064,
    3.983,4.489, 5616.0,1336.0,0.008413,0.0122400,3946.0, -9.0900,3.5070,-0.48380,0.028720,-0.0006243,
    3.469,3.907, 5725.0,1461.0,0.008829,0.0127500,3785.0, -9.4490,3.6350,-0.50010,0.029610,-0.0006421,
    3.519,3.963, 6065.0,1243.0,0.007782,0.0132600,3650.0, -9.8090,3.7630,-0.51640,0.030500,-0.0006600,
    3.140,3.535, 6288.0,1372.0,0.007361,0.0137700,3453.0,-10.1700,3.8910,-0.53270,0.031390,-0.0006779,
    3.553,4.004, 6205.0, 555.1,0.008763,0.0142800,3297.0,-10.5300,4.0190,-0.54900,0.032290,-0.0006957,
    3.696,4.175, 4673.0, 387.8,0.021880,0.0147900,3174.0,-11.1800,4.2520,-0.57910,0.033990,-0.0007314,
    4.210,4.750, 6953.0, 295.2,0.006809,0.0153000,3194.0,-11.5700,4.3940,-0.59800,0.035060,-0.0007537,
    5.041,5.697, 7173.0, 202.6,0.006725,0.0158100,3154.0,-11.9500,4.5370,-0.61690,0.036130,-0.0007759,
    5.554,6.300, 6496.0, 110.0,0.009689,0.0163200,3097.0,-12.3400,4.6800,-0.63580,0.037210,-0.0007981,
    5.323,6.012, 7611.0, 292.5,0.006447,0.0168300,3024.0,-12.7200,4.8230,-0.65470,0.038280,-0.0008203,
    5.847,6.656, 7395.0, 117.5,0.007684,0.0173400,3006.0,-13.1100,4.9650,-0.67350,0.039350,-0.0008425,
    5.611,6.335, 8046.0, 365.2,0.006244,0.0178500,2928.0,-13.4000,5.0830,-0.69060,0.040420,-0.0008675,
    6.411,7.250, 8262.0, 220.0,0.006087,0.0183600,2855.0,-13.6900,5.2000,-0.70760,0.041500,-0.0008925,
    5.694,6.429, 8478.0, 292.9,0.006087,0.0188600,2855.0,-13.9200,5.2660,-0.71400,0.041730,-0.0008943,
    6.339,7.159, 8693.0, 330.3,0.006003,0.0193700,2815.0,-14.1400,5.3310,-0.72050,0.041960,-0.0008962,
    6.407,7.234, 8907.0, 367.8,0.005889,0.0198800,2762.0,-14.3600,5.3970,-0.72690,0.042190,-0.0008980,
    6.734,7.603, 9120.0, 405.2,0.005765,0.0203900,2704.0,-14.5900,5.4630,-0.73330,0.042420,-0.0008998,
    6.902,7.791, 9333.0, 442.7,0.005587,0.0209000,2621.0,-16.2200,6.0940,-0.82250,0.047910,-0.0010240,
    6.425,7.248, 9545.0, 480.2,0.005367,0.0214100,2517.0,-17.8500,6.7250,-0.91160,0.053390,-0.0011480,
    6.799,7.671, 9756.0, 517.6,0.005315,0.0219200,2493.0,-17.9600,6.7520,-0.91350,0.053410,-0.001147,
    6.108,6.887, 9966.0, 555.1,0.005151,0.0224300,2416.0,-18.0700,6.7790,-0.91540,0.053420,-0.0011450,
    5.924,6.677,10180.0, 592.5,0.004919,0.0229400,2307.0,-18.1800,6.8060,-0.91730,0.053430,-0.0011430,
    5.238,5.900,10380.0, 630.0,0.004758,0.0234500,2231.0,-18.2800,6.8330,-0.91920,0.053450,-0.0011420,
    5.623,6.354, 7160.0, 337.6,0.013940,0.0239600,2193.0,-18.3900,6.8600,-0.92110,0.053460,-0.0011400,
    5.814,6.554,10800.0, 355.5,0.004626,0.0244700,2170.0,-18.6200,6.9150,-0.92430,0.053400,-0.0011340,
    6.230,7.024,11010.0, 370.9,0.004540,0.0249800,2129.0,-18.8500,6.9690,-0.92750,0.053350,-0.0011270,
    6.410,7.227,11210.0, 386.4,0.004474,0.0254900,2099.0,-19.0700,7.0240,-0.93080,0.053290,-0.0011210,
    7.500,8.480, 8608.0, 348.0,0.009074,0.0260000,2069.0,-19.5700,7.2250,-0.96030,0.055180,-0.0011650,
    6.979,7.871,11620.0, 392.4,0.004402,0.0265100,2065.0,-20.0700,7.4260,-0.98990,0.057070,-0.0012090,
    7.725,8.716,11830.0, 394.8,0.004376,0.0270200,2052.0,-20.5600,7.6270,-1.01900,0.058960,-0.0012540,
    8.231,9.289,12030.0, 397.3,0.004384,0.0275300,2056.0,-21.0600,7.8280,-1.04900,0.060850,-0.0012980,
    7.287,8.218,12230.0, 399.7,0.004447,0.0280400,2086.0,-20.4000,7.5400,-1.00400,0.057820,-0.0012240,
    7.899,8.911,12430.0, 402.1,0.004511,0.0285500,2116.0,-19.7400,7.2520,-0.95880,0.054790,-0.0011510,
    8.041,9.071,12630.0, 404.5,0.004540,0.0290600,2129.0,-19.0800,6.9640,-0.91360,0.051760,-0.0010770,
    7.489,8.444,12830.0, 406.9,0.004420,0.0295700,2073.0,-18.4300,6.6770,-0.86840,0.048720,-0.0010030,
    7.291,8.219,13030.0, 409.3,0.004298,0.0300800,2016.0,-17.7700,6.3890,-0.82330,0.045690,-0.0009292,
    7.098,8.000,13230.0, 411.8,0.004182,0.0305900,1962.0,-17.1100,6.1010,-0.77810,0.042660,-0.0008553,
    6.910,7.786,13430.0, 414.2,0.004058,0.0311000,1903.0,-16.4500,5.8130,-0.73300,0.039630,-0.0007815,
    6.728,7.580,13620.0, 416.6,0.003976,0.0316100,1865.0,-15.7900,5.5260,-0.68780,0.036600,-0.0007077,
    6.551,7.380,13820.0, 419.0,0.003877,0.0321200,1819.0,-15.1300,5.2380,-0.64260,0.033570,-0.0006339,
    6.739,7.592,14020.0, 421.4,0.003863,0.0326300,1812.0,-14.4700,4.9500,-0.59750,0.030530,-0.0005601,
    6.212,6.996,14210.0, 423.9,0.003725,0.0331400,1747.0,-14.5600,4.9840,-0.60220,0.030820,-0.0005668,
    5.517,6.210,14400.0, 426.3,0.003632,0.0336500,1703.0,-14.6500,5.0180,-0.60690,0.031110,-0.0005734,
    5.219,5.874,14600.0, 428.7,0.003498,0.0341600,1640.0,-14.7400,5.0510,-0.61170,0.031410,-0.0005801,
    5.071,5.706,14790.0, 433.0,0.003405,0.0346700,1597.0,-14.8300,5.0850,-0.61640,0.031700,-0.0005867,
    4.926,5.542,14980.0, 433.5,0.003342,0.0351800,1567.0,-14.9100,5.1190,-0.62110,0.031990,-0.0005933,
    4.787,5.386,15170.0, 435.9,0.003292,0.0356900,1544.0,-15.0000,5.1530,-0.62580,0.032280,-0.0006000,
    4.893,5.505,15360.0, 438.4,0.003243,0.0362000,1521.0,-15.0900,5.1860,-0.63050,0.032570,-0.0006066,
    5.028,5.657,15550.0, 440.8,0.003195,0.0367100,1499.0,-15.1800,5.2200,-0.63530,0.032860,-0.0006133,
    4.738,5.329,15740.0, 443.2,0.003186,0.0372200,1494.0,-15.2700,5.2540,-0.64000,0.033150,-0.0006199,
    4.574,5.144,15930.0, 442.4,0.003144,0.0377300,1475.0,-15.6700,5.3920,-0.65770,0.034180,-0.0006426,
    5.200,5.851,16120.0, 441.6,0.003122,0.0382400,1464.0,-16.0700,5.5290,-0.67550,0.035210,-0.0006654,
    5.070,5.704,16300.0, 440.9,0.003082,0.0387500,1446.0,-16.4700,5.6670,-0.69320,0.036240,-0.0006881,
    4.945,5.563,16490.0, 440.1,0.002965,0.0392600,1390.0,-16.8800,5.8040,-0.71100,0.037270,-0.0007109,
    4.476,5.034,16670.0, 439.3,0.002871,0.0397700,1347.0,-17.2800,5.9420,-0.72870,0.038300,-0.0007336,
    4.856,5.460,18320.0, 438.5,0.002542,0.0402800,1354.0,-17.0200,5.8460,-0.71490,0.037400,-0.0007114,
    4.308,4.843,17040.0, 487.8,0.002882,0.0407900,1352.0,-17.8400,6.1830,-0.76590,0.040760,-0.0007925,
    4.723,5.311,17220.0, 537.0,0.002913,0.0413000,1366.0,-18.6600,6.5200,-0.81690,0.044110,-0.0008737,
    5.319,5.982,17400.0, 586.3,0.002871,0.0418100,1347.0,-19.4800,6.8570,-0.86780,0.047470,-0.0009548,
    5.956,6.700,17800.0, 677.0,0.002660,0.0423200,1336.0,-19.5500,6.8710,-0.86860,0.047480,-0.0009544,
    6.158,6.928,17770.0, 586.3,0.002812,0.0428300,1319.0,-19.6200,6.8840,-0.86940,0.047480,-0.0009540,
    6.204,6.979,17950.0, 586.3,0.002776,0.0433400,1302.0,-19.6900,6.8980,-0.87020,0.047490,-0.0009536,
    6.181,6.954,18120.0, 586.3,0.002748,0.0438500,1289.0,-19.7600,6.9120,-0.87100,0.047490,-0.0009532,
    6.949,7.820,18300.0, 586.3,0.002737,0.0443600,1284.0,-19.8300,6.9260,-0.87180,0.047500,-0.0009528,
    7.506,8.448,18480.0, 586.3,0.002727,0.0448700,1279.0,-19.9000,6.9400,-0.87260,0.047510,-0.0009524,
    7.649,8.609,18660.0, 586.3,0.002697,0.0453800,1265.0,-19.9700,6.9530,-0.87330,0.047510,-0.0009520,
    7.710,8.679,18830.0, 586.3,0.002641,0.0458900,1239.0,-20.0400,6.9670,-0.87410,0.047520,-0.0009516,
    7.407,8.336,19010.0, 586.3,0.002603,0.0464000,1221.0,-20.1100,6.9810,-0.87490,0.047520,-0.0009512,
    7.290,8.204,19180.0, 586.3,0.002573,0.0469100,1207.0,-20.1800,6.9950,-0.87570,0.047530,-0.0009508
  };
  
  if ( E < 10.0 ) {
    ionloss = A[i][0] * sqrt(E) ;
    
  } else if ( E < 1000.0 ) {
    G4double Slow  = A[i][1] * pow(E, 0.45) ;
    G4double Shigh = log( 1.0 + A[i][3]/E + A[i][4]*E ) * A[i][2]/E ;
    ionloss = Slow*Shigh / (Slow + Shigh) ; 
    
  } else {
    G4double le = log(E) ;
    G4double gam = 1.0 + KinEnergy / proton_mass_c2 ;
    G4double beta2 = 1.0 - 1.0/ (gam*gam) ;
    ionloss = ( log(A[i][6]*beta2/(1.0 - beta2)) - beta2 -
		A[i][7] - A[i][8]*le - A[i][9]*le*le - A[i][10]*le*le*le -
		A[i][11]*le*le*le*le ) * A[i][5]/beta2 ;
  }

  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPowerICRU_R49p(G4int iz,
							   G4double
							   KinEnergy,
							   G4String compound)
{
  G4double ionloss ;
  
  G4int i = iz-1 ;  // index of atom
  
  // The data and the fit from: 
  // ICRU Report 49, 1993
  
  G4double KinE = KinEnergy/(keV*ProtonMassAMU) ; // energy in keV/amu 
  
  static G4double A[92][5] = {
    1.254E+0, 1.440E+0, 2.426E+2, 1.200E+4, 1.159E-1,
    1.229E+0, 1.397E+0, 4.845E+2, 5.873E+3, 5.225E-2,
    1.411E+0, 1.600E+0, 7.256E+2, 3.013E+3, 4.578E-2,
    2.248E+0, 2.590E+0, 9.660E+2, 1.538E+2, 3.475E-2,
    2.474E+0, 2.815E+0, 1.206E+3, 1.060E+3, 2.855E-2,
    2.631E+0, 2.601E+0, 1.701E+3, 1.279E+3, 1.638E-2,
    2.954E+0, 3.350E+0, 1.683E+3, 1.900E+3, 2.513E-2,
    2.652E+0, 3.000E+0, 1.920E+3, 2.000E+3, 2.230E-2,
    2.085E+0, 2.352E+0, 2.157E+3, 2.634E+3, 1.816E-2,
    1.951E+0, 2.199E+0, 2.393E+3, 2.699E+3, 1.568E-2,
    2.542E+0, 2.869E+0, 2.628E+3, 1.854E+3, 1.472E-2,
    3.791E+0, 4.293E+0, 2.862E+3, 1.009E+3, 1.397E-2,
    4.154E+0, 4.739E+0, 2.766E+3, 1.645E+2, 2.023E-2,
    4.914E+0, 5.598E+0, 3.193E+3, 2.327E+2, 1.419E-2,
    3.232E+0, 3.647E+0, 3.561E+3, 1.560E+3, 1.267E-2,
    3.447E+0, 3.891E+0, 3.792E+3, 1.219E+3, 1.211E-2,
    5.301E+0, 6.008E+0, 3.969E+3, 6.451E+2, 1.183E-2,
    5.731E+0, 6.500E+0, 4.253E+3, 5.300E+2, 1.123E-2,
    5.152E+0, 5.833E+0, 4.482E+3, 5.457E+2, 1.129E-2,
    5.521E+0, 6.252E+0, 4.710E+3, 5.533E+2, 1.112E-2,
    5.201E+0, 5.884E+0, 4.938E+3, 5.609E+2, 9.995E-3,
    4.858E+0, 5.489E+0, 5.260E+3, 6.511E+2, 8.930E-3,
    4.479E+0, 5.055E+0, 5.391E+3, 9.523E+2, 9.117E-3,
    3.983E+0, 4.489E+0, 5.616E+3, 1.336E+3, 8.413E-3,
    3.469E+0, 3.907E+0, 5.725E+3, 1.461E+3, 8.829E-3,
    3.519E+0, 3.963E+0, 6.065E+3, 1.243E+3, 7.782E-3,
    3.140E+0, 3.535E+0, 6.288E+3, 1.372E+3, 7.361E-3,
    3.553E+0, 4.004E+0, 6.205E+3, 5.551E+2, 8.763E-3,
    3.696E+0, 4.194E+0, 4.649E+3, 8.113E+1, 2.242E-2,
    4.210E+0, 4.750E+0, 6.953E+3, 2.952E+2, 6.809E-3,
    5.041E+0, 5.697E+0, 7.173E+3, 2.026E+2, 6.725E-3,
    5.554E+0, 6.300E+0, 6.496E+3, 1.100E+2, 9.689E-3,
    5.323E+0, 6.012E+0, 7.611E+3, 2.925E+2, 6.447E-3,
    5.874E+0, 6.656E+0, 7.395E+3, 1.175E+2, 7.684E-3,
    6.658E+0, 7.536E+0, 7.694E+3, 2.223E+2, 6.509E-3,
    6.413E+0, 7.240E+0, 1.185E+4, 1.537E+2, 2.880E-3,
    5.694E+0, 6.429E+0, 8.478E+3, 2.929E+2, 6.087E-3,
    6.339E+0, 7.159E+0, 8.693E+3, 3.303E+2, 6.003E-3,
    6.407E+0, 7.234E+0, 8.907E+3, 3.678E+2, 5.889E-3,
    6.734E+0, 7.603E+0, 9.120E+3, 4.052E+2, 5.765E-3,
    6.901E+0, 7.791E+0, 9.333E+3, 4.427E+2, 5.587E-3,
    6.424E+0, 7.248E+0, 9.545E+3, 4.802E+2, 5.376E-3,
    6.799E+0, 7.671E+0, 9.756E+3, 5.176E+2, 5.315E-3,
    6.109E+0, 6.887E+0, 9.966E+3, 5.551E+2, 5.151E-3,
    5.924E+0, 6.677E+0, 1.018E+4, 5.925E+2, 4.919E-3,
    5.238E+0, 5.900E+0, 1.038E+4, 6.300E+2, 4.758E-3,
    5.345E+0, 6.038E+0, 6.790E+3, 3.978E+2, 1.676E-2,
    5.814E+0, 6.554E+0, 1.080E+4, 3.555E+2, 4.626E-3,
    6.229E+0, 7.024E+0, 1.101E+4, 3.709E+2, 4.540E-3,
    6.409E+0, 7.227E+0, 1.121E+4, 3.864E+2, 4.474E-3,
    7.500E+0, 8.480E+0, 8.608E+3, 3.480E+2, 9.074E-3,
    6.979E+0, 7.871E+0, 1.162E+4, 3.924E+2, 4.402E-3,
    7.725E+0, 8.716E+0, 1.183E+4, 3.948E+2, 4.376E-3,
    8.337E+0, 9.425E+0, 1.051E+4, 2.696E+2, 6.206E-3,
    7.287E+0, 8.218E+0, 1.223E+4, 3.997E+2, 4.447E-3,
    7.899E+0, 8.911E+0, 1.243E+4, 4.021E+2, 4.511E-3,
    8.041E+0, 9.071E+0, 1.263E+4, 4.045E+2, 4.540E-3,
    7.488E+0, 8.444E+0, 1.283E+4, 4.069E+2, 4.420E-3,
    7.291E+0, 8.219E+0, 1.303E+4, 4.093E+2, 4.298E-3,
    7.098E+0, 8.000E+0, 1.323E+4, 4.118E+2, 4.182E-3,
    6.909E+0, 7.786E+0, 1.343E+4, 4.142E+2, 4.058E-3,
    6.728E+0, 7.580E+0, 1.362E+4, 4.166E+2, 3.976E-3,
    6.551E+0, 7.380E+0, 1.382E+4, 4.190E+2, 3.877E-3,
    6.739E+0, 7.592E+0, 1.402E+4, 4.214E+2, 3.863E-3,
    6.212E+0, 6.996E+0, 1.421E+4, 4.239E+2, 3.725E-3,
    5.517E+0, 6.210E+0, 1.440E+4, 4.263E+2, 3.632E-3,
    5.220E+0, 5.874E+0, 1.460E+4, 4.287E+2, 3.498E-3,
    5.071E+0, 5.706E+0, 1.479E+4, 4.330E+2, 3.405E-3,
    4.926E+0, 5.542E+0, 1.498E+4, 4.335E+2, 3.342E-3,
    4.788E+0, 5.386E+0, 1.517E+4, 4.359E+2, 3.292E-3,
    4.893E+0, 5.505E+0, 1.536E+4, 4.384E+2, 3.243E-3,
    5.028E+0, 5.657E+0, 1.555E+4, 4.408E+2, 3.195E-3,
    4.738E+0, 5.329E+0, 1.574E+4, 4.432E+2, 3.186E-3,
    4.587E+0, 5.160E+0, 1.541E+4, 4.153E+2, 3.406E-3,
    5.201E+0, 5.851E+0, 1.612E+4, 4.416E+2, 3.122E-3,
    5.071E+0, 5.704E+0, 1.630E+4, 4.409E+2, 3.082E-3,
    4.946E+0, 5.563E+0, 1.649E+4, 4.401E+2, 2.965E-3,
    4.477E+0, 5.034E+0, 1.667E+4, 4.393E+2, 2.871E-3,
    4.844E+0, 5.458E+0, 7.852E+3, 9.758E+2, 2.077E-2,
    4.307E+0, 4.843E+0, 1.704E+4, 4.878E+2, 2.882E-3,
    4.723E+0, 5.311E+0, 1.722E+4, 5.370E+2, 2.913E-3,
    5.319E+0, 5.982E+0, 1.740E+4, 5.863E+2, 2.871E-3,
    5.956E+0, 6.700E+0, 1.780E+4, 6.770E+2, 2.660E-3,
    6.158E+0, 6.928E+0, 1.777E+4, 5.863E+2, 2.812E-3,
    6.203E+0, 6.979E+0, 1.795E+4, 5.863E+2, 2.776E-3,
    6.181E+0, 6.954E+0, 1.812E+4, 5.863E+2, 2.748E-3,
    6.949E+0, 7.820E+0, 1.830E+4, 5.863E+2, 2.737E-3,
    7.506E+0, 8.448E+0, 1.848E+4, 5.863E+2, 2.727E-3,
    7.648E+0, 8.609E+0, 1.866E+4, 5.863E+2, 2.697E-3,
    7.711E+0, 8.679E+0, 1.883E+4, 5.863E+2, 2.641E-3,
    7.407E+0, 8.336E+0, 1.901E+4, 5.863E+2, 2.603E-3,
    7.290E+0, 8.204E+0, 1.918E+4, 5.863E+2, 2.673E-3
  };
  
  static G4double Amol[11][5] = {
    1.187E+1, 1.343E+1, 1.069E+4, 7.723E+2, 2.153E-2, 
    7.802E+0, 8.814E+0, 8.303E+3, 7.446E+2, 7.966E-3, 
    7.294E+0, 8.284E+0, 5.010E+3, 4.544E+2, 8.153E-3, 
    8.646E+0, 9.800E+0, 7.066E+3, 4.581E+2, 9.383E-3, 
    1.286E+1, 1.462E+1, 5.625E+3, 2.621E+3, 3.512E-2, 
    3.229E+1, 3.696E+1, 8.918E+3, 3.244E+3, 1.273E-1, 
    1.604E+1, 1.825E+1, 6.967E+3, 2.307E+3, 3.775E-2, 
    8.049E+0, 9.099E+0, 9.257E+3, 3.846E+2, 1.007E-2, 
    4.015E+0, 4.542E+0, 3.955E+3, 4.847E+2, 7.904E-3, 
    4.571E+0, 5.173E+0, 4.346E+3, 4.779E+2, 8.572E-3, 
    2.631E+0, 2.601E+0, 1.701E+3, 1.279E+3, 1.638E-2 
  };
  
  static G4double Coeff[5] = {
    0.,0.,0.,0.,0.
  };
  
  if (compound == "Ele") {
    for (G4int iCoef=0; iCoef<5; iCoef++) {
      Coeff[iCoef] = A[i][iCoef];
    }
  }
  else if (compound == "Mol") {
    for (G4int iCoef=0; iCoef<5; iCoef++) {
      Coeff[iCoef] = Amol[i][iCoef];
    }
  }
  
  // carbon is a particular case 
  // It should be done correctly, but ICRU doesn't give details of their
  // complete fit
  
  G4String carbonType = " ";  
   
  if ((compound == "Ele") && (i == 5) ) { carbonType = "amorphous";}
  else if  ((compound == "Mol") && (i == 10) ) { carbonType = "graphite";}
    
  if ( carbonType != " ") {  
    // amorphous carbon
    // ICRU - Ziegler are similar below 10 keV
    // the parameter Coeff[5][0] is set to the value from Ziegler 1977
    // taken from Vladimir's work (Ziegler's book is missing at CERN)
    // from 40 to 1000 keV use ICRU coefficients Coeff[5][1] to
    // Coeff[5][5] for Ziegler's function 
    // BUT it should be only up to 600 keV.
    // Linear interpolation between 10 and 40 keV (in log10(E1))
    
    if (KinE < 10.0) {
      ionloss = Coeff[0] * sqrt(KinE) ;
    }
    else if ( KinE <40.0) {
      G4double Elow = Coeff[0] * sqrt(10.) ;
      G4double Slow = Coeff[1] * pow(40., 0.45) ; 
      G4double Shigh = log( 1.0 + Coeff[3]/40. + Coeff[4]*40. ) * Coeff[2]/40.;
      G4double Ehigh = Slow*Shigh / (Slow + Shigh) ; 
      ionloss = (Ehigh-Elow)/(log10(40.)-log10(10.))*(log10(KinE)-log10(10.))+Elow ;
    }
    else if ( KinE <10000.0) {
      G4double Slow  = Coeff[1] * pow(KinE, 0.45) ;
      G4double Shigh = log( 1.0 + Coeff[3]/KinE + Coeff[4]*KinE ) * Coeff[2]/KinE ;
      ionloss = Slow*Shigh / (Slow + Shigh) ; 
    }
    if ( ionloss <= 0.) ionloss = 0. ;
    
    // Graphite is implemented in a very approximate way (scaling 
    // amorphous results according to rough fits to ICRU tables of results:
    // 1-100 keV: *(1+0.023+0.0066*log10(E))
    // 100-700 keV: *(1+0.089-0.0248*log10(E-99.))
    // 700-10000 keV: *(1+0.089-0.0248*log10(700.-99.))
    // continuity is (should!) be garanteed, but not continuity of the
    // first derivative. A better fit is in order!       
    if ( carbonType == "graphite") { 
      if (KinE < .1) {
	ionloss = 0.;
      } 
      else if (KinE < 100.0) {    
	ionloss *= (1+0.023+0.0066*log10(KinE));  
      }
      else if (KinE < 700.0) {   
	ionloss *=(1+0.089-0.0248*log10(KinE-99.));
      } 
      else if (KinE < 10000.0) {    
	ionloss *=(1+0.089-0.0248*log10(700.-99.));
      }     
    }    
    return ionloss;
  }   
  
  if ( KinE < 10.0 ) {
    ionloss = Coeff[0] * sqrt(KinE) ;
    
  } else if ( KinE < 10000.0 ) {
    G4double Slow  = Coeff[1] * pow(KinE, 0.45) ;
    G4double Shigh = log( 1.0 + Coeff[3]/KinE + Coeff[4]*KinE ) * Coeff[2]/KinE ;
    ionloss = Slow*Shigh / (Slow + Shigh) ;     
  } 
  
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977He(G4int iz, G4double KinEnergy)
{
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom
  
  // The He4 data and the fit from: 
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
  
  G4double E = KinEnergy/keV ;  // energy in keV
  
  static G4double A[92][9] = {
    //H-Solid    0.9661,0.4126,  6.92, 8.831,  2.582, 2.371, 0.5462,   -0.07932,-0.006853,
    0.39,0.63,  4.17, 85.55,  19.55, 2.371, 0.5462,   -0.07932,-0.006853,
    //He-Solid   2.027, 0.2931, 26.34, 6.66,   0.3409,2.809, 0.4847,   -0.08756,-0.007281,
    0.58, 0.59, 6.3, 130.0,   44.07, 2.809, 0.4847,   -0.08756,-0.007281,
    1.42,  0.49,   12.25,32.,     9.161, 3.095, 0.4434,   -0.09259,-0.007459,
    2.206, 0.51,   15.32, 0.25,   8.995, 3.28,  0.4188,   -0.09564,-0.007604,
    3.691, 0.4128, 18.48,50.72,   9.,    3.426, 0.4,      -0.09796,-0.007715,
    4.232, 0.3877, 22.99,35.,     7.993, 3.588, 0.3921,   -0.09935,-0.007804,
    //N-solid    2.51,  0.4752, 38.26,13.02,   1.892, 3.759, 0.4094,   -0.09646,-0.007661,
    2.0,  0.548, 29.82,18.11,   4.37, 3.759, 0.4094,   -0.09646,-0.007661,
    //O-solid    1.766, 0.5261, 37.11,15.24,   2.804, 3.782, 0.3734,   -0.1011, -0.007874,
    2.717, 0.4858, 32.88, 25.88, 4.336, 3.782, 0.3734,   -0.1011, -0.007874,
    1.533, 0.531,  40.44,18.41,   2.718, 3.816, 0.3504,   -0.1046, -0.008074,
    //Ne-solid    1.183, 0.55,   39.83,17.49,   4.001, 3.863, 0.3342,   -0.1072, -0.008231,
    2.303, 0.4861,   37.01, 37.96,   5.092, 3.863, 0.3342,   -0.1072, -0.008231,
    9.894, 0.3081, 23.65, 0.384, 92.93,  3.898, 0.3191,   -0.1086, -0.008271,
    4.3,   0.47,   34.3,  3.3,   12.74,  3.961, 0.314,    -0.1091, -0.008297,
    2.5,   0.625,  45.7,  0.1,    4.359, 4.024, 0.3113,   -0.1093, -0.008306,
    2.1,   0.65,   49.34, 1.788,  4.133, 4.077, 0.3074,   -0.1089, -0.008219,
    1.729, 0.6562, 53.41, 2.405,  3.845, 4.124, 0.3023,   -0.1094, -0.00824,
    1.402, 0.6791, 58.98, 3.528,  3.211, 4.164, 0.2964,   -0.1101, -0.008267,
    1.117, 0.7044, 69.69, 3.705,  2.156, 4.21,  0.2936,   -0.1103, -0.00827,
    0.9172,0.724,  79.44, 3.648,  1.646, 4.261, 0.2994,   -0.1085, -0.008145,
    8.554, 0.3817, 83.61,11.84,   1.875, 4.3,   0.2903,   -0.1103, -0.008259,
    6.297, 0.4622, 65.39,10.14,   5.036, 4.334, 0.2897,   -0.1102, -0.008245,
    5.307, 0.4918, 61.74,12.4,    6.665, 4.327, 0.2707,   -0.1127, -0.00837,
    4.71,  0.5087, 65.28, 8.806,  5.948, 4.34,  0.2618,   -0.1138, -0.00842,
    6.151, 0.4524, 83.,  18.31,   2.71,  4.361, 0.2559,   -0.1145, -0.008447,
    6.57,  0.4322, 84.76,15.53,   2.779, 4.349, 0.24,     -0.1166, -0.00855,
    5.738, 0.4492, 84.61,14.18,   3.101, 4.362, 0.2327,   -0.1174, -0.008588,
    5.013, 0.4707, 85.58,16.55,   3.211, 4.375, 0.2253,   -0.1185, -0.008648,
    4.32,  0.4947, 76.14,10.85,   5.441, 4.362, 0.2069,   -0.1214, -0.008815,
    4.652, 0.4571, 80.73,22.,     4.952, 4.346, 0.1857,   -0.1249, -0.009021,
    3.114, 0.5236, 76.67, 7.62,   6.385, 4.355, 0.18,     -0.1255, -0.009045,
    3.114, 0.5236, 76.67, 7.62,   7.502, 4.389, 0.1806,   -0.1253, -0.009028,
    3.114, 0.5236, 76.67, 7.62,   8.514, 4.407, 0.1759,   -0.1258, -0.009054,
    5.746, 0.4662, 79.24, 1.185,  7.993, 4.419, 0.1694,   -0.1267, -0.009094,
    2.792, 0.6346,106.1,  0.2986, 2.331, 4.412, 0.1545,   -0.1289, -0.009202,
    4.667, 0.5095,124.3,  2.102,  1.667, 4.419, 0.1448,   -0.1303, -0.009269,
    2.44,  0.6346,105.,   0.83,   2.851, 4.436, 0.1443,   -0.1299, -0.009229,
    1.491, 0.7118,120.6,  1.101,  1.877, 4.478, 0.1608,   -0.1262, -0.008983,
    11.72, 0.3826,102.8,  9.231,  4.371, 4.489, 0.1517,   -0.1278, -0.009078,
    7.126, 0.4804,119.3,  5.784,  2.454, 4.514, 0.1551,   -0.1268, -0.009005,
    11.61, 0.3955,146.7,  7.031,  1.423, 4.533, 0.1568,   -0.1261, -0.008945,
    10.99, 0.41,  163.9,  7.1,    1.052, 4.548, 0.1572,   -0.1256, -0.008901,
    9.241, 0.4275,163.1,  7.954,  1.102, 4.553, 0.1544,   -0.1255, -0.008883,
    9.276, 0.418, 157.1,  8.038,  1.29,  4.548, 0.1485,   -0.1259, -0.008889,
    3.999, 0.6152, 97.6,  1.297,  5.792, 4.489, 0.1128,   -0.1309, -0.009107,
    4.306, 0.5658, 97.99, 5.514,  5.754, 4.402, 0.06656,  -0.1375, -0.009421,
    3.615, 0.6197, 86.26, 0.333,  8.689, 4.292, 0.01012,  -0.1459, -0.009835,
    5.8,   0.49,  147.2,  6.903,  1.289, 4.187,-0.04539,  -0.1542, -0.01025,
    5.6,   0.49,  130.,  10.,     2.844, 4.577, 0.13,     -0.1285, -0.009067,
    3.55,  0.6068,124.7,  1.112,  3.119, 4.583, 0.1253,   -0.1291, -0.009084,
    3.6,   0.62,  105.8,  0.1692, 6.026, 4.58,  0.1174,   -0.1301, -0.009129,
    5.4,   0.53,  103.1,  3.931,  7.767, 4.581, 0.111,    -0.1309, -0.009161,
    3.97,  0.6459,131.8,  0.2233, 2.723, 4.582, 0.1046,   -0.1317, -0.009193,
    3.65,  0.64,  126.8,  0.6834, 3.411, 4.6,   0.1052,   -0.1315, -0.009178,
    3.118, 0.6519,164.9,  1.208,  1.51,  4.614, 0.1043,   -0.1315, -0.009175,
    2.031, 0.7181,153.1,  1.362,  1.958, 4.619, 0.09769,  -0.1325, -0.009231,
    14.4,  0.3923,152.5,  8.354,  2.597, 4.671, 0.1136,   -0.1298, -0.009078,
    10.99, 0.4599,138.4,  4.811,  3.726, 4.706, 0.1206,   -0.1287, -0.009009,
    16.6,  0.3773,224.1,  6.28,   0.9121,4.732, 0.1244,   -0.128,  -0.008968,
    10.54, 0.4533,159.3,  4.832,  2.529, 4.722, 0.1156,   -0.1292, -0.00903,
    10.33, 0.4502,162.,   5.132,  2.444, 4.71,  0.106,    -0.1305, -0.0091, 
    10.15, 0.4471,165.6,  5.378,  2.328, 4.698, 0.09647,  -0.1319, -0.009169,
    9.976, 0.4439,168.,   5.721,  2.258, 4.681, 0.08536,  -0.1335, -0.009252,
    9.804, 0.4408,176.2,  5.675,  1.997, 4.676, 0.07819,  -0.1345, -0.009302,
    14.22, 0.363, 228.4,  7.024,  1.016, 4.663, 0.06867,  -0.1358, -0.009373,
    9.952, 0.4318,233.5,  5.065,  0.9244,4.676, 0.06861,  -0.1357, -0.009363,
    9.272, 0.4345,210.,   4.911,  1.258, 4.649, 0.05362,  -0.1379, -0.00948,
    10.13, 0.4146,225.7,  5.525,  1.055, 4.634, 0.04335,  -0.1394, -0.009558,
    8.949, 0.4304,213.3,  5.071,  1.221, 4.603, 0.02679,  -0.1418, -0.00969,
    11.94, 0.3783,247.2,  6.655,  0.849, 4.584, 0.01494,  -0.1436, -0.009783,
    8.472, 0.4405,195.5,  4.051,  1.604, 4.576, 0.007043, -0.1447, -0.009841,
    8.301, 0.4399,203.7,  3.667,  1.459, 4.571, 0.0007046,-0.1456, -0.009886,
    6.567, 0.4858,193.,   2.65,   1.66,  4.566,-0.005626, -0.1464, -0.00993,
    5.951, 0.5016,196.1,  2.662,  1.589, 4.561,-0.01197,  -0.1473, -0.009975,
    7.495, 0.4523,251.4,  3.433,  0.8619,4.572,-0.012,    -0.1472, -0.009965,
    6.335, 0.4825,255.1,  2.834,  0.8228,4.569,-0.01755,  -0.148,  -0.01,   
    4.314, 0.5558,214.8,  2.354,  1.263, 4.573,-0.01992,  -0.1482, -0.01001,
    4.02,  0.5681,219.9,  2.402,  1.191, 4.57, -0.02547,  -0.149,  -0.01005,
    3.836, 0.5765,210.2,  2.742,  1.305, 4.528,-0.04613,  -0.1521, -0.01022,
    4.68,  0.5247,244.7,  2.749,  0.8962,4.494,-0.0637,   -0.1548, -0.01037,
    3.223, 0.5883,232.7,  2.954,  1.05,  4.564,-0.027,    -0.1471, -0.009852,
    2.892, 0.6204,208.6,  2.415,  1.416, 4.546,-0.04963,  -0.1523, -0.01022,
    4.728, 0.5522,217.,   3.091,  1.386, 4.594,-0.03339,  -0.1496, -0.01006,
    6.18,  0.52,  170.,   4.,     3.224, 4.608,-0.02886,  -0.1485, -0.00999,
    9.,    0.47,  198.,   3.8,    2.032, 4.624,-0.02639,  -0.1481, -0.009971,
    2.324, 0.6997,216.,   1.599,  1.399, 4.636,-0.02422,  -0.1477, -0.009939,
    1.961, 0.7286,223.,   1.621,  1.296, 4.648,-0.02172,  -0.1471, -0.009903,
    1.75,  0.7427,350.1,  0.9789, 0.5507,4.662,-0.1192,   -0.1752, -0.01196,
    10.31, 0.4613,261.2,  4.738,  0.9899,4.69, -0.009867, -0.1449, -0.009771,
    7.962, 0.519, 235.7,  4.347,  1.313, 4.715,-0.002113, -0.1435, -0.009689,
    6.227, 0.5645,231.9,  3.961,  1.379, 4.729, 0.001392, -0.1428, -0.009644,
    5.246, 0.5947,228.6,  4.027,  1.423, 4.729,-0.0005983,-0.143,  -0.009647,
    5.408, 0.5811,235.7,  3.961,  1.358, 4.738, 0.001075, -0.1425, -0.009618,
    5.218, 0.5828,245.,   3.838,  1.25,  4.751, 0.004244, -0.1419, -0.009576
  };

  if ( E < 1.0 ) {
    G4double Slow  = A[i][0] ;
    G4double Shigh = log( 1.0 + A[i][3]*1000.0 + A[i][4]*0.001 ) 
                   * A[i][2]*1000.0 ;
    ionloss = Slow*Shigh / (Slow + Shigh) ; 
    ionloss *= sqrt(E) ; 
    
  } else if ( E < 10000.0 ) {
    G4double Slow  = A[i][0] * pow(E, A[i][1]) ;
    G4double E1 = E/1000.0 ;
    G4double Shigh = log( 1.0 + A[i][3]/E1 + A[i][4]*E1 ) * A[i][2]/E1 ;
    ionloss = Slow*Shigh / (Slow + Shigh) ; 
    
  } else {
    G4double le = log(1000.0/E) ;
    ionloss = exp( A[i][5] + A[i][6]*le + A[i][7]*le*le + A[i][8]*le*le*le) ;
  }
  
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPowerICRU_R49PowersHe(G4int iz, G4double KinEnergy)
{
  // ********* NOT IMPLEMENTED YET                         *************** 
  // ********* ICRU USES ZIEGLER FOR SOME PARAMETRISATIONS *************** 
  // ********* AND POWERS FOR OTHER MATERIALS              ***************
  // ********* TO BE DONE!                                 *************** 
  // Powers' parametrisations
  
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom
  
  // The data and the fit from: 
  // ICRU Report 49, 1993
  
  G4double KinE = KinEnergy/keV ;  // energy in keV
  
  static G4double A[30][7] = {
      8.0080,   3.6287,  23.0700,  14.9900,  0.8507, 0.60, 2.0
    , 13.3100,  3.7432,  39.4130,  12.1990,  1.0950, 0.38, 1.4
    , 22.7240,  3.6040,  47.1810,  17.5490,  0.9040, 0.40, 1.4
    , 24.4040,  2.4032,  48.9440,  27.9730,  1.2933, 0.40, 1.6
    , 58.4719,  1.5115,  77.6421,  102.490,  1.5811, 0.50, 2.0
    , 60.5408,  1.6297,  91.7601,  94.1260,  1.3662, 0.50, 2.0
    , 48.4480,  6.4323,  59.2890,  18.3810,  0.4937, 0.48, 1.6
    , 59.0346,  5.1305,  47.0866,  30.0857,  0.3500, 0.60, 2.0
    , 71.8691,  2.8250,  51.1658,  57.1235,  0.4477, 0.60, 2.0
    , 78.3520,  4.0961,  136.731,  28.4470,  1.0621, 0.52, 1.2
    , 120.553,  1.5374,  49.8740,  82.2980,  0.8733, 0.45, 1.6
    , 249.896,  0.6996,  -37.274,  248.592,  1.1052, 0.50, 1.5
    , 246.698,  0.6219,  -58.391,  292.921,  0.8186, 0.56, 1.8
    , 248.563,  0.6235,  -36.8968, 306.960,  1.3214, 0.50, 2.0
    , 25.5860,  1.7125,  154.723,  118.620,  2.2580, 0.50, 2.0
    , 138.294,  25.6413, 231.873,  17.3780,  0.3218, 0.58, 1.3
    , 83.2091,  1.1294,  135.7457, 190.865,  2.3461, 0.50, 2.0
    , 263.542,  1.4754,  1541.446, 781.898,  1.9209, 0.40, 2.0
    , 59.5545,  1.5354,  132.1523, 153.3537, 2.0262, 0.50, 2.0
    , 31.7380,  19.820,  125.2100, 6.8910,   0.7242, 0.50, 1.1
    , 31.7549,  1.5682,  97.4777,  106.0774, 2.3204, 0.50, 2.0
    , 230.465,  4.8967,  1845.320, 358.641,  1.0774, 0.46, 1.2
    , 423.444,  5.3761,  1189.114, 319.030,  0.7652, 0.48, 1.5
    , 86.3410,  3.3322,  91.0433,  73.1091,  0.4650, 0.50, 2.0
    , 146.105,  9.4344,  515.1500, 82.8860,  0.6239, 0.55, 1.5
    , 238.050,  5.6901,  372.3575, 146.1835, 0.3992, 0.50, 2.0
    , 124.2338, 2.6730,  133.8175, 99.4109,  0.7776, 0.50, 2.0
    , 221.723,  1.5415,  87.7315,  192.5266, 1.0742, 0.50, 2.0
    , 26.7537,  1.3717,  90.8007,  77.1587,  2.3264, 0.50, 2.0
    , 37.6121,  1.8052,  73.0250,  66.2070,  1.4038, 0.50, 2.0 
    
  };
  
  if ( KinE < 20000.0 ) {
    G4double E1 = KinE/1000.0 ;
    G4double A1 = 1-exp(-A[i][1]*pow(E1,-2+A[i][5])) ;
    G4double A2 = (A[i][0]*log(E1)/E1+A[i][2]/E1)*exp(-A[i][4]*pow(E1,-A[i][6]))+A[i][3]/(E1*E1) ;
    ionloss = A1*A2 ; 
    
  }
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPowerICRU_R49He(G4int iz, G4double KinEnergy)
{
  
  G4double ionloss ;
  G4int i = iz-1 ;  // index of atom
  
  // The data and the fit from: 
  // ICRU Report 49, 1993
  // Ziegler's parametrisations
  
  G4double KinE = KinEnergy/MeV ;  // energy in MeV
  
  static G4double A[92][5] = {
      0.35485, 0.6456, 6.01525,  20.8933, 4.3515
    , 0.58,    0.59,   6.3,	 130.0,   44.07
    , 1.42,    0.49,   12.25,    32.0,    9.161
    , 2.1895,  0.47183,7.2362,   134.30,  197.96
    , 3.691,   0.4128, 18.48,    50.72,   9.0
    , 3.83523, 0.42993,12.6125,  227.41,  188.97
    , 1.9259,  0.5550, 27.15125, 26.0665, 6.2768
    , 2.81015, 0.4759, 50.0253,  10.556,  1.0382
    , 1.533,   0.531,  40.44,    18.41,   2.718
    , 2.303,   0.4861, 37.01,    37.96,   5.092
    , 9.894,   0.3081, 23.65,    0.384,   92.93
    , 4.3,     0.47,   34.3,     3.3,     12.74
    , 2.5,     0.625,  45.7,     0.1,     4.359
    , 2.1,     0.65,   49.34,    1.788,   4.133
    , 1.729,   0.6562, 53.41,    2.405,   3.845
    , 1.402,   0.6791, 58.98,    3.528,   3.211
    , 1.117,   0.7044, 69.69,    3.705,    2.156
    , 2.291,   0.6284, 73.88,    4.478,    2.066
    , 8.554,   0.3817, 83.61,    11.84,    1.875
    , 6.297,   0.4622, 65.39,    10.14,    5.036
    , 5.307,   0.4918, 61.74,    12.4,	   6.665
    , 4.71,    0.5087, 65.28,    8.806,    5.948
    , 6.151,   0.4524, 83.0,	 18.31,    2.71
    , 6.57,    0.4322, 84.76,    15.53,    2.779
    , 5.738,   0.4492, 84.6,	 14.18,    3.101
    , 5.013,   0.4707, 85.8,	 16.55,    3.211
    , 4.32,    0.4947, 76.14,    10.85,    5.441
    , 4.652,   0.4571, 80.73,    22.0,	   4.952
    , 3.114,   0.5236, 76.67,    7.62,	   6.385
    , 3.114,   0.5236, 76.67,    7.62,	   7.502
    , 3.114,   0.5236, 76.67,    7.62,	   8.514
    , 5.746,   0.4662, 79.24,    1.185,    7.993
    , 2.792,   0.6346, 106.1,    0.2986,   2.331
    , 4.667,   0.5095, 124.3,    2.102,    1.667
    , 2.44,    0.6346, 105.0,    0.83,	   2.851
    , 1.413,   0.7377, 147.9,    1.466,    1.016
    , 11.72,   0.3826, 102.8,    9.231,    4.371
    , 7.126,   0.4804, 119.3,    5.784,    2.454
    , 11.61,   0.3955, 146.7,    7.031,    1.423
    , 10.99,   0.41,   163.9,	 7.1,	   1.052
    , 9.241,   0.4275, 163.1,    7.954,    1.102
    , 9.276,   0.418,  157.1,	 8.038,    1.29
    , 3.999,   0.6152, 97.6,	 1.297,    5.792
    , 4.306,   0.5658, 97.99,    5.514,    5.754
    , 3.615,   0.6197, 86.26,    0.333,    8.689
    , 5.8,     0.49,   147.2,	 6.903,    1.289
    , 5.6,     0.49,   130.0,	 10.0,	   2.844
    , 3.55,    0.6068, 124.7,    1.112,    3.119
    , 3.6,     0.62,   105.8,	 0.1692,   6.026
    , 5.4,     0.53,   103.1,	 3.931,    7.767
    , 3.97,    0.6459, 131.8,    0.2233,   2.723 
    , 3.65,    0.64,   126.8,	 0.6834,   3.411
    , 3.118,   0.6519, 164.9,    1.208,    1.51
    , 3.949,   0.6209, 200.5,    1.878,    0.9126
    , 14.4,    0.3923, 152.5,    8.354,    2.597
    , 10.99,   0.4599, 138.4,    4.811,    3.726
    , 16.6,    0.3773, 224.1,    6.28,	   0.9121 
    , 10.54,   0.4533, 159.3,	 4.832,	   2.529
    , 10.33,   0.4502, 162.0,	 5.132,	   2.444
    , 10.15,   0.4471, 165.6,	 5.378,	   2.328
    , 9.976,   0.4439, 168.0,	 5.721,	   2.258
    , 9.804,   0.4408, 176.2,	 5.675,	   1.997
    , 14.22,   0.363,  228.4,	 7.024,	   1.016
    , 9.952,   0.4318, 233.5,	 5.065,	   0.9244
    , 9.272,   0.4345, 210.0,	 4.911,	   1.258
    , 10.13,   0.4146, 225.7,	 5.525,	   1.055
    , 8.949,   0.4304, 213.3,	 5.071,	   1.221
    , 11.94,   0.3783, 247.2,	 6.655,	   0.849
    , 8.472,   0.4405, 195.5,	 4.051,	   1.604
    , 8.301,   0.4399, 203.7,	 3.667,	   1.459
    , 6.567,   0.4858, 193.0,	 2.65,	   1.66
    , 5.951,   0.5016, 196.1,	 2.662,	   1.589
    , 7.495,   0.4523, 251.4,	 3.433,	   0.8619
    , 6.335,   0.4825, 255.1,	 2.834,	   0.8228
    , 4.314,   0.5558, 214.8,	 2.354,	   1.263
    , 4.02,    0.5681, 219.9,	 2.402,	   1.191
    , 3.836,   0.5765, 210.2,	 2.742,	   1.305
    , 4.68,    0.5247, 244.7,	 2.749,	   0.8962
    , 3.223,   0.5883, 232.7,	 2.954,	   1.05
    , 2.892,   0.6204, 208.6,	 2.415,	   1.416
    , 4.728,   0.5522, 217.0,	 3.091,	   1.386
    , 6.18,    0.52,   170.0,	 4.0,	   3.224
    , 9.0,     0.47,   198.0,	 3.8,	   2.032
    , 2.324,   0.6997, 216.0,	 1.599,	   1.399
    , 1.961,   0.7286, 223.0,	 1.621,	   1.296
    , 1.75,    0.7427, 350.1,	 0.9789,   0.5507
    , 10.31,   0.4613, 261.2,	 4.738,	   0.9899
    , 7.962,   0.519,  235.7,	 4.347,	   1.313
    , 6.227,   0.5645, 231.9,	 3.961,	   1.379
    , 5.246,   0.5947, 228.6,	 4.027,	   1.432
    , 5.408,   0.5811, 235.7,	 3.961,	   1.358
    , 5.218,   0.5828, 245.0,	 3.838,	   1.25	
  };
  
   if ( KinE < 0.001 ) {
    G4double Slow  = A[i][0] ;
    G4double Shigh = log( 1.0 + A[i][3]*1000.0 + A[i][4]*0.001 ) 
                   * A[i][2]*1000.0 ;
    ionloss = Slow*Shigh / (Slow + Shigh) ; 
    ionloss *= sqrt(KinE*1000.0) ; 
    
  } else {
    G4double Slow  = A[i][0] * pow((KinE*1000.0), A[i][1]) ;
    G4double Shigh = log( 1.0 + A[i][3]/KinE + A[i][4]*KinE ) * A[i][2]/KinE ;
    ionloss = Slow*Shigh / (Slow + Shigh) ; 
    
  }
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1977n(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double KinEnergy)
{
  // The fit of nuclear stopping from: 
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
  
  G4double E = KinEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (M1 + M2) * sqrt( pow(Z1, 0.667) + pow(Z2, 0.667) ) ;
  
  G4double er = 32.53 * M2 * E / ( Z1 * Z2 * rm ) ;  // reduced energy
  
  if ( er < 0.01 ) {
    ionloss = sqrt(er) * 1.593 ; 
    
  } else if ( E < 10.0 ) {
    ionloss = 1.7 * sqrt(er) * log(er + exp(1.0)) / 
      (1.0 + 6.8 * er + 3.4 * pow(er, 1.5)) ; 
    
  } else {
    ionloss = log(0.47 * er) * 0.5 / er  ;
  }
  
  ionloss *= 8.462 * Z1 * Z2 * M1 / rm ; // Return to [ev/(10^15 atoms/cm^2]
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPower1985n(G4double Z1, G4double Z2, 
						       G4double M1, G4double M2, G4double KinEnergy)
{
  // The fit of nuclear stopping from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  
  G4double E = KinEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  G4double rm = (M1 + M2) * sqrt( pow(Z1, .23) + pow(Z2, .23) ) ;
  
  G4double er = 32.536 * M2 * E / ( Z1 * Z2 * rm ) ;  // reduced energy
  
  if ( er <= 30 ) {
    ionloss = 0.5*log(1+1.1383*er)/(er+0.01312*pow(er,0.21226)+0.19593*pow(er,0.5)) ; 
    
  } else {
    ionloss = 0.5*log(er)/er ; 
  }
  
  ionloss *= 8.462 * Z1 * Z2 * M1 / rm ; // Return to [ev/(10^15 atoms/cm^2]
  
  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetStoppingPowerMoliere(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double KinEnergy)
{
  // The fit of nuclear stopping from: 
  // G.Moliere "Theorie der Streuung schneller geladener Teilchen I;
  // Einzelstreuungam abbgeschirmten Coulomb-Feld" Z. f. Naturforsch, A2,
  // 133 (1947).
  // See also ICRU Report 49, 1993 
  
  G4double E = KinEnergy/keV ;  // energy in keV
  G4double ionloss ;
  
  static G4double A[104][2] = {
    1.0E+8,	5.831E-8,
    8.0E+7,	7.288E-8,
    6.0E+7,	9.719E-8,
    5.0E+7,	1.166E-7,
    4.0E+7,	1.457E-7,
    3.0E+7,	1.942E-7,
    2.0E+7,	2.916E-7,
    1.5E+7,	3.887E-7,
    1.0E+7,	5.833E-7,
    8.0E+6,	7.287E-7,
    6.0E+6,	9.712E-7,
    5.0E+6,	1.166E-6,
    4.0E+6,	1.457E-6,
    3.0E+6,	1.941E-6,
    2.0E+6,	2.911E-6,
    1.5E+6,	3.878E-6,
    1.0E+6,	5.810E-6,
    8.0E+5,	7.262E-6,
    6.0E+5,	9.663E-6,
    5.0E+5,	1.157E-5,
    4.0E+5,	1.442E-5,
    3.0E+5,	1.913E-5,
    2.0E+5,	2.845E-5,
    1.5E+5,	3.762E-5,
    1.0E+5,	5.554E-5,
    8.0E+4,	6.866E-5,
    6.0E+4,	9.020E-5,
    5.0E+4,	1.070E-4,
    4.0E+4,	1.319E-4,
    3.0E+4,	1.722E-4,
    2.0E+4,	2.499E-4,
    1.5E+4,	3.248E-4,
    1.0E+4,	4.688E-4,
    8.0E+3,	5.729E-4,
    6.0E+3,	7.411E-4,
    5.0E+3,	8.718E-4,
    4.0E+3,	1.063E-3,
    3.0E+3,	1.370E-3,
    2.0E+3,	1.955E-3,
    1.5E+3,	2.511E-3,
    1.0E+1,	1.210E-1,
    8.0E+0,	1.377E-1,
    6.0E+0,	1.611E-1,
    5.0E+0,	1.768E-1,
    4.0E+0,	1.968E-1,
    3.0E+0,	2.235E-1,
    2.0E+0,	2.613E-1,
    1.5E+0,	2.871E-1,
    1.0E+0,	3.199E-1,
    8.0E-1,	3.354E-1,
    6.0E-1,	3.523E-1,
    5.0E-1,	3.609E-1,
    4.0E-1,	3.693E-1,
    3.0E-1,	3.766E-1,
    2.0E-1,	3.803E-1,
    1.5E-1,	3.788E-1,
    1.0E-1,	3.711E-1,
    8.0E-2,	3.644E-1,
    6.0E-2,	3.530E-1,
    5.0E-2,	3.444E-1,
    4.0E-2,	3.323E-1,
    3.0E-2,	3.144E-1,
    2.0E-2,	2.854E-1,
    1.5E-2,	2.629E-1,
    1.0E-2,	2.298E-1,
    8.0E-3,	2.115E-1,
    6.0E-3,	1.883E-1,
    5.0E-3,	1.741E-1,
    4.0E-3,	1.574E-1,
    3.0E-3,	1.372E-1,
    2.0E-3,	1.116E-1,
    1.5E-3,	9.559E-2,
    1.0E-3,	7.601E-2,
    8.0E-4,	6.668E-2,
    6.0E-4,	5.605E-2,
    5.0E-4,	5.008E-2,
    4.0E-4,	4.352E-2,
    3.0E-4,	3.617E-2,
    2.0E-4,	2.768E-2,
    1.5E-4,	2.279E-2,
    1.0E+3,	3.563E-3,
    8.0E+2,	4.314E-3,
    6.0E+2,	5.511E-3,
    5.0E+2,	6.430E-3,
    4.0E+2,	7.756E-3,
    3.0E+2,	9.855E-3,
    2.0E+2,	1.375E-2,
    1.5E+2,	1.736E-2,
    1.0E+2,	2.395E-2,
    8.0E+1,	2.850E-2,
    6.0E+1,	3.552E-2,
    5.0E+1,	4.073E-2,
    4.0E+1,	4.802E-2,
    3.0E+1,	5.904E-2,
    1.5E+1,	9.426E-2,
    1.0E-4,	1.723E-2,
    8.0E-5,	1.473E-2,
    6.0E-5,	1.200E-2,
    5.0E-5,	1.052E-2,
    4.0E-5,	8.950E-3,
    3.0E-5,	7.246E-3,
    2.0E-5,	5.358E-3,
    1.5E-5,	4.313E-3,
    1.0E-5,	3.166E-3
  };
  
  G4double rm = (M1 + M2) * sqrt( pow(Z1, .23) + pow(Z2, .23) ) ;
  
  G4double er = 32.536 * M2 * E / ( Z1 * Z2 * rm ) ;  // reduced energy
  
  for (G4int i=0; 104; i++)
    {
      if (er > A[i][0]) {
	ionloss =
	  (er-A[i][0])/(A[i-1][0]-A[i][0])*(A[i-1][1]-A[i][1])+A[i][1];   
	break;
      }
    }
  
  ionloss *= 8.462 * Z1 * Z2 * M1 / rm ; // Return to [ev/(10^15 atoms/cm^2]

  if ( ionloss <= 0.) ionloss = 0. ;
  
  return ionloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::MolecIsInZiegler1988(const G4Material* material) 

{
  // The list of molecules from
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

  // If the meterial is in the table then the Stopping Power at 125 keV exist
  // In that case the return value ExpStopPower125 > 0
  G4double ExpStopPower125 = -1.0;
  
  const G4String chFormula = material->GetChemicalFormula() ;
  if (" " == chFormula ) return ExpStopPower125 ;
  
  //  There are no evidence for difference of stopping power depended on
  //  phase of the compound except for water. The stopping power of the 
  //  water in gas phase can be predicted using Bragg's rule.
  //  
  //  Do not specify chemical formula for water-gas!!! 
   
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && "H_2O" == chFormula) return ExpStopPower125 ;
    
  const size_t NumberOfMolecula = 53 ;
  const G4double HeEff = 2.8735 ;    // The coffecient from Table.4 of Ziegler & Manoyan
  
  static G4String Name[NumberOfMolecula] = {
    "H_2O",      "C_2H_4O",    "C_3H_6O",  "C_2H_2",             "C_H_3OH",
    "C_2H_5OH",  "C_3H_7OH",   "C_3H_4",   "NH_3",               "C_14H_10",
    "C_6H_6",    "C_4H_10",    "C_4H_6",   "C_4H_8O",            "CCl_4",
    "CF_4",      "C_6H_8",     "C_6H_12",  "C_6H_10O",           "C_6H_10",
    "C_8H_16",   "C_5H_10",    "C_5H_8",   "C_3H_6-Cyclopropane","C_2H_4F_2",
    "C_2H_2F_2", "C_4H_8O_2",  "C_2H_6",   "C_2F_6",             "C_2H_6O",
    "C_3H_6O",   "C_4H_10O",   "C_2H_4",   "C_2H_4O",            "C_2H_4S",
    "SH_2",      "CH_4",       "CCLF_3",   "CCl_2F_2",           "CHCl_2F",
    "(CH_3)_2S", "N_2O",       "C_5H_10O", "C_8H_6",             "(CH_2)_N",
    "(C_3H_6)_N","(C_8H_8)_N", "C_3H_8",   "C_3H_6-Propylene",   "C_3H_6O",
    "C_3H_6S",   "C_4H_4S",    "C_7H_8"
  } ;
    
  static G4double ExpStopping[NumberOfMolecula] = {
   66.1,  190.4, 258.7,  42.2, 141.5, 
  210.9,  279.6, 198.8,  31.0, 267.5,
  122.8,  311.4, 260.3, 328.9, 391.3,
  206.6,  374.0, 422.0, 432.0, 398.0,
  554.0,  353.0, 326.0,  74.6, 220.5,
  197.4,  362.0, 170.0, 330.5, 211.3,
  262.3,  349.6,  51.3, 187.0, 236.9,
  121.9,   35.8, 247.0, 292.6, 268.0,
  262.3,   49.0, 398.9, 444.0,  22.91,
   68.0,  155.0,  84.0,  74.2, 254.7,
  306.8,  324.4, 420.0
  } ;

  static G4double ExpCharge[NumberOfMolecula] = {
  HeEff, HeEff, HeEff,   1.0, HeEff, 
  HeEff, HeEff, HeEff,   1.0,   1.0,
    1.0, HeEff, HeEff, HeEff, HeEff,
  HeEff, HeEff, HeEff, HeEff, HeEff,
  HeEff, HeEff, HeEff,   1.0, HeEff,
  HeEff, HeEff, HeEff, HeEff, HeEff,
  HeEff, HeEff,   1.0, HeEff, HeEff,
  HeEff,   1.0, HeEff, HeEff, HeEff,
  HeEff,   1.0, HeEff, HeEff,   1.0,
    1.0,   1.0,   1.0,   1.0, HeEff,
  HeEff, HeEff, HeEff
  } ;

  static G4double NumberOfAtomsPerMolecula[NumberOfMolecula] = {
   3.0,  7.0, 10.0,  4.0,  6.0,  
   9.0, 12.0,  7.0,  4.0, 24.0,
  12.0, 14.0, 10.0, 13.0,  5.0,
   5.0, 14.0, 18.0, 17.0, 17.0,
  24.0, 15.0, 13.0,  9.0,  8.0,
   6.0, 14.0,  8.0,  8.0,  9.0,
  10.0, 15.0,  6.0,  7.0,  7.0,
   3.0,  5.0,  5.0,  5.0,  5.0,
   9.0,  3.0, 16.0, 14.0,  3.0,
   9.0, 16.0, 11.0,  9.0, 10.0,
  10.0,  9.0, 15.0
  } ;

  // Search for the compaund in the table
  for (G4int i=0; i<NumberOfMolecula; i++)
    { 
      if(chFormula == Name[i]) { 
        ExpStopPower125 = ExpStopping[i] * 
                          (material-> GetTotNbOfAtomsPerVolume()) *
                          ZieglerFactor / 
                          (ExpCharge[i] * NumberOfAtomsPerMolecula[i]) ;
        return ExpStopPower125 ;
      }
    }
  
  return ExpStopPower125;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetChemicalFactor(const G4double ExpStopPower125, 
                                                   const G4double KinEnergy,
                                                   const G4double BraggStopPower125)
{
  // Approximation of Chemical Factor according to
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
  
  G4double gamma    = 1.0 + KinEnergy/proton_mass_c2 ;    
  G4double gamma25  = 1.0 + 25.0*keV /proton_mass_c2 ;
  G4double gamma125 = 1.0 + 125.0*keV/proton_mass_c2 ;
  G4double beta     = sqrt(1.0 - 1.0/(gamma*gamma)) ;
  G4double beta25   = sqrt(1.0 - 1.0/(gamma25*gamma25)) ;
  G4double beta125  = sqrt(1.0 - 1.0/(gamma125*gamma125)) ;
  
  G4double factor = 1.0 + (ExpStopPower125/BraggStopPower125 - 1.0) *
       (1.0 + exp( 1.48 * ( beta125/beta25 - 7.0 ) ) ) /
       (1.0 + exp( 1.48 * ( beta/beta25    - 7.0 ) ) ) ;
  
  return factor ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetHeEffChargeSquare(const G4int iz,
                                                      const G4double HeKinEnergy)

{
  // The aproximation of He effective charge from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  static G4double C[6] = {0.2865,  0.1266, -0.001429,
                          0.02402,-0.01135, 0.001475} ;

  G4double E = log( G4std::max( 1.0, HeKinEnergy/(keV*HeMassAMU) ) ) ; 
  G4double x = C[0] ;
  G4double y = 1.0 ;
    for (G4int i=1; i<6; i++) {
      y *= E ;
      x += y * C[i] ;
    }
  G4double z = 7.6 -  E ;
  z = 1.0 + (0.007 + 0.00005*iz) * exp( -z*z ) ;
 
  G4double w = 4.0 * (1.0 - exp(-x)) * z * z ;

  return w ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetIonEffChargeSquare(const G4Material* aMaterial,
                                                       const G4double KinEnergy,
                                                       const G4double IonCharge)
 
{
  // The aproximation of ion effective charge from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  // Fast ions or hadrons
  G4double ReducedEnergy = KinEnergy * MassRatio ;
  if( (ReducedEnergy > IonCharge * 10.0 * MeV) || 
      (IonCharge < 1.5) ) return IonCharge*IonCharge ;

  static G4double VFermi[92] = {
   1.0309,  0.15976, 0.59782, 1.0781,  1.0486,  1.0,     1.058,   0.93942, 0.74562, 0.3424,
   0.45259, 0.71074, 0.90519, 0.97411, 0.97184, 0.89852, 0.70827, 0.39816, 0.36552, 0.62712,
   0.81707, 0.9943,  1.1423,  1.2381,  1.1222,  0.92705, 1.0047,  1.2,     1.0661,  0.97411,
   0.84912, 0.95,    1.0903,  1.0429,  0.49715, 0.37755, 0.35211, 0.57801, 0.77773, 1.0207,
   1.029,   1.2542,  1.122,   1.1241,  1.0882,  1.2709,  1.2542,  0.90094, 0.74093, 0.86054,
   0.93155, 1.0047,  0.55379, 0.43289, 0.32636, 0.5131,  0.695,   0.72591, 0.71202, 0.67413,
   0.71418, 0.71453, 0.5911,  0.70263, 0.68049, 0.68203, 0.68121, 0.68532, 0.68715, 0.61884,
   0.71801, 0.83048, 1.1222,  1.2381,  1.045,   1.0733,  1.0953,  1.2381,  1.2879,  0.78654,
   0.66401, 0.84912, 0.88433, 0.80746, 0.43357, 0.41923, 0.43638, 0.51464, 0.73087, 0.81065,
   1.9578,  1.0257} ;

  static G4double LFactor[92] = {
   1.0,  1.0,  1.1,  1.06, 1.01, 1.03, 1.04, 0.99, 0.95, 0.9,
   0.82, 0.81, 0.83, 0.88, 1.0,  0.95, 0.97, 0.99, 0.98, 0.97,
   0.98, 0.97, 0.96, 0.93, 0.91, 0.9,  0.88, 0.9,  0.9,  0.9,
   0.9,  0.85, 0.9,  0.9,  0.91, 0.92, 0.9,  0.9,  0.9,  0.9,
   0.9,  0.88, 0.9,  0.88, 0.88, 0.9,  0.9,  0.88, 0.9,  0.9,
   0.9,  0.9,  0.96, 1.2,  0.9,  0.88, 0.88, 0.85, 0.9,  0.9, 
   0.92, 0.95, 0.99, 1.03, 1.05, 1.07, 1.08, 1.1,  1.08, 1.08,
   1.08, 1.08, 1.09, 1.09, 1.1,  1.11, 1.12, 1.13, 1.14, 1.15,
   1.17, 1.2,  1.18, 1.17, 1.17, 1.16, 1.16, 1.16, 1.16, 1.16,
   1.16, 1.16} ; 

  static G4double C[6] = {0.2865,  0.1266, -0.001429,
                          0.02402,-0.01135, 0.001475} ;

  // get elements in the actual material,
  const G4ElementVector* theElementVector = aMaterial->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = aMaterial->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = aMaterial->GetNumberOfElements() ;
  
  //  loop for the elements in the material
  //  to find out average values Z, vF, lF
  G4double Z = 0.0, vF =0.0, lF = 0.0 , Norm = 0.0 ; 

  if( 1 == NumberOfElements ) {
    Z = aMaterial->GetZ() ;
    G4int iz = G4int(Z) - 1 ;
    if(iz < 0) iz = 0 ;
    else if(iz > 91) iz =91 ;
    vF   = VFermi[iz] ;
    lF   = LFactor[iz] ;

  } else {
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        const G4Element* element = (*theElementVector)(iel) ;
        G4double Z2 = element->GetZ() ;
        const G4double W2 = theAtomicNumDensityVector[iel] ;
        Norm += W2 ;
        Z    += Z2 * W2 ;
        G4int iz = G4int(Z2) - 1 ;
        if(iz < 0) iz = 0 ;
        else if(iz > 91) iz =91 ;
        vF   += VFermi[iz] * W2 ;
        lF   += LFactor[iz] * W2 ;
      }
    Z  /= Norm ;
    vF /= Norm ;
    lF /= Norm ;
  }

  G4double w, Q ;

  // Helium ion case
  if( IonCharge < 2.5 ) {

    G4double E = log( G4std::max( 1.0, KinEnergy / (keV*HeMassAMU) ) ) ; 
    G4double x = C[0] ;
    G4double y = 1.0 ;
      for (G4int i=1; i<6; i++) {
        y *= E ;
        x += y * C[i] ;
      }
    Q = 7.6 -  E ; 
    Q = 1.0 + ( 0.007 + 0.00005 * Z ) * exp( -Q*Q ) ;
    return  4.0 * Q * Q * (1.0 - exp(-x)) ;

  // Heavy ion case
  } else {

    // v1 is ion velocity in vF unit
    G4double v1 = sqrt( ReducedEnergy / (25.0 * keV) )/ vF ;
    G4double y ;
    G4double Z13 = pow(IonCharge, 0.3333) ;

    // Faster than Fermi velocity
    if ( v1 > 1.0 ) {
      y = vF * v1 * ( 1.0 + 0.2 / (v1*v1) ) / (Z13*Z13) ;

    // Slower than Fermi velocity
    } else {
      y = 0.75 * vF * (1.0 + 2.0*v1*v1/3.0 + v1*v1*v1*v1/15.0) / (Z13*Z13) ;
    }

    G4double y3 = pow(y, 0.3) ;
    G4double q = 1.0 - exp( 0.803*y3 - 1.3167*y3*y3 - 0.38157*y - 0.008983*y*y ) ;     
    if( q < 0.0 ) q = 0.0 ;

    Q = 7.6 -  log(G4std::max(1.0, ReducedEnergy/keV)) ; 
    Q = 1.0 + ( 0.18 + 0.0015 * Z ) * exp( -Q*Q )/ (IonCharge*IonCharge) ;

    // Screen length according to
    // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
    // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

    G4double Lambda = 10.0 * vF * pow(1.0-q, 0.6667) / (Z13 * (6.0 + q)) ;
    G4double Qeff   = IonCharge * Q *
                      ( q + 0.5*(1.0-q) * log(1.0 + Lambda*Lambda) / (vF*vF) ) ;
    if( 1.0 > Qeff ) Qeff = 1.0 ; 
    return Qeff*Qeff ;    
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n         Good description above the mean excitation energy.\n";
  comments += "         delta ray energy sampled from  differential Xsection.";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << LowestKineticEnergy / eV << " eV " 
         << " to " << HighestKineticEnergy / TeV << " TeV "
         << " in " << TotBin << " bins."
         << "\n        Low energy losses approximation is taken from  " << DEDXtable
         << "\n        from " << ParamLowEnergy / keV << " keV "
         << " to " << ParamHighEnergy / MeV << " MeV " << "." << G4endl ;
  if(pbarStop){
    G4cout << "        Parametrization of the Barkas effect is switched on." << G4endl ;
    G4cout << "        Quantal Harmonic Oscillator Model is switched on " << G4endl ;
  }
  if(nStopping) {
    G4cout << "        Simulation of nuclear stopping is switched on." << G4endl ; 
  }
}











































