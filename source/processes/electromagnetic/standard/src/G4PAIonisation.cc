// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIonisation.cc,v 1.3 1999-12-15 14:51:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4PAIonisation physics process -----------
//                 modified by V.Grichine 27.11.97 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of charged hadrons.
// **************************************************************
//
// 08-04-98: remove 'traking cut' of the ionizing particle, MMa
// 30-11-97: V. Grichine
// 
 

#include "G4PAIonisation.hh"
#include "G4PAIxSection.hh"

const G4double G4PAIonisation:: LowestKineticEnergy = 100.0*MeV  ;
const G4double G4PAIonisation::HighestKineticEnergy = 10.*TeV  ;
G4int G4PAIonisation::TotBin = 100  ;  // 50

      // create physics vector and fill it
      
G4PhysicsLogVector* 
G4PAIonisation::fProtonEnergyVector = new G4PhysicsLogVector(LowestKineticEnergy,
							   HighestKineticEnergy,
							   TotBin);
  



//////////////////////////////////////////////////////////////////////////////
//
// constructor and destructor
//
 
G4PAIonisation::G4PAIonisation( const G4String& materialName,
                                const G4String& processName)
   : G4PAIenergyLoss(processName),
     theElectron ( G4Electron::Electron() )
{
  G4int numberOfMat, iMat ;
  theMeanFreePathTable  = NULL; 
  lastCutInRange = 0. ;
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  numberOfMat = theMaterialTable->length() ;
  for(iMat=0;iMat<numberOfMat;iMat++)
  {
    if(materialName == (*theMaterialTable)[iMat]->GetName() )
    {
      fMatIndex = (*theMaterialTable)[iMat]->GetIndex() ;
      break ;
    }
  }
  if(iMat == numberOfMat)
  {
    G4Exception("Invalid material name in G4PAIonisation constructor") ;
  }
  ComputeSandiaPhotoAbsCof() ;


  // G4cout<<"G4PAIonisation constructor is called"<<G4endl ;
  // BuildPAIonisationTable() ;
}

///////////////////////////////////////////////////////////////////////////
//
//

G4PAIonisation::~G4PAIonisation() 
{
     if (theMeanFreePathTable)
     {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }

}
 
/////////////////////////////////////////////////////////////////////////
//
//


void G4PAIonisation::ComputeSandiaPhotoAbsCof()
{
   G4int i, j, numberOfElements ;

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // 

   G4SandiaTable thisMaterialSandiaTable(fMatIndex) ;
   numberOfElements = (*theMaterialTable)[fMatIndex]->
                                              GetNumberOfElements() ;
   G4int* thisMaterialZ = new G4int[numberOfElements] ;
   for(i=0;i<numberOfElements;i++)
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[fMatIndex]->
                                      GetElement(i)->GetZ() ;
   }
   fSandiaIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
   fSandiaIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                             (*theMaterialTable)[fMatIndex]->GetFractionVector() ,
        		     numberOfElements,fSandiaIntervalNumber) ;
   
   fSandiaPhotoAbsCof = new G4double*[fSandiaIntervalNumber] ;

   for(i=0;i<fSandiaIntervalNumber;i++)
   {
     fSandiaPhotoAbsCof[i] = new G4double[5] ;
   }
   for(i=0;i<fSandiaIntervalNumber;i++)
   {
      fSandiaPhotoAbsCof[i][0] = thisMaterialSandiaTable.
                                  GetPhotoAbsorpCof(i+1,0) ; // keV ;

                                               // G4double energyCof = keV ;
      for(j=1;j<5;j++)
      {
           fSandiaPhotoAbsCof[i][j] = thisMaterialSandiaTable.
	                              GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex]->GetDensity() ;
	    // *(cm2/g)*energyCof ;
	    // energyCof *= keV ;
      }
   }
   delete[] thisMaterialZ ;
}


 
////////////////////////////////////////////////////////////////////////
//
//  just call BuildLossTable+BuildLambdaTable
//

void 
G4PAIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)

{
    G4double Charge = aParticleType.GetPDGCharge();
    G4double Chargesquare = Charge*Charge ;     
    CutInRange = aParticleType.GetLengthCuts(); 

       BuildLossTable(aParticleType) ;
 
    if(Charge>0.)
    {
       RecorderOfpProcess[CounterOfpProcess] = (*this).theLossTable ;
       CounterOfpProcess++;
    }
    else
    {
       RecorderOfpbarProcess[CounterOfpbarProcess] = (*this).theLossTable ;
       CounterOfpbarProcess++;
    }
    if(CutInRange != lastCutInRange)
    {
       lastCutInRange = CutInRange ;
       BuildLambdaTable(aParticleType) ;
    }
    //  G4PAIenergyLoss::BuildDEDXTable(aParticleType) ;
}

////////////////////////////////////////////////////////////////////////////
//
// Build tables for the ionization energy loss
//  the tables are built for MATERIALS
//                           *********

void 
G4PAIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
//G4PAIonisation::BuildPAIonisationTable()
{
   G4double  Charge = aParticleType.GetPDGCharge() ;

   G4double LowEdgeEnergy , ionloss ;
   G4double ParticleMass , RateMass ;
   G4bool isOutRange ;
   static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   const G4double SmallIonLoss = DBL_MIN ;
   const G4double twoln10 = 2.*log(10.) ;
   const G4double Factor = twopi_mc2_rcl2 ;
   const G4double bg2lim = 0.0169 , taulim = 8.4146e-3 ;

   // cuts for p/pbar and electron 

   //   /* *********************************************
   if(Charge>0.)
   {
      ParticleCutInKineticEnergy = theProton->GetCutsInEnergy() ;
   }
   else
   {
      ParticleCutInKineticEnergy = theAntiProton->GetCutsInEnergy() ;
   }
   DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
//   **************************************************  */

   ParticleMass = proton_mass_c2;
   RateMass = electron_mass_c2/ParticleMass ;
   G4int numOfMaterials = theMaterialTable->length();  //  create table

   if ( theLossTable)
   {
      theLossTable->clearAndDestroy();
      delete theLossTable;
   }
      theLossTable = new G4PhysicsTable(numOfMaterials);

   //   theLossTable = new G4PhysicsTable(1);
   
   if( fPAItransferBank )    
   {
     fPAItransferBank->clearAndDestroy() ;
     delete fPAItransferBank ;
   }
   fPAItransferBank = new G4PhysicsTable(TotBin) ;


   for (G4int J=0; J<numOfMaterials; J++)   //  loop for materials
   {
     if( J != fMatIndex ) continue ;   // skip another material 
     
       //create physics vector then fill it ....

     G4PhysicsLogVector* aVector = new G4PhysicsLogVector( LowestKineticEnergy, 
							   HighestKineticEnergy,
							   TotBin               ) ;

      // get material parameters needed for the energy loss calculation

      G4double ElectronDensity, Eexc, Eexc2, Cden, Mden, Aden, X0den, X1den, taul ;
      G4double* ShellCorrectionVector ;
   
      const G4Material* material= (*theMaterialTable)[J];

      ElectronDensity = material->GetElectronDensity();
      Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
      Eexc2 = Eexc*Eexc ;
      Cden = material->GetIonisation()->GetCdensity();
      Mden = material->GetIonisation()->GetMdensity();
      Aden = material->GetIonisation()->GetAdensity();
      X0den = material->GetIonisation()->GetX0density();
      X1den = material->GetIonisation()->GetX1density();
      taul = material->GetIonisation()->GetTaul() ;
      ShellCorrectionVector = material->GetIonisation()
                                      ->GetShellCorrectionVector();

      // get elements in the actual material,
      // they are needed for the low energy part ....

      const G4ElementVector* theElementVector = material->GetElementVector() ;
      const G4double* theAtomicNumDensityVector =
                      material->GetAtomicNumDensityVector() ;
      const G4int NumberOfElements = material->GetNumberOfElements() ;
 
      // get electron cut in kin. energy for the material

      //      DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

      // From gas detector experience

      DeltaCutInKineticEnergyNow = 100*keV ;

      // some local variables 

      G4double tau,tau0,Tmax,gamma,bg2,beta2,rcut,delta,x,sh ;

      // G4cout<<"Material no. = "<<J<<"\t"<<"TotBin = "<<TotBin<<G4endl ;

      for (G4int i = 0 ; i < TotBin ; i++)  //The loop for the kinetic energy 
      {
  	 G4PhysicsFreeVector* transferVector ; 
         LowEdgeEnergy = fProtonEnergyVector->GetLowEdgeEnergy(i) ;
         tau = LowEdgeEnergy/ParticleMass ;



          // high energy part , <dE/dx> according PAI cross section
         {
	    if(tau < 0.01)  // was 0.11, 0.05
	    {
	       tau = 0.01 ;
	    }
	    gamma = tau +1. ;
	    G4cout<<"gamma = "<<gamma<<G4endl ;
            bg2 = tau*(tau+2.) ;
            beta2 = bg2/(gamma*gamma) ;
            Tmax = 2.*electron_mass_c2*bg2
                   /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

	    if ( DeltaCutInKineticEnergyNow > Tmax)         // was <
	    {
               DeltaCutInKineticEnergyNow = Tmax ;
	    }
            G4PAIxSection protonPAI(J,DeltaCutInKineticEnergyNow,bg2,
                                   fSandiaPhotoAbsCof,fSandiaIntervalNumber) ;
	    
	    ionloss = protonPAI.GetMeanEnergyLoss() ;   //  total <dE/dx>

	    G4cout<<"ionloss = "<<ionloss*cm/keV<<" keV/cm"<<G4endl ;
  G4cout<<"n1 = "<<protonPAI.GetIntegralPAIxSection(1)*cm<<" 1/cm"<<G4endl ;
	    // G4cout<<"protonPAI.GetSplineSize() = "<<
            // protonPAI.GetSplineSize()<<G4endl ;

            transferVector = new 
                             G4PhysicsFreeVector(protonPAI.GetSplineSize()) ;

            for(G4int k=0;k<protonPAI.GetSplineSize();k++)
	    {
              transferVector->PutValue( k ,
                                        protonPAI.GetSplineEnergy(k+1),
                                        protonPAI.GetIntegralPAIxSection(k+1) ) ;
	    }
         }
         if ( ionloss <= 0.)
	 {
            ionloss = SmallIonLoss ;
	 }
         aVector->PutValue(i,ionloss) ;

         fPAItransferBank->insertAt(i,transferVector) ;

            // delete[] transferVector ;
      }                                        // end of Tkin loop
      theLossTable->insert(aVector);
   }                                           // end of material loop
   // G4cout<<"G4PAIonisation::BuildPAIonisationTable() have been called"<<G4endl ;
   // G4cout<<"G4PAIonisation::BuildLossTable() have been called"<<G4endl ;
}

///////////////////////////////////////////////////////////////////////
//
// Build mean free path tables for the delta ray production process
//     tables are built for MATERIALS 
//

void 
G4PAIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
    G4double LowEdgeEnergy , Value ,sigma ;
    G4bool isOutRange ;
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    const G4double BigValue = DBL_MAX ;

    G4int numOfMaterials = theMaterialTable->length();    //create table

    if (theMeanFreePathTable) 
    {
       theMeanFreePathTable->clearAndDestroy();
       delete theMeanFreePathTable;
    }
    theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

    // get electron and particle cuts in kinetic energy

    DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
    ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

    for (G4int J=0 ; J < numOfMaterials; J++)  // loop for materials 
    { 
       //create physics vector then fill it ....

       G4PhysicsLogVector* aVector = new G4PhysicsLogVector( LowestKineticEnergy, 
							     HighestKineticEnergy,
							     TotBin             ) ;

       // compute the (macroscopic) cross section first
 
       const G4Material* material= (*theMaterialTable)[J] ;
        
       const G4ElementVector* theElementVector= material->GetElementVector() ;
       const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector();
       const G4int NumberOfElements = material->GetNumberOfElements() ;
 
       // get the electron kinetic energy cut for the actual material,
       //  it will be used in ComputeMicroscopicCrossSection
       // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )

       DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

       for ( G4int i = 0 ; i < TotBin ; i++ )
       {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
          sigma = 0. ;
          for (G4int iel=0; iel<NumberOfElements; iel++ )
          {
              sigma +=  theAtomicNumDensityVector[iel]*
                        ComputeMicroscopicCrossSection(aParticleType,
                        LowEdgeEnergy,
                       (*theElementVector)(iel)->GetZ() ) ;
          }
          // mean free path = 1./macroscopic cross section

           Value = sigma <= 0 ? BigValue: 1./sigma ;     

           aVector->PutValue(i, Value) ;
        }
        theMeanFreePathTable->insert(aVector);
    }
}

////////////////////////////////////////////////////////////////////////////
//
// Cross section formula is OK for spin=0 and 1/2 only !
// Calculates the microscopic cross section in GEANT4 internal units
// ( it is called for elements , AtomicNumber = Z )
//

G4double 
G4PAIonisation::
ComputeMicroscopicCrossSection( const G4ParticleDefinition& aParticleType,
                                      G4double KineticEnergy ,
                                      G4double AtomicNumber               )
{
    G4double TotalEnergy, ParticleMass, betasquare,
             MaxKineticEnergyTransfer,
	     TotalCrossSection, tempvar ;
    const G4double SmallCrossSection = DBL_MIN;

    ParticleMass=aParticleType.GetPDGMass() ; // get particle data 
    TotalEnergy=KineticEnergy + ParticleMass;

    betasquare = KineticEnergy*(TotalEnergy+ParticleMass)   //  kinematics
                 /(TotalEnergy*TotalEnergy);

    tempvar = ParticleMass+electron_mass_c2;

    MaxKineticEnergyTransfer = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+ParticleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

    //  total cross section

    if( MaxKineticEnergyTransfer > DeltaCutInKineticEnergyNow )
    {
       tempvar=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer;
       TotalCrossSection = (1.-tempvar*(1.-betasquare*log(tempvar)))
                           /DeltaCutInKineticEnergyNow;
                           // +term for spin=1/2 particle
       if(aParticleType.GetPDGSpin() == 1)
       {
          TotalCrossSection +=  0.5
                            *(MaxKineticEnergyTransfer-DeltaCutInKineticEnergyNow)
                            /(TotalEnergy*TotalEnergy);
       }
       TotalCrossSection = twopi_mc2_rcl2 * AtomicNumber
                           *TotalCrossSection/betasquare;
    }
    else
    {
       TotalCrossSection=SmallCrossSection ;
    }
    return TotalCrossSection ;
}
 
///////////////////////////////////////////////////////////////////////////
//
// Units are expressed in GEANT4 internal units.
//
 
G4VParticleChange* 
G4PAIonisation::PostStepDoIt( const G4Track& trackData,   
                              const G4Step& stepData          )         
{
   const G4DynamicParticle* aParticle ;
   G4Material* aMaterial;
   G4double KineticEnergy, TotalEnergy, ParticleMass, TotalMomentum,
           betasquare, MaxKineticEnergyTransfer, DeltaKineticEnergy, 
	   DeltaTotalMomentum, costheta, sintheta, phi, dirx, diry, 
	   dirz, finalKineticEnergy, finalPx, finalPy, finalPz, x, xc, 
	   te2, grej, Psquare, Esquare, summass, rate,grejc,finalMomentum ;
  
    G4double Charge;
    
    aParticleChange.Initialize(trackData) ;
    aMaterial = trackData.GetMaterial() ;
    aParticle = trackData.GetDynamicParticle() ;

    Charge=aParticle->GetDefinition()->GetPDGCharge();
    KineticEnergy=aParticle->GetKineticEnergy();
    ParticleMass=aParticle->GetDefinition()->GetPDGMass();
    TotalEnergy=KineticEnergy + ParticleMass ;
    Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
    Esquare=TotalEnergy*TotalEnergy ;
    summass = ParticleMass + electron_mass_c2 ;    
    G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;

    //  get kinetic energy cut for the electron....

    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];
 
    betasquare=Psquare/Esquare ;  //  kinematics

    MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);

    // sampling kinetic energy of the delta ray 

    if( MaxKineticEnergyTransfer <= DeltaCutInKineticEnergyNow ) // no change at all
    {
       //return &aParticleChange;
       return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
    }
    else // normal case 
    {
       xc=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer ;
       rate=MaxKineticEnergyTransfer/TotalEnergy ;

       if(aParticle->GetDefinition()->GetPDGSpin() == 1) 
       {
          te2 = 0.5*rate*rate ;
       }
       else
       {
          te2 = 0.0 ;
       }
       grejc=1.-betasquare*xc+te2*xc*xc ;   // sampling follows ...

       do 
       {
          x=xc/(1.-(1.-xc)*G4UniformRand());
          grej=(1.-x*(betasquare-x*te2))/grejc ;
       }
       while( G4UniformRand()>grej );
    }
    DeltaKineticEnergy = x * MaxKineticEnergyTransfer ;
    if(DeltaKineticEnergy <= 0.)
    {
       return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
    }
    DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                               2. * electron_mass_c2 )) ;
    TotalMomentum = sqrt(Psquare) ;
    costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
            /(DeltaTotalMomentum * TotalMomentum) ;

   

    if ( costheta < -1. )  //  protection against costheta > 1 or < -1 
    {
        costheta = -1. ;
    }
    if ( costheta > +1. ) 
    {
        costheta = +1. ; 
    }                    //  direction of the delta electron  ........

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
    theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z()); 
    theDeltaRay->SetDefinition(G4Electron::Electron());

    // fill aParticleChange
 
    finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ; 

    if (finalKineticEnergy > 0.)
    {
        // changed energy and momentum of the actual particle
       finalMomentum=sqrt(finalKineticEnergy*
                         (finalKineticEnergy+2.*ParticleMass)) ;

       finalPx = (TotalMomentum*ParticleDirection.x()
                -DeltaTotalMomentum*DeltaDirection.x())/finalMomentum ; 
       finalPy = (TotalMomentum*ParticleDirection.y()
                -DeltaTotalMomentum*DeltaDirection.y())/finalMomentum ; 
       finalPz = (TotalMomentum*ParticleDirection.z()
                -DeltaTotalMomentum*DeltaDirection.z())/finalMomentum ; 

       aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    }
    else
    {
       finalKineticEnergy = 0. ;
       if (aParticle->GetDefinition()->GetParticleName() == "proton")
       {
                aParticleChange.SetStatusChange(fStopAndKill);
       }
       else     aParticleChange.SetStatusChange(fStopButAlive);
    }
    aParticleChange.SetEnergyChange( finalKineticEnergy );
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary( theDeltaRay );
    aParticleChange.SetLocalEnergyDeposit (0.);
      
   // ResetNumberOfInteractionLengthLeft;

    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}


/////////////////////////////////////////////////////////////////////////
//
// compute the energy loss after a Step 
//

G4VParticleChange* G4PAIonisation::AlongStepDoIt( const G4Track& trackData,
                                                   const G4Step& stepData    ) 
{
  //  G4cout<<"G4PAIonisation::AlongStepDoIt is called"<<G4endl ;
  
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4bool isOut;
  G4double E,ScaledE,finalT,Step,Tbin,rangebin ;
  const G4double smallLoss=DBL_MIN;
  const G4double BigRange = DBL_MAX ;
  G4int index ;
  G4double cc,discr ; 

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;
  index = aMaterial->GetIndex() ;

  // get the actual (true) Step length from stepData 
  // there is no loss for Step=0. !

  Step = stepData.GetStepLength() ;

  if( Step == 0. || index != fMatIndex ) return &aParticleChange ; 
 
  G4cout<<"step = "<<Step/mm<<" mm"<<G4endl ;




  // get particle and material pointers from trackData
 
  aParticle = trackData.GetDynamicParticle() ;


  E = aParticle->GetKineticEnergy() ;

  G4double Charge = aParticle->GetDefinition()->GetPDGCharge() ;

  G4double Chargesquare = Charge*Charge ;

  G4double MassRatio = proton_mass_c2/aParticle->GetDefinition()->GetPDGMass() ;

  ScaledE = E*MassRatio ;


  ParticleCutInKineticEnergyNow =
               (aParticle->GetDefinition()->GetEnergyCuts())[index] ;

  if(Step >= BigRange)
  {
    finalT = E ;
    fMeanLoss = 0. ;
  }
  else  // here comes the 'real' energy loss calculation (material is NOT vacuum)
  {
	  //  fMeanLoss = ScaledE-0.5*(discr-RangeCoeffB)/RangeCoeffA ;

          //  now the loss with fluctuation

          finalT = E-GetLossWithFluct(Step,aParticle,aMaterial)*Chargesquare ;

          if (finalT<0.) finalT = 0. ;

          fMeanLoss *= Chargesquare ;
    
  }
  //  kill the particle if the kinetic energy <= 0  

  if (finalT <= 0. )
    {
      finalT = 0.;
      if (aParticle->GetDefinition()->GetParticleName() == "proton")
               aParticleChange.SetStatusChange(fStopAndKill);
      else     aParticleChange.SetStatusChange(fStopButAlive); 
    } 

  aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;

}

///////////////////////////////////////////////////////////////////////
//
//

G4double  
G4PAIonisation::GetLossWithFluct( G4double Step,
                                   const G4DynamicParticle* aParticle,
                                   G4Material* aMaterial               )
{  
  G4int iTkin, iTransfer  ;
  G4long iCollision, numOfCollisions ;
  G4int       index = aMaterial->GetIndex() ;
  G4bool isOutRange ;

  // G4cout<<"G4PAIenergyLoss::GetLossWithFluct"<<G4endl ;

  G4double loss = 0.0 ;
  G4double transfer, position, E1, E2, W1, W2, W, firstMu, secondMu ;
  G4double      Tkin = aParticle->GetKineticEnergy() ;
  G4double MassRatio = proton_mass_c2/aParticle->GetDefinition()->GetPDGMass() ;
  G4double TkinScaled = Tkin*MassRatio ;
  G4PhysicsLogVector* 
  aLogVector = new G4PhysicsLogVector( G4PAIonisation::GetMinKineticEnergy(),
                                       G4PAIonisation::GetMaxKineticEnergy(),
                                       G4PAIonisation::GetBinNumber()        ) ;

  for(iTkin=0;iTkin<G4PAIonisation::GetBinNumber();iTkin++)
  {
    if(TkinScaled < aLogVector->GetLowEdgeEnergy(iTkin)) // <= ?
    {
      break ;
    } 
  }
  G4int iPlace = iTkin - 1 ; // index*(G4PAIonisation::GetBinNumber()) +

  G4cout<<"iPlace = "<<iPlace<<G4endl ;

  G4PhysicsVector*  firstVector = (*fPAItransferBank)(iPlace)     ;
  G4PhysicsVector* secondVector = (*fPAItransferBank)(iPlace + 1) ;

  if(iTkin == G4PAIonisation::GetBinNumber()) // Fermi plato, try from left
  {
    numOfCollisions = RandPoisson::shoot((*(*fPAItransferBank)(iPlace))(0)*Step) ;

      G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

    while(numOfCollisions)
    {
      position = (*(*fPAItransferBank)(iPlace))(0)*G4UniformRand() ;

      for(iTransfer=0;;iTransfer++)
      {
        if(position >= (*(*fPAItransferBank)(iPlace))(iTransfer)) break ;
      }
      loss += (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
      numOfCollisions-- ;
    }
  }
  else
  {
    if(iTkin == 0) // Tkin is too small, trying from right only
    {
      numOfCollisions = RandPoisson::
                        shoot((*(*fPAItransferBank)(iPlace+1))(0)*Step) ;

      G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

      while(numOfCollisions)
      {
        position = (*(*fPAItransferBank)(iPlace+1))(0)*G4UniformRand() ;

        for(iTransfer=0;;iTransfer++)
        {
          if(position >= (*(*fPAItransferBank)(iPlace+1))(iTransfer)) break ;
        }
        loss += (*fPAItransferBank)(iPlace+1)->GetLowEdgeEnergy(iTransfer) ;
        numOfCollisions-- ;
      }
    } 
    else // general case: Tkin between two vectors of the material
    {
      E1 = aLogVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = aLogVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;

      // G4cout<<"(*(*fPAItransferBank)(iPlace))(0) = "<<
      //   (*(*fPAItransferBank)(iPlace))(0)<<G4endl ;
      // G4cout<<"(*(*fPAItransferBank)(iPlace+1))(0) = "<<
      //     (*(*fPAItransferBank)(iPlace+1))(0)<<G4endl ;

      numOfCollisions = RandPoisson::shoot(
                     ( (*(*fPAItransferBank)(iPlace))(0)*W1 + 
                     (*(*fPAItransferBank)(iPlace+1))(0)*W2 )*Step) ;

      G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

      while(numOfCollisions)
      {
        position =( (*(*fPAItransferBank)(iPlace))(0)*W1 + 
                    (*(*fPAItransferBank)(iPlace+1))(0)*W2 )*G4UniformRand() ;

        // G4cout<<position<<"\t" ;

        for(iTransfer=0;;iTransfer++)
        {
          if( position >=
          ( (*(*fPAItransferBank)(iPlace))(iTransfer)*W1 + 
            (*(*fPAItransferBank)(iPlace+1))(iTransfer)*W2) )
          {
	      break ;
	  }
        }
        loss += (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ; 
        numOfCollisions-- ;    
      }
    }
  } 
  G4cout<<"PAI loss = "<<loss/keV<<" keV"<<G4endl ; 
  return loss ;
}


//
//
/////////////////////////////////////////////////////////////////////////
