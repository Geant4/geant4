//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MuIonisation.cc,v 1.25 2002-12-04 14:51:54 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------- G4MuIonisation physics process ------------------------------
//                 by Laszlo Urban, September 1997 
// -----------------------------------------------------------------------------
//
// 08-04-98 remove 'tracking cut' of the ionizing particle (mma)
// 26-10-98 new stuff from R.Kokoulin + cleanup , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 23-03-01 R.Kokoulin's correction is commented out, L.Urban
// 29-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma) 
// 28-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 26-09-01 completion of RetrievePhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)  
// 07-11-01 correction(Tmax+xsection computation) L.Urban
// 08-11-01 particleMass becomes a local variable (mma)
// 04-12-02 fix misprint in majorant in PostStep (VI)
// -----------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4MuIonisation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuIonisation::LowerBoundLambda = 1.*keV;
G4double G4MuIonisation::UpperBoundLambda = 1000000.*TeV;
G4int	 G4MuIonisation::NbinLambda = 150;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4MuIonisation::G4MuIonisation(const G4String& processName)
   : G4VMuEnergyLoss(processName),
     theMeanFreePathTable(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     
G4MuIonisation::~G4MuIonisation() 
{
  if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::SetLowerBoundLambda(G4double val)
 {LowerBoundLambda = val;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::SetUpperBoundLambda(G4double val)
 {UpperBoundLambda = val;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::SetNbinLambda(G4int n)
 {NbinLambda = n;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	    
G4double G4MuIonisation::GetLowerBoundLambda()
 { return LowerBoundLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuIonisation::GetUpperBoundLambda()
 { return UpperBoundLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4MuIonisation::GetNbinLambda()
 {return NbinLambda;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::BuildPhysicsTable(const G4ParticleDefinition& ParticleType)
// just call BuildLossTable+BuildLambdaTable
{
  // get bining from EnergyLoss
  LowestKineticEnergy  = GetLowerBoundEloss();
  HighestKineticEnergy = GetUpperBoundEloss();
  TotBin               = GetNbinEloss();

  BuildLossTable(ParticleType);

  if (ParticleType.GetPDGCharge() > 0.)
    {
      RecorderOfmuplusProcess[CounterOfmuplusProcess]   = (*this).theLossTable;
      CounterOfmuplusProcess++;
    }
  else
    {
      RecorderOfmuminusProcess[CounterOfmuminusProcess] = (*this).theLossTable;
      CounterOfmuminusProcess++;
    }
 
  if( !EqualCutVectors(G4Electron::Electron()
                               ->GetLengthCuts(), lastelectronCutInRange))  
     BuildLambdaTable(ParticleType);
 
  G4VMuEnergyLoss::BuildDEDXTable(ParticleType);

  if(&ParticleType == G4MuonPlus::MuonPlus())  PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
 // Build tables of dE/dx due to  the ionization process
 // the tables are built for *MATERIALS*

 // create table
 //
 const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
 G4int numOfMaterials = G4Material::GetNumberOfMaterials();

 if (theLossTable) {theLossTable->clearAndDestroy(); delete theLossTable;}
 theLossTable = new G4PhysicsTable(numOfMaterials);
  
 // get delta cut in energy 
 G4double* DeltaCutInKineticEnergy = (G4Electron::Electron())->GetEnergyCuts();
  
 //  loop for materials
 //
 for (G4int J=0; J<numOfMaterials; J++)
  {
   // create physics vector and fill it
   G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                       LowestKineticEnergy, HighestKineticEnergy, TotBin);
 
   const G4Material* material= (*theMaterialTable)[J];
 
  // get  electron cut in kinetic energy for the material
  G4double DeltaThreshold = DeltaCutInKineticEnergy[J];

  // now comes the loop for the kinetic energy values
  //
  for (G4int i = 0 ; i < TotBin ; i++)
    {
      G4double dEdx = ComputeRestrictedMeandEdx(aParticleType,
	                                        aVector->GetLowEdgeEnergy(i),
	                                        material,
	                                        DeltaThreshold);
      aVector->PutValue(i,dEdx);
    }
   theLossTable->insert(aVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
 // Build mean free path tables for the delta ray production process
 //     tables are built for MATERIALS 

 //create table
 //
 
 const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
 G4int numOfMaterials = G4Material::GetNumberOfMaterials();

 if (theMeanFreePathTable)
   {theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}

 theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

 // get electron cut in kinetic energy
 G4double* DeltaCutInKineticEnergy = (G4Electron::Electron())->GetEnergyCuts();
 
 // loop for materials 

 for (G4int J=0 ; J < numOfMaterials; J++)
    { 
     //create physics vector then fill it ....
     G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowerBoundLambda,UpperBoundLambda,NbinLambda);

     // compute the (macroscopic) cross section first
 
     const G4Material* material= (*theMaterialTable)[J];       
     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* NbOfAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
     G4int NumberOfElements = material->GetNumberOfElements();
 
     // get the electron kinetic energy cut for the actual material,
     //  it will be used in ComputeCrossSectionPerAtom
     // ( --> it will be the same for all the elements in this material)
     G4double DeltaThreshold =DeltaCutInKineticEnergy[J];

     for ( G4int i = 0 ; i < NbinLambda ; i++ )
        {
          G4double LowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
          G4double sigma = 0.;          
          for (G4int iel=0; iel<NumberOfElements; iel++ )
            {
             sigma +=  NbOfAtomsPerVolume[iel]*
                       ComputeCrossSectionPerAtom(aParticleType,
                                                  LowEdgeEnergy,
		               (*theElementVector)[iel]->GetZ(),
		                               DeltaThreshold);
            }

          // mean free path = 1./macroscopic cross section
          G4double Value = sigma > DBL_MIN ? 1./sigma : DBL_MAX;     
          aVector->PutValue(i, Value);
        }
     theMeanFreePathTable->insert(aVector);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuIonisation::ComputeRestrictedMeandEdx (
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
				 const G4Material* material,
				 G4double DeltaThreshold)
{
 // calculate the dE/dx due to the ionization process (Geant4 internal units)
 // Bethe-Bloch formula
 //
 G4double particleMass = aParticleType.GetPDGMass();     
 
 G4double ElectronDensity = material->GetElectronDensity();
 G4double Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
 G4double Eexc2 = Eexc*Eexc;
 
 G4double tau = KineticEnergy/particleMass;
 G4double gamma = tau + 1., bg2 = tau*(tau+2.), beta2 = bg2/(gamma*gamma);
 G4double RateMass = electron_mass_c2/particleMass;
 G4double Tmax=2.*electron_mass_c2*bg2/(1.+2.*gamma*RateMass+RateMass*RateMass);

 G4double taul = material->GetIonisation()->GetTaul();

 G4double dEdx = 0.;
 //
 // high energy part , Bethe-Bloch formula 
 // 
 if (tau > taul)
   {
     G4double rcut = G4std::min(DeltaThreshold/Tmax, 1.);
     dEdx = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)
            +log(rcut)-(1.+rcut)*beta2;
	    
     //density correction 
     G4double Cden   = material->GetIonisation()->GetCdensity();
     G4double Mden   = material->GetIonisation()->GetMdensity();
     G4double Aden   = material->GetIonisation()->GetAdensity();
     G4double X0den  = material->GetIonisation()->GetX0density();
     G4double X1den  = material->GetIonisation()->GetX1density();
  
     const G4double twoln10 = 2.*log(10.); 
     G4double  x = log(bg2)/twoln10;
     G4double delta;
     if (x < X0den) delta = 0.;
     else          {delta = twoln10*x - Cden;
                    if (x < X1den) delta += Aden*pow((X1den-x),Mden);
                   } 

     // shell correction 
     G4double* ShellCorrectionVector = material->GetIonisation()->
                                       GetShellCorrectionVector();
     const G4double bg2lim = 0.0169, taulim = 8.4146e-3;
     G4double sh = 0., xs = 1.;
     if (bg2 > bg2lim) for (G4int k=0; k<3; k++)
                          {xs *= bg2; sh += ShellCorrectionVector[k]/xs;}
     else { for (G4int k=0; k<3; k++)
                       {xs *= bg2lim; sh += ShellCorrectionVector[k]/xs;}
            sh *= log(tau/taul)/log(taulim/taul);     
          }

     // now you can compute the total ionization loss
     dEdx -= (delta + sh); dEdx /= beta2;
     
     // correction of R. Kokoulin  // has been taken out *************** 
     //  G4double E = KineticEnergy+particleMass;
     //  G4double epmax = RateMass*E*E/(RateMass*E+particleMass);
     //  G4double apar = log(2.*epmax/electron_mass_c2);
     //  dEdx += fine_structure_const*(log(2.*E/particleMass)-apar/3.)*
     //                                  apar*apar/twopi; 
     
     dEdx *= twopi_mc2_rcl2*ElectronDensity;
     if (dEdx < 0.) dEdx = 0.;
   }
 //   
 //  low energy part , parametrized energy loss formulae
 //         
 if (tau <= taul)
   {   
     // get elements in the actual material, 
     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* NbOfAtomsPerVolume=material->GetVecNbOfAtomsPerVolume();
     G4int NumberOfElements = material->GetNumberOfElements();
     
     //  loop for the elements in the material
     dEdx = 0.;
     for (G4int iel=0; iel<NumberOfElements; iel++)
        {
          const G4Element* element = (*theElementVector)[iel];         
          if (tau < element->GetIonisation()->GetTau0())  
            dEdx += NbOfAtomsPerVolume[iel]
                       *(element->GetIonisation()->GetAlow()*sqrt(tau)
                       + element->GetIonisation()->GetBlow()*tau);
          else
            dEdx += NbOfAtomsPerVolume[iel]
                       * element->GetIonisation()->GetClow()/sqrt(tau);
        }
     G4double deltaloss = 0.;
     if (DeltaThreshold < Tmax)
       {
         deltaloss = log(Tmax/DeltaThreshold)-
                      beta2*(1.-DeltaThreshold/Tmax) ; 
         if (aParticleType.GetPDGSpin() == 0.5)
            deltaloss += 0.25*(Tmax-DeltaThreshold)*(Tmax-DeltaThreshold)/
                 (KineticEnergy*KineticEnergy+proton_mass_c2*proton_mass_c2);
            deltaloss *= twopi_mc2_rcl2*ElectronDensity/beta2;
       }
     dEdx -= deltaloss;
     if (dEdx < 0.) dEdx = 0.;
   }
 return dEdx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuIonisation::ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber,
				 G4double DeltaThreshold)
{
 // calculates the totalcross section per atom in GEANT4 internal units
 //    ( it is called for elements , AtomicNumber = Z )
 //
     
 G4double particleMass = aParticleType.GetPDGMass();     
 G4double TotalEnergy = KineticEnergy + particleMass;
 G4double tempvar = particleMass+electron_mass_c2;
 G4double KnockonMaxEnergy = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+particleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

 G4double TotalCrossSection = 0.;
 if (KnockonMaxEnergy <= DeltaThreshold) return TotalCrossSection;
    
 const G4double xgi[] = {0.06943,0.33001,0.66999,0.93057};
 const G4double wgi[] = {0.17393,0.32607,0.32607,0.17393};
 const G4double ak1 = 4.6;
 const G4int k2 = 2;

 G4double aaa = log(DeltaThreshold);
 G4double bbb = log(KnockonMaxEnergy);
 G4int    kkk = int((bbb-aaa)/ak1)+k2;
 G4double hhh = (bbb-aaa)/kkk;
 G4double step = exp(hhh);
 G4double ymax = 1./KnockonMaxEnergy;
      
 for (G4int k=0; k<kkk; k++)
    {
      G4double ymin = ymax;
      ymax = ymin*step;
      G4double hhy = ymax-ymin;
      for (G4int i=0; i<4; i++)
         {
           G4double y = ymin+hhy*xgi[i];
           G4double ep = 1./y ;
           TotalCrossSection += ep*ep*wgi[i]*hhy*
                                ComputeDifCrossSectionPerAtom(
                                aParticleType,KineticEnergy,
                                AtomicNumber,ep);
         }
    }
 return TotalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuIonisation::ComputeDifCrossSectionPerAtom(
                                 const G4ParticleDefinition& ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double KnockonEnergy)
 // Calculates the differential cross section per atom
 //   using the cross section formula of R.P. Kokoulin (10/98)
{
  const G4double alphaprime = fine_structure_const/twopi;
  G4double particleMass = ParticleType.GetPDGMass();      
  G4double TotalEnergy = KineticEnergy + particleMass;
  G4double betasquare = KineticEnergy*(TotalEnergy+particleMass)
                      /(TotalEnergy*TotalEnergy);
  G4double tempvar = particleMass+electron_mass_c2;
  G4double KnockonMaxEnergy = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+particleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

  G4double DifCrossSection = 0.;
  if(KnockonEnergy >=  KnockonMaxEnergy)  return DifCrossSection;

  G4double v = KnockonEnergy/TotalEnergy;
  DifCrossSection = twopi_mc2_rcl2*AtomicNumber*
                   (1.-betasquare*KnockonEnergy/KnockonMaxEnergy+0.5*v*v)/
                   (betasquare*KnockonEnergy*KnockonEnergy);
  G4double a1 = log(1.+2.*KnockonEnergy/electron_mass_c2);
  G4double a3 = log(4.*TotalEnergy*(TotalEnergy-KnockonEnergy)/
                    (particleMass*particleMass));
  DifCrossSection *= (1.+alphaprime*a1*(a3-a1)); 

  return DifCrossSection;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4MuIonisation::PostStepDoIt(const G4Track& trackData,   
                                               const G4Step&  stepData)         
{
 aParticleChange.Initialize(trackData);
  
 G4Material* aMaterial = trackData.GetMaterial();
 const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle();

 G4double particleMass = aParticle->GetDefinition()->GetPDGMass();
 G4double KineticEnergy = aParticle->GetKineticEnergy();
 G4double TotalEnergy = KineticEnergy + particleMass;
 G4double Psquare = KineticEnergy*(TotalEnergy+particleMass);
 G4double Esquare = TotalEnergy*TotalEnergy;
 G4double betasquare=Psquare/Esquare; 
 G4double summass = particleMass + electron_mass_c2;
 G4double MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);
 G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();
 
 // get electron cut in kinetic energy
 G4double* DeltaCutInKineticEnergy = (G4Electron::Electron())->GetEnergyCuts();
 G4double DeltaThreshold = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

 // sampling kinetic energy of the delta ray 
 //
  if (MaxKineticEnergyTransfer <= DeltaThreshold)
    // pathological case (it should not happen, there is no change at all)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

 // normal case 
 G4double xc = DeltaThreshold/MaxKineticEnergyTransfer;
 G4double rate = MaxKineticEnergyTransfer/TotalEnergy;
 G4double te2 = 0.5*rate*rate;
 
 // sampling follows ...
 G4double x,twoep,a1,grej;
 const G4double alphaprime = fine_structure_const/twopi; 
 G4double a0=log(2.*TotalEnergy/particleMass); 
 G4double grejc=(1.-xc*(betasquare-xc*te2))*(1.+ alphaprime*a0*a0);
 do { x=xc/(1.-(1.-xc)*G4UniformRand());
      twoep = 2.*x*MaxKineticEnergyTransfer;
      a1    = log(1.+twoep/electron_mass_c2);
      grej  = (1.-x*(betasquare-x*te2))*(1.+alphaprime*a1*
           (a0+log((2.*TotalEnergy-twoep)/particleMass)-a1))/grejc ;
    } while(G4UniformRand() > grej);
    
 G4double  DeltaKineticEnergy = x * MaxKineticEnergyTransfer;

 if (DeltaKineticEnergy <= 0.)
   return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

 G4double DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                               2. * electron_mass_c2 ));
 G4double TotalMomentum = sqrt(Psquare);
 G4double costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
            /(DeltaTotalMomentum * TotalMomentum);

 if (costheta < -1.) costheta = -1.;
 if (costheta > +1.) costheta = +1.;

 //  direction of the delta electron
 //  
 G4double phi = twopi * G4UniformRand(); 
 G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
 G4double dirx = sintheta*cos(phi), diry = sintheta*sin(phi), dirz = costheta;

 G4ThreeVector DeltaDirection(dirx,diry,dirz);
 DeltaDirection.rotateUz(ParticleDirection);

 // create G4DynamicParticle object for delta ray
 //
 G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
 theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
 theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z()); 
 theDeltaRay->SetDefinition(G4Electron::Electron());

 // fill aParticleChange
 // 
 G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy;

 if (finalKineticEnergy > 0.)
   {
    G4double finalPx = TotalMomentum*ParticleDirection.x()
                       - DeltaTotalMomentum*DeltaDirection.x();
    G4double finalPy = TotalMomentum*ParticleDirection.y()
                       - DeltaTotalMomentum*DeltaDirection.y();
    G4double finalPz = TotalMomentum*ParticleDirection.z()
                       - DeltaTotalMomentum*DeltaDirection.z();
    G4double finalMomentum =
              sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
    finalPx /= finalMomentum;
    finalPy /= finalMomentum;
    finalPz /= finalMomentum;
    aParticleChange.SetMomentumChange(finalPx,finalPy,finalPz);
   }
 else
   {
    finalKineticEnergy = 0.;
    aParticleChange.SetStatusChange(fStopButAlive);
   }

 aParticleChange.SetEnergyChange( finalKineticEnergy );
 aParticleChange.SetNumberOfSecondaries(1);   
 aParticleChange.AddSecondary(theDeltaRay);
 aParticleChange.SetLocalEnergyDeposit (0.);
      
 //ResetNumberOfInteractionLengthLeft();
return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MuIonisation::StorePhysicsTable(G4ParticleDefinition* particle,
				              const G4String& directory, 
				              G4bool          ascii)
{
  G4String filename;
  
  // store stopping power table
  filename = GetPhysicsTableFileName(particle,directory,"StoppingPower",ascii);
  if ( !theLossTable->StorePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theLossTable->StorePhysicsTable in " << filename
           << G4endl;
    return false;
  }
  // store mean free path table
  filename = GetPhysicsTableFileName(particle,directory,"MeanFreePath",ascii);
  if ( !theMeanFreePathTable->StorePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theMeanFreePathTable->StorePhysicsTable in " << filename
           << G4endl;
    return false;
  }
  
  G4cout << GetProcessName() << " for " << particle->GetParticleName()
         << ": Success to store the PhysicsTables in "  
         << directory << G4endl;
	 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MuIonisation::RetrievePhysicsTable(G4ParticleDefinition* particle,
					         const G4String& directory, 
				                 G4bool          ascii)
{
  // delete theLossTable and theMeanFreePathTable
  if (theLossTable != 0) {
    theLossTable->clearAndDestroy();
    delete theLossTable; 
  }   
  if (theMeanFreePathTable != 0) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  
  // get bining from EnergyLoss
  LowestKineticEnergy  = GetLowerBoundEloss();
  HighestKineticEnergy = GetUpperBoundEloss();
  TotBin               = GetNbinEloss();
  
  G4String filename;
  
  // retreive stopping power table
  filename = GetPhysicsTableFileName(particle,directory,"StoppingPower",ascii);
  theLossTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());
  if ( !theLossTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theLossTable0->RetrievePhysicsTable in " << filename
           << G4endl;  
    return false;
  }
  
  // retreive mean free path table
  filename = GetPhysicsTableFileName(particle,directory,"MeanFreePath",ascii);
  theMeanFreePathTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());
  if ( !theMeanFreePathTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theMeanFreePathTable->RetrievePhysicsTable in " << filename
           << G4endl;  
    return false;
  }
  
  G4cout << GetProcessName() << " for " << particle->GetParticleName()
         << ": Success to retrieve the PhysicsTables from "
         << directory << G4endl;
	 
  if (particle->GetPDGCharge() > 0.)
    {
      RecorderOfmuplusProcess[CounterOfmuplusProcess]   = (*this).theLossTable;
      CounterOfmuplusProcess++;
    }
  else
    {
      RecorderOfmuminusProcess[CounterOfmuminusProcess] = (*this).theLossTable;
      CounterOfmuminusProcess++;
    }
 

  G4VMuEnergyLoss::BuildDEDXTable(*particle);
  	 
  return true;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . "
            "\n         Good description above the mean excitation energy.\n"
            "         delta ray energy sampled from  differential Xsection.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowerBoundLambda,
                                                  "Energy")
         << " to " << G4BestUnit(UpperBoundLambda,"Energy")
         << " in " << TotBin << " bins. \n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
