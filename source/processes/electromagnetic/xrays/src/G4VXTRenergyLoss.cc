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
// $Id: G4VXTRenergyLoss.cc,v 1.3 2002-01-18 17:26:21 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "G4VXTRenergyLoss.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"
#include "g4std/complex"

#include "G4LogicalVolume.hh"

#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"

// Initialization of local constants

G4double G4VXTRenergyLoss::fTheMinEnergyTR   =    1.0*keV  ;
G4double G4VXTRenergyLoss::fTheMaxEnergyTR   =  100.0*keV  ;
G4double G4VXTRenergyLoss::fTheMaxAngle    =      1.0e-3   ;
G4double G4VXTRenergyLoss::fTheMinAngle    =      5.0e-6   ;
G4int    G4VXTRenergyLoss::fBinTR            =   50        ;

G4double G4VXTRenergyLoss::fMinProtonTkin = 100.0*GeV  ;
G4double G4VXTRenergyLoss::fMaxProtonTkin = 100.0*TeV  ;
G4int    G4VXTRenergyLoss::fTotBin        =  50        ;
// Proton energy vector initialization

G4PhysicsLogVector* G4VXTRenergyLoss::
fProtonEnergyVector = new G4PhysicsLogVector(fMinProtonTkin,
                                             fMaxProtonTkin,
                                                    fTotBin  ) ;

G4double G4VXTRenergyLoss::fPlasmaCof = 4.0*pi*fine_structure_const*
                                       hbarc*hbarc*hbarc/electron_mass_c2 ;

G4double G4VXTRenergyLoss::fCofTR     = fine_structure_const/pi ;





////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4VXTRenergyLoss::G4VXTRenergyLoss(G4LogicalVolume *anEnvelope, 
				   G4Material* foilMat,G4Material* gasMat,
                                    G4double a, G4double b,
                                    G4int n,const G4String& processName) :
  G4VContinuousProcess(processName)
{
  fEnvelope = anEnvelope ;
  //  fPlateNumber = fEnvelope->GetNoDaughters() ;
  fPlateNumber = n ;
  G4cout<<"the number of TR radiator plates = "<<fPlateNumber<<G4endl ;
  if(fPlateNumber == 0)
  {
    G4Exception("No plates in X-ray TR radiator") ;
  }
  // Mean thicknesses of plates and gas gaps

  fPlateThick = a ;
  fGasThick   = b ;

  fTotalDist  = fPlateNumber*(fPlateThick+fGasThick) ;  
  G4cout<<"total radiator thickness = "<<fTotalDist/cm<<" cm"<<G4endl ;

  // index of plate material
  fMatIndex1 = foilMat->GetIndex()  ;
  G4cout<<"plate material = "<<foilMat->GetName()<<G4endl ;

  // index of gas material
  fMatIndex2 = gasMat->GetIndex()  ;
  G4cout<<"gas material = "<<gasMat->GetName()<<G4endl ;

  // plasma energy squared for plate material

  fSigma1 = fPlasmaCof*foilMat->GetElectronDensity()  ;
  //  fSigma1 = (20.9*eV)*(20.9*eV) ;
  G4cout<<"plate plasma energy = "<<sqrt(fSigma1)/eV<<" eV"<<G4endl ;

  // plasma energy squared for gas material

  fSigma2 = fPlasmaCof*gasMat->GetElectronDensity()  ;
  G4cout<<"gas plasma energy = "<<sqrt(fSigma2)/eV<<" eV"<<G4endl ;

  // Compute cofs for preparation of linear photo absorption

  ComputePlatePhotoAbsCof() ;
  ComputeGasPhotoAbsCof() ;

}

///////////////////////////////////////////////////////////////////////////

G4VXTRenergyLoss::~G4VXTRenergyLoss()
{
   G4int i ;
   for(i=0;i<fGasIntervalNumber;i++)
   {
     delete[] fGasPhotoAbsCof[i] ;
   }
   delete[] fGasPhotoAbsCof ;

   for(i=0;i<fPlateIntervalNumber;i++)
   {
     delete[] fPlatePhotoAbsCof[i] ;
   }
   delete[] fPlatePhotoAbsCof ; 
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns condition for application of the model depending on particle type


G4bool G4VXTRenergyLoss::IsApplicable(const G4ParticleDefinition& particle)
{
  return  ( particle.GetPDGCharge() != 0.0 ) ; 
}

//////////////////////////////////////////////////////////////////////////////////
//
// GetContinuousStepLimit
//

G4double 
G4VXTRenergyLoss::GetContinuousStepLimit(const G4Track& aTrack,
				         G4double  ,
				         G4double  ,
                                         G4double& )
{
	G4double StepLimit = DBL_MAX;

	return StepLimit;
}

//////////////////////////////////////////////////////////////////////////
//
// Build integral energy distribution of XTR photons

void G4VXTRenergyLoss::BuildTable() 
{
  G4int iTkin, iTR, iPlace ;
  G4double radiatorCof = 1.0 ;           // for tuning of XTR yield

  fEnergyDistrTable = new G4PhysicsTable(fTotBin) ;

  fGammaTkinCut = 0.0 ;
  
  // setting of min/max TR energies 
  
  if(fGammaTkinCut > fTheMinEnergyTR)  fMinEnergyTR = fGammaTkinCut ;
  else                                 fMinEnergyTR = fTheMinEnergyTR ;
	
  if(fGammaTkinCut > fTheMaxEnergyTR) fMaxEnergyTR = 2.0*fGammaTkinCut ;  
  else                                fMaxEnergyTR = fTheMaxEnergyTR ;

  G4cout.precision(4) ;
  G4Timer timer ;
  timer.Start() ;
	
  for( iTkin = 0 ; iTkin < fTotBin ; iTkin++ )      // Lorentz factor loop
  {
     G4PhysicsLogVector* energyVector = new G4PhysicsLogVector( fMinEnergyTR,
                                                                fMaxEnergyTR,
                                                                fBinTR         ) ;

     fGamma = 1.0 + (fProtonEnergyVector->
                            GetLowEdgeEnergy(iTkin)/proton_mass_c2) ;

     fMaxThetaTR = 25.0/(fGamma*fGamma) ;  // theta^2

     fTheMinAngle = 1.0e-6 ; // was 5.e-6, e-5, e-4
 
     if( fMaxThetaTR > fTheMaxAngle )    fMaxThetaTR = fTheMaxAngle ; 
     else
     {
       if( fMaxThetaTR < fTheMinAngle )  fMaxThetaTR = fTheMinAngle ;
     }

     G4PhysicsLinearVector* angleVector = new G4PhysicsLinearVector(        0.0,
                                                                    fMaxThetaTR,
                                                                    fBinTR      ) ;

     G4double energySum = 0.0 ;
     G4double angleSum  = 0.0 ;
     G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral ;
     energyVector->PutValue(fBinTR-1,energySum) ;
     angleVector->PutValue(fBinTR-1,angleSum)   ;

     for( iTR = fBinTR - 2 ; iTR >= 0 ; iTR-- )
     {
        energySum += radiatorCof*fCofTR*integral.Legendre10(
		     this,&G4VXTRenergyLoss::SpectralXTRdEdx,
                     energyVector->GetLowEdgeEnergy(iTR),
                     energyVector->GetLowEdgeEnergy(iTR+1) ) ; 

	//    angleSum  += fCofTR*integral.Legendre96(
	//       this,&G4VXTRenergyLoss::AngleXTRdEdx,
	//       angleVector->GetLowEdgeEnergy(iTR),
	//       angleVector->GetLowEdgeEnergy(iTR+1) ) ;

        energyVector->PutValue(iTR,energySum/fTotalDist) ;
        //  angleVector ->PutValue(iTR,angleSum)   ;
     }
     G4cout<<iTkin<<"\t"
           <<"fGamma = "<<fGamma<<"\t"  //  <<"  fMaxThetaTR = "<<fMaxThetaTR
           <<"sumE = "<<energySum      // <<" ; sumA = "<<angleSum
           <<G4endl ;
     iPlace = iTkin ;
     fEnergyDistrTable->insertAt(iPlace,energyVector) ;
     //  fAngleDistrTable->insertAt(iPlace,angleVector) ;
  }     
  timer.Stop() ;
  G4cout.precision(6) ;
  G4cout<<G4endl ;
  G4cout<<"total time for build X-ray TR energy loss tables = "
        <<timer.GetUserElapsed()<<" s"<<G4endl ;
  return ;
}

//////////////////////////////////////////////////////////////////////////
//
//

void G4VXTRenergyLoss::BuildEnergyTable()
{
  return ;
}

////////////////////////////////////////////////////////////////////////
//
//

void G4VXTRenergyLoss::BuildAngleTable()
{
  return ;
} 


//////////////////////////////////////////////////////////////////////////////
//
// The main function which is responsible for the treatment of a particle passage
// trough G4Envelope

G4VParticleChange* G4VXTRenergyLoss::AlongStepDoIt( const G4Track& aTrack, 
		                                    const G4Step&  aStep   )
{
  G4int iTkin, iPlace,   numOfTR, iTR ;
  G4double energyTR, meanNumOfTR, theta, phi, dirX, dirY, dirZ, rand ;
  G4double W, W1, W2, E1, E2 ;
    G4cout<<"Start of G4VXTRenergyLoss::AlongStepDoIt "<<G4endl ;

  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) 
  {
    G4cout<<"Go out from G4VXTRenergyLoss::AlongStepDoIt: wrong volume "<<G4endl ;
    return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
  }
  aParticleChange.Initialize(aTrack);
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
	
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  
 
  // Now we are ready to Generate TR photons

  G4double chargeSq  = charge*charge ;
  G4double kinEnergy = aParticle->GetKineticEnergy() ;
  G4double mass      = aParticle->GetDefinition()->GetPDGMass() ;
  G4double gamma     = 1.0 + kinEnergy/mass ;

  if(verboseLevel > 0 )
  {
    G4cout<<"gamma = "<<gamma<<G4endl ;
  }
  G4double massRatio = proton_mass_c2/mass ;
  G4double TkinScaled = kinEnergy*massRatio ;

  G4ThreeVector      startPos  = pPreStepPoint->GetPosition();
  G4double           startTime = pPreStepPoint->GetGlobalTime();
  G4ParticleMomentum direction = aParticle->GetMomentumDirection();
  G4double           distance  = aStep.GetStepLength() ;


  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break ;    
  }
  iPlace = iTkin - 1 ;

  if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
  {
    G4cout<<"Go out from G4VXTRenergyLoss::AlongStepDoIt: iTkin=0 "<<G4endl ;
    return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
  } 
  else          // general case: Tkin between two vectors of the material
  {
    if(iTkin == fTotBin) 
    {
      meanNumOfTR = (*(*fEnergyDistrTable)(iPlace))(0)*chargeSq*distance ;
      numOfTR = RandPoisson::shoot(meanNumOfTR) ;
    }
    else
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;
      meanNumOfTR = ( (*(*fEnergyDistrTable)(iPlace  ))(0)*W1+
                      (*(*fEnergyDistrTable)(iPlace+1))(0)*W2 )*chargeSq*distance ;
      
      if(verboseLevel > 0 )
      {
        G4cout<<iTkin<<" mean TR number = "<<meanNumOfTR
              <<" or mean over energy-angle tables "
              <<(((*(*fEnergyDistrTable)(iPlace))(0)+
                  (*(*fAngleDistrTable)(iPlace))(0))*W1 + 
                 ((*(*fEnergyDistrTable)(iPlace + 1))(0)+
                  (*(*fAngleDistrTable)(iPlace + 1))(0))*W2)*chargeSq*0.5
              <<endl ;
      }
      numOfTR = RandPoisson::shoot( meanNumOfTR ) ;
    }
    if( numOfTR == 0 ) // no change, return 
    {
      aParticleChange.SetNumberOfSecondaries(0);
    G4cout<<"Go out from G4VXTRenergyLoss::AlongStepDoIt: numOfTR=0 "<<G4endl ;
      return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep); 
    }
    else
    {
      G4cout<<"Number of X-ray TR photons = "<<numOfTR<<endl ;

      aParticleChange.SetNumberOfSecondaries(numOfTR);

      G4double sumEnergyTR = 0.0 ;

        for(iTR=0;iTR<numOfTR;iTR++)
        {

      //    energyPos = ((*(*fEnergyDistrTable)(iPlace))(0)*W1+
      //          (*(*fEnergyDistrTable)(iPlace + 1))(0)*W2)*G4UniformRand() ;
      //  for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
      //	{
      //    if(energyPos >= ((*(*fEnergyDistrTable)(iPlace))(iTransfer)*W1+
      //                 (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer)*W2)) break ;
      //	}
      //   energyTR = ((*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer))*W1+
      //     ((*fEnergyDistrTable)(iPlace + 1)->GetLowEdgeEnergy(iTransfer))*W2 ;

      energyTR = GetXTRrandomEnergy(TkinScaled,iTkin) ;

	   G4cout<<"energyTR = "<<energyTR/keV<<"keV"<<endl ;

        sumEnergyTR += energyTR ;

        theta = abs(G4RandGauss::shoot(0.0,pi/gamma)) ;

        if( theta >= 0.1 ) theta = 0.1 ;

	// G4cout<<" : theta = "<<theta<<endl ;

        phi = twopi*G4UniformRand() ;

        dirX = sin(theta)*cos(phi)  ;
        dirY = sin(theta)*sin(phi)  ;
        dirZ = cos(theta)           ;

        G4ThreeVector directionTR(dirX,dirY,dirZ) ;
        directionTR.rotateUz(direction) ;
        directionTR.unit() ;

        G4DynamicParticle* aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(),
                                                              directionTR,energyTR) ;

	// A XTR photon is set along the particle track and is not moved to 
	// the G4Envelope surface as in standard X-ray TR models

	rand = G4UniformRand();
        G4double delta = rand*distance ;
	G4double deltaTime = delta /
                       ((pPreStepPoint->GetVelocity()+
                         pPostStepPoint->GetVelocity())/2.);

        G4double aSecondaryTime = startTime + deltaTime;

        G4ThreeVector positionTR = startPos + delta*direction ;


        G4Track* aSecondaryTrack = new G4Track( aPhotonTR, 
		                                aSecondaryTime,positionTR ) ;
        aSecondaryTrack->SetTouchableHandle(aStep.GetPostStepPoint()
                                                  ->GetTouchableHandle());
        aSecondaryTrack->SetParentID(aTrack.GetTrackID());

	aParticleChange.AddSecondary(aSecondaryTrack);
      }
      kinEnergy -= sumEnergyTR ;
      aParticleChange.SetEnergyChange(kinEnergy) ;
    }
  }
  // return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
  return &aParticleChange;
}

///////////////////////////////////////////////////////////////////////
//
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2, or 2->1)
// varAngle =2* (1 - cos(theta)) or approximately = theta*theta
//

G4complex G4VXTRenergyLoss::OneInterfaceXTRdEdx( G4double energy,
                                           G4double gamma,
                                           G4double varAngle ) 
{
  G4complex Z1    = GetPlateComplexFZ(energy,gamma,varAngle) ;
  G4complex Z2    = GetGasComplexFZ(energy,gamma,varAngle) ;

  G4complex zOut  = (Z1 - Z2)*(Z1 - Z2)
                    * (varAngle*energy/hbarc/hbarc) ;  
  return    zOut  ;

}


//////////////////////////////////////////////////////////////////////////////
//
// For photon energy distribution tables. Integrate first over angle
//

G4double G4VXTRenergyLoss::SpectralAngleXTRdEdx(G4double varAngle)
{
  G4double result =  GetStackFactor(fEnergy,fGamma,varAngle) ;
  if(result < 0.0) result = 0.0 ;
  return result ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double G4VXTRenergyLoss::SpectralXTRdEdx(G4double energy)
{
  fEnergy = energy ;
  G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4VXTRenergyLoss::SpectralAngleXTRdEdx,
                             0.0,0.3*fMaxThetaTR) +
         integral.Legendre96(this,&G4VXTRenergyLoss::SpectralAngleXTRdEdx,
	                     0.3*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double G4VXTRenergyLoss::AngleSpectralXTRdEdx(G4double energy)
{
  G4double result =  GetStackFactor(energy,fGamma,fVarAngle) ;
  if(result < 0) result = 0.0 ;
  return result ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4VXTRenergyLoss::AngleXTRdEdx(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre10(this,&G4VXTRenergyLoss::AngleSpectralXTRdEdx,
			     fMinEnergyTR,fMaxEnergyTR) ;
}















//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for plates. Omega is energy !!!

G4double G4VXTRenergyLoss::GetPlateFormationZone( G4double omega ,
                                                G4double gamma ,
                                                G4double varAngle    ) 
{
  G4double cof, lambda ;
  lambda = 1.0/gamma/gamma + varAngle + fSigma1/omega/omega ;
  cof = 2.0*hbarc/omega/lambda ;
  return cof ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates complex formation zone for plates. Omega is energy !!!

G4complex G4VXTRenergyLoss::GetPlateComplexFZ( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle    ) 
{
  G4double cof, length,delta, real, image ;

  length = 0.5*GetPlateFormationZone(omega,gamma,varAngle) ;
  delta  = length*GetPlateLinearPhotoAbs(omega) ;
  cof    = 1.0/(1.0 + delta*delta) ;

  real   = length*cof ;
  image  = real*delta ;

  G4complex zone(real,image); 
  return zone ;
}

////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// plate material

void G4VXTRenergyLoss::ComputePlatePhotoAbsCof() 
{
   G4int i, j, numberOfElements ;
   static const G4MaterialTable* 
   theMaterialTable = G4Material::GetMaterialTable();

   G4SandiaTable thisMaterialSandiaTable(fMatIndex1) ;
   numberOfElements = (*theMaterialTable)[fMatIndex1]->GetNumberOfElements() ;
   G4int* thisMaterialZ = new G4int[numberOfElements] ;

   for(i=0;i<numberOfElements;i++)
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[fMatIndex1]->
                                      GetElement(i)->GetZ() ;
   }
   fPlateIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
   fPlateIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                           (*theMaterialTable)[fMatIndex1]->GetFractionVector() ,
        		     numberOfElements,fPlateIntervalNumber) ;
   
   fPlatePhotoAbsCof = new G4double*[fPlateIntervalNumber] ;

   for(i=0;i<fPlateIntervalNumber;i++)
   {
     fPlatePhotoAbsCof[i] = new G4double[5] ;
   }
   for(i=0;i<fPlateIntervalNumber;i++)
   {
      fPlatePhotoAbsCof[i][0] = thisMaterialSandiaTable.
                                GetPhotoAbsorpCof(i+1,0) ; 
                              
      for(j=1;j<5;j++)
      {
           fPlatePhotoAbsCof[i][j] = thisMaterialSandiaTable.
	                             GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex1]->GetDensity() ;
      }
   }
   delete[] thisMaterialZ ;
   return ;
}

//////////////////////////////////////////////////////////////////////
//
// Returns the value of linear photo absorption coefficient (in reciprocal 
// length) for plate for given energy of X-ray photon omega

G4double G4VXTRenergyLoss::GetPlateLinearPhotoAbs(G4double omega) 
{
  G4int i ;
  G4double omega2, omega3, omega4 ; 

  omega2 = omega*omega ;
  omega3 = omega2*omega ;
  omega4 = omega2*omega2 ;

  for(i=0;i<fPlateIntervalNumber;i++)
  {
    if( omega < fPlatePhotoAbsCof[i][0] ) break ;
  }
  if( i == 0 )
  { 
    G4Exception("Invalid (<I1) energy in G4VXTRenergyLoss::GetPlateLinearPhotoAbs");
  }
  else i-- ;
  
  return fPlatePhotoAbsCof[i][1]/omega  + fPlatePhotoAbsCof[i][2]/omega2 + 
         fPlatePhotoAbsCof[i][3]/omega3 + fPlatePhotoAbsCof[i][4]/omega4  ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for gas. Omega is energy !!!

G4double G4VXTRenergyLoss::GetGasFormationZone( G4double omega ,
                                              G4double gamma ,
                                              G4double varAngle   ) 
{
  G4double cof, lambda ;
  lambda = 1.0/gamma/gamma + varAngle + fSigma2/omega/omega ;
  cof = 2.0*hbarc/omega/lambda ;
  return cof ;

}


//////////////////////////////////////////////////////////////////////
//
// Calculates complex formation zone for gas gaps. Omega is energy !!!

G4complex G4VXTRenergyLoss::GetGasComplexFZ( G4double omega ,
                                           G4double gamma ,
                                           G4double varAngle    ) 
{
  G4double cof, length,delta, real, image ;

  length = 0.5*GetGasFormationZone(omega,gamma,varAngle) ;
  delta  = length*GetGasLinearPhotoAbs(omega) ;
  cof    = 1.0/(1.0 + delta*delta) ;

  real   = length*cof ;
  image  = real*delta ;

  G4complex zone(real,image); 
  return zone ;
}



////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// gas material

void G4VXTRenergyLoss::ComputeGasPhotoAbsCof() 
{
   G4int i, j, numberOfElements ;
   static const G4MaterialTable* 
   theMaterialTable = G4Material::GetMaterialTable();

   G4SandiaTable thisMaterialSandiaTable(fMatIndex2) ;
   numberOfElements = (*theMaterialTable)[fMatIndex2]->GetNumberOfElements() ;
   G4int* thisMaterialZ = new G4int[numberOfElements] ;

   for(i=0;i<numberOfElements;i++)
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[fMatIndex2]->
                                      GetElement(i)->GetZ() ;
   }
   fGasIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
   fGasIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                           (*theMaterialTable)[fMatIndex2]->GetFractionVector() ,
        		     numberOfElements,fGasIntervalNumber) ;
   
   fGasPhotoAbsCof = new G4double*[fGasIntervalNumber] ;

   for(i=0;i<fGasIntervalNumber;i++)
   {
     fGasPhotoAbsCof[i] = new G4double[5] ;
   }
   for(i=0;i<fGasIntervalNumber;i++)
   {
      fGasPhotoAbsCof[i][0] = thisMaterialSandiaTable.
                                GetPhotoAbsorpCof(i+1,0) ; 
                              
      for(j=1;j<5;j++)
      {
           fGasPhotoAbsCof[i][j] = thisMaterialSandiaTable.
	                             GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex2]->GetDensity() ;
      }
   }
   delete[] thisMaterialZ ;
   return ;
}

//////////////////////////////////////////////////////////////////////
//
// Returns the value of linear photo absorption coefficient (in reciprocal 
// length) for gas

G4double G4VXTRenergyLoss::GetGasLinearPhotoAbs(G4double omega) 
{
  G4int i ;
  G4double omega2, omega3, omega4 ; 

  omega2 = omega*omega ;
  omega3 = omega2*omega ;
  omega4 = omega2*omega2 ;

  for(i=0;i<fGasIntervalNumber;i++)
  {
    if( omega < fGasPhotoAbsCof[i][0] ) break ;
  }
  if( i == 0 )
  { 
   G4Exception("Invalid (<I1) energy in G4VXTRenergyLoss::GetGasLinearPhotoAbs");
  }
  else i-- ;
  
  return fGasPhotoAbsCof[i][1]/omega  + fGasPhotoAbsCof[i][2]/omega2 + 
         fGasPhotoAbsCof[i][3]/omega3 + fGasPhotoAbsCof[i][4]/omega4  ;

}

//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// Omega is energy !!!

G4double G4VXTRenergyLoss::GetPlateZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetPlateFormationZone(omega,gamma,varAngle)*GetPlateLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// G4cout and output in file in some energy range.

void G4VXTRenergyLoss::GetPlateZmuProduct() 
{
  G4std::ofstream outPlate("plateZmu.dat", G4std::ios::out ) ;
  outPlate.setf( G4std::ios::scientific, G4std::ios::floatfield );

  G4int i ;
  G4double omega, varAngle, gamma ;
  gamma = 10000. ;
  varAngle = 1/gamma/gamma ;
  G4cout<<"energy, keV"<<"\t"<<"Zmu for plate"<<G4endl ;
  for(i=0;i<100;i++)
  {
    omega = (1.0 + i)*keV ;
    G4cout<<omega/keV<<"\t"<<GetPlateZmuProduct(omega,gamma,varAngle)<<"\t" ;
    outPlate<<omega/keV<<"\t\t"<<GetPlateZmuProduct(omega,gamma,varAngle)<<G4endl ;
  }
  return  ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for gas. 
// Omega is energy !!!

G4double G4VXTRenergyLoss::GetGasZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetGasFormationZone(omega,gamma,varAngle)*GetGasLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof byformation zone for gas. 
// G4cout and output in file in some energy range.

void G4VXTRenergyLoss::GetGasZmuProduct() 
{
  G4std::ofstream outGas("gasZmu.dat", G4std::ios::out ) ;
  outGas.setf( G4std::ios::scientific, G4std::ios::floatfield );
  G4int i ;
  G4double omega, varAngle, gamma ;
  gamma = 10000. ;
  varAngle = 1/gamma/gamma ;
  G4cout<<"energy, keV"<<"\t"<<"Zmu for gas"<<G4endl ;
  for(i=0;i<100;i++)
  {
    omega = (1.0 + i)*keV ;
    G4cout<<omega/keV<<"\t"<<GetGasZmuProduct(omega,gamma,varAngle)<<"\t" ;
    outGas<<omega/keV<<"\t\t"<<GetGasZmuProduct(omega,gamma,varAngle)<<G4endl ;
  }
  return  ;
}

///////////////////////////////////////////////////////////////////////
//
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2, or 2->1)
// varAngle =2* (1 - cos(theta)) or approximately = theta*theta
//

G4double
G4VXTRenergyLoss::OneBoundaryXTRNdensity( G4double energy,G4double gamma,
                                         G4double varAngle ) const
{
  G4double  formationLength1, formationLength2 ;
  formationLength1 = 1.0/
  (1.0/(gamma*gamma)
  + fSigma1/(energy*energy)
  + varAngle) ;
  formationLength2 = 1.0/
  (1.0/(gamma*gamma)
  + fSigma2/(energy*energy)
  + varAngle) ;
  return (varAngle/energy)*(formationLength1 - formationLength2)
              *(formationLength1 - formationLength2)  ;

}


//////////////////////////////////////////////////////////////////////////////
//
// For photon energy distribution tables. Integrate first over angle
//

G4double G4VXTRenergyLoss::XTRNSpectralAngleDensity(G4double varAngle)
{
  return OneBoundaryXTRNdensity(fEnergy,fGamma,varAngle)*
         GetStackFactor(fEnergy,fGamma,varAngle)             ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double G4VXTRenergyLoss::XTRNSpectralDensity(G4double energy)
{
  fEnergy = energy ;
  G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4VXTRenergyLoss::XTRNSpectralAngleDensity,
                             0.0,0.2*fMaxThetaTR) +
         integral.Legendre10(this,&G4VXTRenergyLoss::XTRNSpectralAngleDensity,
	                     0.2*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double G4VXTRenergyLoss::XTRNAngleSpectralDensity(G4double energy)
{
  return OneBoundaryXTRNdensity(energy,fGamma,fVarAngle)*
         GetStackFactor(energy,fGamma,fVarAngle)             ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4VXTRenergyLoss::XTRNAngleDensity(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4VXTRenergyLoss::XTRNAngleSpectralDensity,
			     fMinEnergyTR,fMaxEnergyTR) ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Check number of photons for a range of Lorentz factors from both energy 
// and angular tables

void G4VXTRenergyLoss::GetNumberOfPhotons()
{
  G4int iTkin ;
  G4double gamma, numberE ;

  G4std::ofstream outEn("numberE.dat", G4std::ios::out ) ;
  outEn.setf( G4std::ios::scientific, G4std::ios::floatfield );

  G4std::ofstream outAng("numberAng.dat", G4std::ios::out ) ;
  outAng.setf( G4std::ios::scientific, G4std::ios::floatfield );

  for(iTkin=0;iTkin<fTotBin;iTkin++)      // Lorentz factor loop
  {
     gamma = 1.0 + (fProtonEnergyVector->
                            GetLowEdgeEnergy(iTkin)/proton_mass_c2) ;
     numberE = (*(*fEnergyDistrTable)(iTkin))(0) ;
     //  numberA = (*(*fAngleDistrTable)(iTkin))(0) ;
     G4cout<<gamma<<"\t\t"<<numberE<<"\t"    //  <<numberA
           <<G4endl ; 
     outEn<<gamma<<"\t\t"<<numberE<<G4endl ; 
     //  outAng<<gamma<<"\t\t"<<numberA<<G4endl ; 
  }
  return ;
}  

/////////////////////////////////////////////////////////////////////////
//
// Returns randon energy of a X-ray TR photon for given scaled kinetic energy
// of a charged particle

G4double G4VXTRenergyLoss::GetXTRrandomEnergy( G4double scaledTkin, G4int iTkin )
{
  G4int iTransfer, iPlace  ;
  G4double transfer = 0.0, position, E1, E2, W1, W2, W ;

  iPlace = iTkin - 1 ;

  //  G4cout<<"iPlace = "<<iPlace<<endl ;

  if(iTkin == fTotBin) // relativistic plato, try from left
  {
      position = (*(*fEnergyDistrTable)(iPlace))(0)*G4UniformRand() ;

      for(iTransfer=0;;iTransfer++)
      {
        if(position >= (*(*fEnergyDistrTable)(iPlace))(iTransfer)) break ;
      }
      transfer = GetXTRenergy(iPlace,position,iTransfer);
  }
  else
  {
    E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
    E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
    W  = 1.0/(E2 - E1) ;
    W1 = (E2 - scaledTkin)*W ;
    W2 = (scaledTkin - E1)*W ;

    position =( (*(*fEnergyDistrTable)(iPlace))(0)*W1 + 
                    (*(*fEnergyDistrTable)(iPlace+1))(0)*W2 )*G4UniformRand() ;

        // G4cout<<position<<"\t" ;

    for(iTransfer=0;;iTransfer++)
    {
          if( position >=
          ( (*(*fEnergyDistrTable)(iPlace))(iTransfer)*W1 + 
            (*(*fEnergyDistrTable)(iPlace+1))(iTransfer)*W2) ) break ;
    }
    transfer = GetXTRenergy(iPlace,position,iTransfer);
    
  } 
  //  G4cout<<"XTR transfer = "<<transfer/keV<<" keV"<<endl ; 
  if(transfer < 0.0 ) transfer = 0.0 ;
  return transfer ;
}

////////////////////////////////////////////////////////////////////////
//
// Returns approximate position of X-ray photon energy during random sampling
// over integral energy distribution

G4double G4VXTRenergyLoss::GetXTRenergy( G4int    iPlace, 
                                       G4double position, 
                                       G4int    iTransfer )
{
  G4double x1, x2, y1, y2, result ;

  if(iTransfer == 0)
  {
    result = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  }  
  else
  {
    y1 = (*(*fEnergyDistrTable)(iPlace))(iTransfer-1) ;
    y2 = (*(*fEnergyDistrTable)(iPlace))(iTransfer) ;

    x1 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
    x2 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;

    if ( x1 == x2 )    result = x2 ;
    else
    {
      if ( y1 == y2  ) result = x1 + (x2 - x1)*G4UniformRand() ;
      else
      {
        result = x1 + (position - y1)*(x2 - x1)/(y2 - y1) ;
      }
    }
  }
  return result ;
}



//
//
///////////////////////////////////////////////////////////////////////

