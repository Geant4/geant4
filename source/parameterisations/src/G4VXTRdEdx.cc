// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VXTRdEdx.cc,v 1.3 2001-02-27 15:24:17 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "G4VXTRdEdx.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"
#include "g4std/complex"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"




////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4VXTRdEdx::G4VXTRdEdx(G4Envelope *anEnvelope, G4double a, G4double b) :
  G4VXrayTRmodel(anEnvelope,a,b)
{
  ;
}

///////////////////////////////////////////////////////////////////////////

G4VXTRdEdx::~G4VXTRdEdx()
{
  ;
}


//////////////////////////////////////////////////////////////////////////////
//
// The main function which is responsible for the treatment of a particle passage
// trough G4Envelope

void G4VXTRdEdx::DoIt( const G4FastTrack& fastTrack , 
		                          G4FastStep&  fastStep         )
{
  G4int iTkin, iPlace,  numOfTR, iTR, iTransfer ;
  G4double energyPos, energyTR, theta, phi, dirX, dirY, dirZ ;
  G4double W, W1, W2, E1, E2 ;

  G4double charge = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGCharge() ;
 
  // Now we are ready to Generate TR photons

  G4double chargeSq = charge*charge ;
  G4double kinEnergy     = fastTrack.GetPrimaryTrack()->GetKineticEnergy() ;
  G4double mass = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGMass() ;
  G4double gamma = 1.0 + kinEnergy/mass ;
  //  G4cout<<"gamma = "<<gamma<<G4endl ;
  G4double massRatio = proton_mass_c2/mass ;
  G4double TkinScaled = kinEnergy*massRatio ;

  G4ParticleMomentum direction(fastTrack.GetPrimaryTrackLocalDirection());

  G4double distance = fastTrack.GetEnvelopeSolid()->
                      DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),
		                    direction) ;

  G4ThreeVector position = fastTrack.GetPrimaryTrackLocalPosition() + 
                           distance*direction ;

  // Set final position:

  fastStep.SetPrimaryTrackFinalPosition(position);


  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break ;    
  }
  iPlace = iTkin - 1 ;

  //  G4ParticleMomentum particleDir = fastTrack.GetPrimaryTrack()->
  //                     GetMomentumDirection() ;

  if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
  {
      return ;
  } 
  else          // general case: Tkin between two vectors of the material
  {
    if(iTkin == fTotBin) 
    {
      numOfTR = RandPoisson::shoot( (*(*fEnergyDistrTable)(iPlace))(0)*chargeSq ) ;
    }
    else
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;
      numOfTR = RandPoisson::shoot( ( (*(*fEnergyDistrTable)(iPlace))(0)*W1+
                                      (*(*fEnergyDistrTable)(iPlace+1))(0)*W2 )
                                      *chargeSq ) ;
    }

  // G4cout<<iTkin<<" mean TR number = "<<(((*(*fEnergyDistrTable)(iPlace))(0)+
  // (*(*fAngleDistrTable)(iPlace))(0))*W1 + 
  //                                ((*(*fEnergyDistrTable)(iPlace + 1))(0)+
  // (*(*fAngleDistrTable)(iPlace + 1))(0))*W2)
  //                                    *chargeSq*0.5<<endl ;

    if( numOfTR == 0 ) // no change, return 
    {
       return ;  
    }
    else
    {
      // G4cout<<"Number of X-ray TR photons = "<<numOfTR<<endl ;

      fastStep.SetNumberOfSecondaries(numOfTR);

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

	  // G4cout<<"energyTR = "<<energyTR/keV<<"keV"<<endl ;

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

        G4DynamicParticle aPhotonTR(G4Gamma::Gamma(),directionTR,energyTR) ;

	// A XTR photon is set along the particle track and is not moved to 
	// the G4Envelope surface as in standard X-ray TR models

        G4ThreeVector positionTR = fastTrack.GetPrimaryTrackLocalPosition() +
	                           G4UniformRand()*distance*direction ;


        fastStep.CreateSecondaryTrack( aPhotonTR, 
                                       positionTR, 
		                       fastTrack.GetPrimaryTrack()->
                                       GetGlobalTime()                 ) ;
      }
      kinEnergy -= sumEnergyTR ;
      fastStep.SetPrimaryTrackFinalKineticEnergy(kinEnergy) ;
    }
  }
  return ;
}


//////////////////////////////////////////////////////////////////////////
//
// Build integral energy distribution of XTR photons

void G4VXTRdEdx::BuildTable() 
{
  G4int iMat, jMat, iTkin, iTR, iPlace ;
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
     G4Integrator<G4VXTRdEdx,G4double(G4VXTRdEdx::*)(G4double)> integral ;
     energyVector->PutValue(fBinTR-1,energySum) ;
     angleVector->PutValue(fBinTR-1,angleSum)   ;

     for( iTR = fBinTR - 2 ; iTR >= 0 ; iTR-- )
     {
        energySum += radiatorCof*fCofTR*integral.Legendre10(
		     this,&G4VXTRdEdx::SpectralXTRdEdx,
                     energyVector->GetLowEdgeEnergy(iTR),
                     energyVector->GetLowEdgeEnergy(iTR+1) ) ; 

	//    angleSum  += fCofTR*integral.Legendre96(
	//       this,&G4VXTRdEdx::AngleXTRdEdx,
	//       angleVector->GetLowEdgeEnergy(iTR),
	//       angleVector->GetLowEdgeEnergy(iTR+1) ) ;

        energyVector->PutValue(iTR,energySum) ;
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
  G4cout<<"total time for build X-ray TR dE/dx tables = "
        <<timer.GetUserElapsed()<<" s"<<G4endl ;
  return ;
}

//////////////////////////////////////////////////////////////////////////
//
//

void G4VXTRdEdx::BuildEnergyTable()
{
  return ;
}

////////////////////////////////////////////////////////////////////////
//
//

void G4VXTRdEdx::BuildAngleTable()
{
  return ;
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

G4complex G4VXTRdEdx::OneInterfaceXTRdEdx( G4double energy,
                                           G4double gamma,
                                           G4double varAngle ) 
{
  G4complex Z1    = GetPlateComplexFZ(energy,gamma,varAngle) ;
  G4complex Z2    = GetGasComplexFZ(energy,gamma,varAngle) ;

  G4complex zOut  = (Z1 - Z2)*(Z1 - Z2) ;
            zOut *= varAngle*energy/hbarc/hbarc ;  
  return    zOut  ;

}


//////////////////////////////////////////////////////////////////////////////
//
// For photon energy distribution tables. Integrate first over angle
//

G4double G4VXTRdEdx::SpectralAngleXTRdEdx(G4double varAngle)
{
  G4double result =  GetStackFactor(fEnergy,fGamma,varAngle) ;
  if(result < 0.0) result = 0.0 ;
  return result ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double G4VXTRdEdx::SpectralXTRdEdx(G4double energy)
{
  fEnergy = energy ;
  G4Integrator<G4VXTRdEdx,G4double(G4VXTRdEdx::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4VXTRdEdx::SpectralAngleXTRdEdx,
                             0.0,0.2*fMaxThetaTR) +
         integral.Legendre10(this,&G4VXTRdEdx::SpectralAngleXTRdEdx,
	                     0.2*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double G4VXTRdEdx::AngleSpectralXTRdEdx(G4double energy)
{
  G4double result =  GetStackFactor(energy,fGamma,fVarAngle) ;
  if(result < 0) result = 0.0 ;
  return result ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4VXTRdEdx::AngleXTRdEdx(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator<G4VXTRdEdx,G4double(G4VXTRdEdx::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4VXTRdEdx::AngleSpectralXTRdEdx,
			     fMinEnergyTR,fMaxEnergyTR) ;
}


//
//
///////////////////////////////////////////////////////////////////////




