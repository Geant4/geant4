// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ForwardXrayTR.cc,v 1.3 1999-12-15 14:52:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4ForwardXrayTR class -- implementation file

// GEANT 4 class implementation file --- Copyright CERN 1995
// CERN Geneva Switzerland

// For information related to this code, please, contact
// CERN, CN Division, ASD Group
// History:
// 1st version 11.09.97 V. Grichine (Vladimir.Grichine@cern.ch )
// 2nd version 17.12.97 V. Grichine


#include <math.h>

// #include "G4ios.hh"
// #include <fstream.h>
// #include <stdlib.h>

#include "G4ForwardXrayTR.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"

// Table initialization

// G4PhysicsTable* G4ForwardXrayTR::fAngleDistrTable  = NULL ;
// G4PhysicsTable* G4ForwardXrayTR::fEnergyDistrTable = NULL ;


// Initialization of local constants

G4int    G4ForwardXrayTR::fSympsonNumber =  100       ;

G4double G4ForwardXrayTR::fTheMinEnergyTR   =    1.0*keV  ;
G4double G4ForwardXrayTR::fTheMaxEnergyTR   =  100.0*keV  ;
G4double G4ForwardXrayTR::fTheMaxAngle    =      1.0e-3   ;
G4double G4ForwardXrayTR::fTheMinAngle    =      5.0e-6   ;
G4int    G4ForwardXrayTR::fBinTR            =   50        ;

G4double G4ForwardXrayTR::fMinProtonTkin = 100.0*GeV  ;
G4double G4ForwardXrayTR::fMaxProtonTkin = 100.0*TeV  ;
G4int    G4ForwardXrayTR::fTotBin        =  50        ;
// Proton energy vector initialization

G4PhysicsLogVector* G4ForwardXrayTR::
fProtonEnergyVector = new G4PhysicsLogVector(fMinProtonTkin,
                                             fMaxProtonTkin,
                                                    fTotBin  ) ;

G4double G4ForwardXrayTR::fPlasmaCof = 4.0*pi*fine_structure_const*
                                       hbarc*hbarc*hbarc/electron_mass_c2 ;

G4double G4ForwardXrayTR::fCofTR     = fine_structure_const/pi ;


/*   ************************************************************************


///////////////////////////////////////////////////////////////////////
//
// Constructor for preparation tables with angle and energy TR distributions
// in all materials involved in test program. Lorentz factors correspond to
// kinetic energies of protons between 100*GeV and 100*TeV, ~ 10^2-10^5
//
// Recommended only for use in applications with 
// few light materials involved                     !!!!!!!!!!!!!!

G4ForwardXrayTR::G4ForwardXrayTR()
  : G4TransitionRadiation("XrayTR")
{
  G4int iMat, jMat, iTkin, iTR, iPlace ;
  static 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
  G4int numOfMat = theMaterialTable->length() ;
  fGammaCutInKineticEnergy = new G4double[numOfMat] ;
  fGammaCutInKineticEnergy = fPtrGamma->GetCutsInEnergy() ;
  fMatIndex1 = -1 ;
  fMatIndex2 = -1 ;
  fAngleDistrTable  = new G4PhysicsTable(numOfMat*(numOfMat - 1)*fTotBin) ;
  fEnergyDistrTable = new G4PhysicsTable(numOfMat*(numOfMat - 1)*fTotBin) ;

        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(fMinProtonTkin,
                                                             fMaxProtonTkin,
                                                                    fTotBin  ) ;

  for(iMat=0;iMat<numOfMat;iMat++) // loop over pairs of different materials
  {
    for(jMat=0;jMat<numOfMat;jMat++)  // transition iMat -> jMat !!!
    {
      if(iMat == jMat)   continue ;      // no TR !!
      else
      {
        const G4Material* mat1 = (*theMaterialTable)[iMat] ;
        const G4Material* mat2 = (*theMaterialTable)[jMat] ;

        fSigma1 = fPlasmaCof*(mat1->GetElectronDensity()) ;
        fSigma2 = fPlasmaCof*(mat2->GetElectronDensity()) ;

//        fGammaTkinCut = fGammaCutInKineticEnergy[jMat] ; // TR photon in jMat !
          fGammaTkinCut = 0.0 ;

        if(fGammaTkinCut > fTheMinEnergyTR) // setting of min/max TR energies 
	{
          fMinEnergyTR = fGammaTkinCut ;
	}
        else
	{
          fMinEnergyTR = fTheMinEnergyTR ;
	}
        if(fGammaTkinCut > fTheMaxEnergyTR)
	{
          fMaxEnergyTR = 2.0*fGammaTkinCut ;    // usually very low TR rate 
	}
        else
	{
          fMaxEnergyTR = fTheMaxEnergyTR ;
	}
        for(iTkin=0;iTkin<fTotBin;iTkin++)      // Lorentz factor loop
	{
          G4PhysicsLogVector* 
                    energyVector = new G4PhysicsLogVector(fMinEnergyTR,
                                                             fMaxEnergyTR,
                                                                   fBinTR  ) ;
          G4PhysicsLinearVector* 
                     angleVector = new G4PhysicsLinearVector(        0.0,
                                                             fMaxThetaTR,
                                                                  fBinTR  ) ;
          G4double energySum = 0.0 ;
          G4double angleSum  = 0.0 ;
          fGamma = 1.0 +   (aVector->GetLowEdgeEnergy(iTkin)/proton_mass_c2) ;
          fMaxThetaTR = 10000.0/(fGamma*fGamma) ;
          if(fMaxThetaTR > fTheMaxAngle)
          { 
            fMaxThetaTR = fTheMaxAngle ;
	  }
          else
	  {
            if(fMaxThetaTR < fTheMinAngle)
	    {
              fMaxThetaTR = fTheMinAngle ;
	    }
	  }
          energyVector->PutValue(fBinTR-1,energySum) ;
          angleVector->PutValue(fBinTR-1,angleSum)   ;

          for(iTR=fBinTR-2;iTR>=0;iTR--)
	  {
            energySum += fCofTR*EnergySum(energyVector->GetLowEdgeEnergy(iTR),
                                        energyVector->GetLowEdgeEnergy(iTR+1)) ; 

            angleSum  += fCofTR*AngleSum(angleVector->GetLowEdgeEnergy(iTR),
                                         angleVector->GetLowEdgeEnergy(iTR+1)) ;
            energyVector->PutValue(iTR,energySum) ;
            angleVector->PutValue(iTR,angleSum)   ;
	  }
          if(jMat < iMat)
	  {
            iPlace = (iMat*(numOfMat-1)+jMat)*fTotBin+iTkin ;
	  }
          else   // jMat > iMat right part of matrices (jMat-1) ! 
	  {
            iPlace = (iMat*(numOfMat-1)+jMat-1)*fTotBin+iTkin ;
	  } 
          fEnergyDistrTable->insertAt(iPlace,energyVector) ;
          fAngleDistrTable->insertAt(iPlace,angleVector) ;
	}    //                      iTkin    
      }      //         jMat != iMat
    }        //     jMat
  }          // iMat
}


****************************************************************  */


//////////////////////////////////////////////////////////////////////
//
// Constructor for creation of physics tables (angle and energy TR 
// distributions) for a couple of selected materials.
//
// Recommended for use in applications with many materials involved,
// when only few (usually couple) materials are interested for generation
// of TR on the interface between them


G4ForwardXrayTR::
G4ForwardXrayTR( const G4String& matName1,   //  G4Material* pMat1,
		 const G4String& matName2,    //  G4Material* pMat2,
                 const G4String& processName                          )
  :        G4TransitionRadiation(processName) 
{
  //  fMatIndex1 = pMat1->GetIndex() ;
  //  fMatIndex2 = pMat2->GetIndex() ;
  G4int iMat, jMat ;
  static 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  G4int numOfMat = theMaterialTable->length() ;

  for(iMat=0;iMat<numOfMat;iMat++)    // check first material name
  {
    if( matName1 == (*theMaterialTable)[iMat]->GetName() )
    {
      fMatIndex1 = (*theMaterialTable)[iMat]->GetIndex() ;
      break ;
    }
  }
  if(iMat == numOfMat)
  {
    G4Exception("Invalid first material name in G4ForwardXrayTR constructor") ;
  }

  for(iMat=0;iMat<numOfMat;iMat++)    // check second material name
  {
    if( matName2 == (*theMaterialTable)[iMat]->GetName() )
    {
      fMatIndex2 = (*theMaterialTable)[iMat]->GetIndex() ;
      break ;
    }
  }
  if(iMat == numOfMat)
  {
    G4Exception("Invalid second material name in G4ForwardXrayTR constructor") ;
  }
  //  G4cout<<"G4ForwardXray constructor is called"<<G4endl ;
  BuildXrayTRtables() ;
}


//////////////////////////////////////////////////////////////////////
//
// Destructor
//

G4ForwardXrayTR::~G4ForwardXrayTR()
{
	;
}

//////////////////////////////////////////////////////////////////////////////
//
// Build physics tables for energy and angular distributions of X-ray TR photon

void G4ForwardXrayTR::BuildXrayTRtables()
{
  G4int iMat, jMat, iTkin, iTR, iPlace ;
  static 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  G4int numOfMat = theMaterialTable->length() ;

  fGammaCutInKineticEnergy = new G4double[numOfMat] ;
  fGammaCutInKineticEnergy = fPtrGamma->GetCutsInEnergy() ;

  fAngleDistrTable  = new G4PhysicsTable(2*fTotBin) ;
  fEnergyDistrTable = new G4PhysicsTable(2*fTotBin) ;


  for(iMat=0;iMat<numOfMat;iMat++)     // loop over pairs of different materials
  {
    if( iMat != fMatIndex1 && iMat != fMatIndex2 ) continue ;

    for(jMat=0;jMat<numOfMat;jMat++)  // transition iMat -> jMat !!!
    {
      if( iMat == jMat || ( jMat != fMatIndex1 && jMat != fMatIndex2 ) )   
      { 
        continue ; 
      } 
      else
      {
        const G4Material* mat1 = (*theMaterialTable)[iMat] ;
        const G4Material* mat2 = (*theMaterialTable)[jMat] ;

        fSigma1 = fPlasmaCof*(mat1->GetElectronDensity()) ;
        fSigma2 = fPlasmaCof*(mat2->GetElectronDensity()) ;

        // fGammaTkinCut = fGammaCutInKineticEnergy[jMat] ; // TR photon in jMat !

        fGammaTkinCut = 0.0 ;

        if(fGammaTkinCut > fTheMinEnergyTR)    // setting of min/max TR energies 
	{
          fMinEnergyTR = fGammaTkinCut ;
	}
        else
	{
          fMinEnergyTR = fTheMinEnergyTR ;
	}
        if(fGammaTkinCut > fTheMaxEnergyTR)
	{
          fMaxEnergyTR = 2.0*fGammaTkinCut ;    // usually very low TR rate 
	}
        else
	{
          fMaxEnergyTR = fTheMaxEnergyTR ;
	}
        for(iTkin=0;iTkin<fTotBin;iTkin++)      // Lorentz factor loop
	{
          G4PhysicsLogVector* 
                    energyVector = new G4PhysicsLogVector( fMinEnergyTR,
                                                           fMaxEnergyTR,
                                                           fBinTR         ) ;

          fGamma = 1.0 +   (fProtonEnergyVector->
                            GetLowEdgeEnergy(iTkin)/proton_mass_c2) ;

          fMaxThetaTR = 10000.0/(fGamma*fGamma) ;
 
          if(fMaxThetaTR > fTheMaxAngle)
          { 
            fMaxThetaTR = fTheMaxAngle ;
	  }
          else
	  {
            if(fMaxThetaTR < fTheMinAngle)
	    {
              fMaxThetaTR = fTheMinAngle ;
	    }
	  }
   // G4cout<<G4endl<<"fGamma = "<<fGamma<<"  fMaxThetaTR = "<<fMaxThetaTR<<G4endl ;
          G4PhysicsLinearVector* 
                     angleVector = new G4PhysicsLinearVector(        0.0,
                                                             fMaxThetaTR,
                                                                  fBinTR  ) ;
          G4double energySum = 0.0 ;
          G4double angleSum  = 0.0 ;

          energyVector->PutValue(fBinTR-1,energySum) ;
          angleVector->PutValue(fBinTR-1,angleSum)   ;

          for(iTR=fBinTR-2;iTR>=0;iTR--)
	  {
            energySum += fCofTR*EnergySum(energyVector->GetLowEdgeEnergy(iTR),
                                        energyVector->GetLowEdgeEnergy(iTR+1)) ; 

            angleSum  += fCofTR*AngleSum(angleVector->GetLowEdgeEnergy(iTR),
                                         angleVector->GetLowEdgeEnergy(iTR+1)) ;

            energyVector->PutValue(iTR,energySum) ;
            angleVector ->PutValue(iTR,angleSum)   ;
	  }
	  // G4cout<<"sumE = "<<energySum<<" ; sumA = "<<angleSum<<G4endl ;

          if(jMat < iMat)
	  {
            iPlace = fTotBin+iTkin ;   // (iMat*(numOfMat-1)+jMat)*
	  }
          else   // jMat > iMat right part of matrices (jMat-1) ! 
	  {
            iPlace = iTkin ;  // (iMat*(numOfMat-1)+jMat-1)*fTotBin+
	  } 
          fEnergyDistrTable->insertAt(iPlace,energyVector) ;
          fAngleDistrTable->insertAt(iPlace,angleVector) ;
	}    //                      iTkin    
      }      //         jMat != iMat
    }        //     jMat
  }          // iMat
  //  G4cout<<"G4ForwardXrayTR::BuildXrayTRtables have been called"<<G4endl ;
}

///////////////////////////////////////////////////////////////////////
//
// This function returns the spectral and angle density of TR quanta
// in X-ray energy region generated forward when a relativistic
// charged particle crosses interface between two materials.
// The high energy small theta approximation is applied.
// (matter1 -> matter2)
// varAngle =2* (1 - cos(Theta)) or approximately = Theta*Theta
//

G4double
G4ForwardXrayTR::SpectralAngleTRdensity( G4double energy,
                                         G4double varAngle ) const
{
  G4double  formationLength1, formationLength2 ;
  formationLength1 = 1.0/
  (1.0/(fGamma*fGamma)
  + fSigma1/(energy*energy)
  + varAngle) ;
  formationLength2 = 1.0/
  (1.0/(fGamma*fGamma)
  + fSigma2/(energy*energy)
  + varAngle) ;
  return (varAngle/energy)*(formationLength1 - formationLength2)
              *(formationLength1 - formationLength2)  ;

}


//////////////////////////////////////////////////////////////////
//
// Analytical formula for angular density of X-ray TR photons
// 

G4double G4ForwardXrayTR::AngleDensity( G4double energy,
                                        G4double varAngle ) const
{
  G4double x, x2, a, b, c, d, f, a2, b2, a4, b4 ;
  G4double cof1, cof2, cof3 ;
  x = 1.0/energy ;
  x2 = x*x ;
  c = 1.0/fSigma1 ;
  d = 1.0/fSigma2 ;
  f = (varAngle + 1.0/(fGamma*fGamma)) ;
  a2 = c*f ;
  b2 = d*f ;
  a4 = a2*a2 ;
  b4 = b2*b2 ;
  a = sqrt(a2) ;
  b = sqrt(b2) ;
  cof1 = c*c*(0.5/(a2*(x2 +a2)) +0.5*log(x2/(x2 +a2))/a4) ;
  cof3 = d*d*(0.5/(b2*(x2 +b2)) +0.5*log(x2/(x2 +b2))/b4) ;
  cof2 = -c*d*(log(x2/(x2 +b2))/b2 - log(x2/(x2 +a2))/a2)/(a2 - b2)   ;
  return -varAngle*(cof1 + cof2 + cof3) ;
}

/////////////////////////////////////////////////////////////////////
//
// Definite integral of X-ray TR spectral-angle density from energy1
// to energy2
//

G4double G4ForwardXrayTR::EnergyInterval( G4double energy1,
                                          G4double energy2,
                                          G4double varAngle ) const
{
  return     AngleDensity(energy2,varAngle)
           - AngleDensity(energy1,varAngle) ;
}

//////////////////////////////////////////////////////////////////////
//
// Integral angle distribution of X-ray TR photons based on analytical
// formula for angle density
//

G4double G4ForwardXrayTR::AngleSum( G4double varAngle1,
                                    G4double varAngle2     )   const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(varAngle2 - varAngle1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
   sumEven += EnergyInterval(fMinEnergyTR,fMaxEnergyTR,varAngle1 + 2*i*h   ) ;
   sumOdd  += EnergyInterval(fMinEnergyTR,fMaxEnergyTR,
                                                   varAngle1 + (2*i - 1)*h ) ;
  }
  sumOdd += EnergyInterval(fMinEnergyTR,fMaxEnergyTR,
                        varAngle1 + (2*fSympsonNumber - 1)*h    ) ;

  return h*(EnergyInterval(fMinEnergyTR,fMaxEnergyTR,varAngle1)
          + EnergyInterval(fMinEnergyTR,fMaxEnergyTR,varAngle2)
          + 4.0*sumOdd + 2.0*sumEven                          )/3.0 ;
}

/////////////////////////////////////////////////////////////////////
//
// Analytical Expression for   spectral density of Xray TR photons
// x = 2*(1 - cos(Theta)) ~ Theta^2
//

G4double G4ForwardXrayTR::SpectralDensity( G4double energy,
                                           G4double      x  ) const
{
  G4double a, b ;
  a =  1.0/(fGamma*fGamma)
     + fSigma1/(energy*energy)  ;
  b =  1.0/(fGamma*fGamma)
     + fSigma2/(energy*energy)  ;
  return ( (a + b)*log((x + b)/(x + a))/(a - b)
          + a/(x + a) + b/(x + b) )/energy ;

}

////////////////////////////////////////////////////////////////////
//
//  The spectral density in some angle interval from varAngle1 to
//  varAngle2
//

G4double G4ForwardXrayTR::AngleInterval( G4double    energy,
                                         G4double varAngle1,
                                         G4double varAngle2   ) const
{
  return     SpectralDensity(energy,varAngle2)
           - SpectralDensity(energy,varAngle1) ;
}

////////////////////////////////////////////////////////////////////
//
// Integral spectral distribution of X-ray TR photons based on
// analytical formula for spectral density
//

G4double G4ForwardXrayTR::EnergySum( G4double energy1,
                                     G4double energy2     )   const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(energy2 - energy1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
   sumEven += AngleInterval(energy1 + 2*i*h,0.0,fMaxThetaTR);
   sumOdd  += AngleInterval(energy1 + (2*i - 1)*h,0.0,fMaxThetaTR) ;
  }
  sumOdd += AngleInterval(energy1 + (2*fSympsonNumber - 1)*h,
                          0.0,fMaxThetaTR) ;

  return h*(  AngleInterval(energy1,0.0,fMaxThetaTR)
            + AngleInterval(energy2,0.0,fMaxThetaTR)
            + 4.0*sumOdd + 2.0*sumEven                          )/3.0 ;
}

/////////////////////////////////////////////////////////////////////////
//
// PostStepDoIt function for creation of forward X-ray photons in TR process
// on boubndary between two materials with really different plasma energies
//

G4VParticleChange* G4ForwardXrayTR::PostStepDoIt(const G4Track& aTrack,
				                const G4Step&  aStep)
{
  aParticleChange.Initialize(aTrack);
  //  G4cout<<"call G4ForwardXrayTR::PostStepDoIt"<<G4endl ;
  G4int iMat, jMat, iTkin, iPlace, numOfMat, numOfTR, iTR, iTransfer ;

  G4double energyPos, anglePos, energyTR, theta, phi, dirX, dirY, dirZ ;
  G4double W, W1, W2, E1, E2 ;
  static 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
  numOfMat = theMaterialTable->length() ;

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  if (pPostStepPoint->GetStepStatus() != fGeomBoundary)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if (aTrack.GetStepLength() <= kCarTolerance*0.5)

  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  // Come on boundary, so begin to try TR

  iMat = pPreStepPoint ->GetPhysicalVolume()->
			 GetLogicalVolume()->GetMaterial()->GetIndex() ;
  jMat = pPostStepPoint->GetPhysicalVolume()->
			 GetLogicalVolume()->GetMaterial()->GetIndex() ;

  // The case of equal or approximate (in terms of plasma energy) materials
  // No TR photons ?!

  if (     iMat == jMat
      || (    (fMatIndex1 >= 0 && fMatIndex1 >= 0) 
           && ( iMat != fMatIndex1 && iMat != fMatIndex2 )
           && ( jMat != fMatIndex1 && jMat != fMatIndex2 )  )

      || (*theMaterialTable)(iMat)->GetState() ==
         (*theMaterialTable)(jMat)->GetState()
  
      ||(    (*theMaterialTable)(iMat)->GetState() == kStateSolid
          && (*theMaterialTable)(jMat)->GetState() == kStateLiquid )
 
      ||(    (*theMaterialTable)(iMat)->GetState() == kStateLiquid
          && (*theMaterialTable)(jMat)->GetState() == kStateSolid  )   )
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep) ;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
 
  if(charge == 0.0) // Uncharged particle doesn't Generate TR photons 
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  // Now we are ready to Generate TR photons

  G4double chargeSq = charge*charge ;
  G4double kinEnergy     = aParticle->GetKineticEnergy() ;
  G4double massRatio = proton_mass_c2/aParticle->GetDefinition()->GetPDGMass() ;
  G4double TkinScaled = kinEnergy*massRatio ;
  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin)) // <= ?
    {
      break ;
    }     
  }
  if(jMat < iMat)
  {
    iPlace = fTotBin + iTkin - 1 ; // (iMat*(numOfMat - 1) + jMat)*
  }
  else
  {
    iPlace = iTkin - 1 ;  // (iMat*(numOfMat - 1) + jMat - 1)*fTotBin + 
  }
  //  G4PhysicsVector*  energyVector1 = (*fEnergyDistrTable)(iPlace)     ;
  //  G4PhysicsVector*  energyVector2 = (*fEnergyDistrTable)(iPlace + 1) ;

  //  G4PhysicsVector*   angleVector1 = (*fAngleDistrTable)(iPlace)      ;
  //  G4PhysicsVector*   angleVector2 = (*fAngleDistrTable)(iPlace + 1)  ;

  G4ParticleMomentum particleDir = aParticle->GetMomentumDirection() ;

  if(iTkin == fTotBin)                 // TR plato, try from left
  {
 // G4cout<<iTkin<<" mean TR number = "<<( (*(*fEnergyDistrTable)(iPlace))(0) +
 //                                   (*(*fAngleDistrTable)(iPlace))(0) )
 //      *chargeSq*0.5<<G4endl ;

    numOfTR = RandPoisson::shoot( ( (*(*fEnergyDistrTable)(iPlace))(0) +
                                    (*(*fAngleDistrTable)(iPlace))(0) )
                                    *chargeSq*0.5 ) ;
    if(numOfTR == 0)
    {
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    else
    {
      // G4cout<<"Number of X-ray TR photons = "<<numOfTR<<G4endl ;

      aParticleChange.SetNumberOfSecondaries(numOfTR);

      for(iTR=0;iTR<numOfTR;iTR++)
      {
        energyPos = (*(*fEnergyDistrTable)(iPlace))(0)*G4UniformRand() ;
        for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
	{
          if(energyPos >= (*(*fEnergyDistrTable)(iPlace))(iTransfer)) break ;
	}
        energyTR = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;

	// G4cout<<"energyTR = "<<energyTR/keV<<"keV"<<G4endl ;

        kinEnergy -= energyTR ; 
        aParticleChange.SetEnergyChange(kinEnergy);

        anglePos = (*(*fAngleDistrTable)(iPlace))(0)*G4UniformRand() ;
        for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
	{
          if(anglePos > (*(*fAngleDistrTable)(iPlace))(iTransfer)) break ;
	}
        theta = sqrt((*fAngleDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1)) ;

	// G4cout<<iTransfer<<" :  theta = "<<theta<<G4endl ;

        phi = twopi*G4UniformRand() ;
        dirX = sin(theta)*cos(phi)  ;
        dirY = sin(theta)*sin(phi)  ;
        dirZ = cos(theta)           ;
        G4ThreeVector directionTR(dirX,dirY,dirZ) ;
        directionTR.rotateUz(particleDir) ;
        G4DynamicParticle* aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(),
                                                             directionTR,
                                                             energyTR     ) ;
        aParticleChange.AddSecondary(aPhotonTR) ;
      }
    }
  }
  else
  {
    if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
    {
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    } 
    else          // general case: Tkin between two vectors of the material
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;

  // G4cout<<iTkin<<" mean TR number = "<<(((*(*fEnergyDistrTable)(iPlace))(0)+
  // (*(*fAngleDistrTable)(iPlace))(0))*W1 + 
  //                                ((*(*fEnergyDistrTable)(iPlace + 1))(0)+
  // (*(*fAngleDistrTable)(iPlace + 1))(0))*W2)
  //                                    *chargeSq*0.5<<G4endl ;

      numOfTR = RandPoisson::shoot((((*(*fEnergyDistrTable)(iPlace))(0)+
(*(*fAngleDistrTable)(iPlace))(0))*W1 + 
                                    ((*(*fEnergyDistrTable)(iPlace + 1))(0)+
(*(*fAngleDistrTable)(iPlace + 1))(0))*W2)
                                    *chargeSq*0.5 ) ;
      if(numOfTR == 0)
      {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      else
      {
        // G4cout<<"Number of X-ray TR photons = "<<numOfTR<<G4endl ;

        aParticleChange.SetNumberOfSecondaries(numOfTR);
        for(iTR=0;iTR<numOfTR;iTR++)
        {
          energyPos = ((*(*fEnergyDistrTable)(iPlace))(0)*W1+
                       (*(*fEnergyDistrTable)(iPlace + 1))(0)*W2)*G4UniformRand() ;
          for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
  	  {
            if(energyPos >= ((*(*fEnergyDistrTable)(iPlace))(iTransfer)*W1+
                       (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer)*W2)) break ;
   	  }
          energyTR = ((*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer))*W1+
               ((*fEnergyDistrTable)(iPlace + 1)->GetLowEdgeEnergy(iTransfer))*W2 ;

	  // G4cout<<"energyTR = "<<energyTR/keV<<"keV"<<G4endl ;

          kinEnergy -= energyTR ; 
          aParticleChange.SetEnergyChange(kinEnergy);

          anglePos = ((*(*fAngleDistrTable)(iPlace))(0)*W1+
                       (*(*fAngleDistrTable)(iPlace + 1))(0)*W2)*G4UniformRand() ;
          for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
	  {
            if(anglePos > ((*(*fAngleDistrTable)(iPlace))(iTransfer)*W1+
                      (*(*fAngleDistrTable)(iPlace + 1))(iTransfer)*W2)) break ;
	  }
          theta = sqrt(((*fAngleDistrTable)(iPlace)->
                        GetLowEdgeEnergy(iTransfer-1))*W1+
                  ((*fAngleDistrTable)(iPlace + 1)->
                        GetLowEdgeEnergy(iTransfer-1))*W2) ;

	  // G4cout<<iTransfer<<" : theta = "<<theta<<G4endl ;

          phi = twopi*G4UniformRand() ;
          dirX = sin(theta)*cos(phi)  ;
          dirY = sin(theta)*sin(phi)  ;
          dirZ = cos(theta)           ;
          G4ThreeVector directionTR(dirX,dirY,dirZ) ;
          directionTR.rotateUz(particleDir) ;
          G4DynamicParticle* aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(),
                                                               directionTR,
                                                               energyTR     ) ;
          aParticleChange.AddSecondary(aPhotonTR) ;
        }
      }
    }
  }
  return &aParticleChange ;
}

////////////////////////////////////////////////////////////////////////////
//
// Test function for checking of PostStepDoIt random preparation of TR photon
// energy
//

G4double 
G4ForwardXrayTR::GetEnergyTR(G4int iMat, G4int jMat, G4int iTkin) const
{
  G4int  iPlace, numOfMat, numOfTR, iTR, iTransfer ;
  G4double energyTR = 0.0 ; // return this value for no TR photons
  G4double energyPos  ;
  G4double  W1, W2, E1, E2 ;

  static 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
  numOfMat = theMaterialTable->length() ;


  // The case of equal or approximate (in terms of plasma energy) materials
  // No TR photons ?!


  if (     iMat == jMat

      || (*theMaterialTable)(iMat)->GetState() ==
         (*theMaterialTable)(jMat)->GetState()
  
      ||(    (*theMaterialTable)(iMat)->GetState() == kStateSolid
          && (*theMaterialTable)(jMat)->GetState() == kStateLiquid )
 
      ||(    (*theMaterialTable)(iMat)->GetState() == kStateLiquid
          && (*theMaterialTable)(jMat)->GetState() == kStateSolid  )   )
  {
    return energyTR ;
  }

  if(jMat < iMat)
  {
    iPlace = (iMat*(numOfMat - 1) + jMat)*fTotBin + iTkin - 1 ;
  }
  else
  {
    iPlace = (iMat*(numOfMat - 1) + jMat - 1)*fTotBin + iTkin - 1 ;
  }
  G4PhysicsVector*  energyVector1 = (*fEnergyDistrTable)(iPlace)     ;
  G4PhysicsVector*  energyVector2 = (*fEnergyDistrTable)(iPlace + 1) ;

  if(iTkin == fTotBin)                 // TR plato, try from left
  {
    numOfTR = RandPoisson::shoot( (*energyVector1)(0)  ) ;
    if(numOfTR == 0)
    {
      return energyTR ;
    }
    else
    {
      for(iTR=0;iTR<numOfTR;iTR++)
      {
        energyPos = (*energyVector1)(0)*G4UniformRand() ;
        for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
	{
          if(energyPos >= (*energyVector1)(iTransfer)) break ;
	}
        energyTR += energyVector1->GetLowEdgeEnergy(iTransfer) ;
      }
    }
  }
  else
  {
    if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
    {
      return energyTR ;
    } 
    else          // general case: Tkin between two vectors of the material
    {             // use trivial mean half/half
      W1 = 0.5 ; 
      W2 = 0.5 ;
     numOfTR = RandPoisson::shoot( (*energyVector1)(0)*W1 +
                                   (*energyVector2)(0)*W2  ) ;
      if(numOfTR == 0)
      {
        return energyTR ;
      }
      else
      {
  G4cout<<"It is still OK in GetEnergyTR(int,int,int)"<<G4endl;
        for(iTR=0;iTR<numOfTR;iTR++)
        {
          energyPos = ((*energyVector1)(0)*W1+
                       (*energyVector2)(0)*W2)*G4UniformRand() ;
          for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
  	  {
            if(energyPos >= ((*energyVector1)(iTransfer)*W1+
                             (*energyVector2)(iTransfer)*W2)) break ;
   	  }
          energyTR += (energyVector1->GetLowEdgeEnergy(iTransfer))*W1+
                      (energyVector2->GetLowEdgeEnergy(iTransfer))*W2 ;

        }
      }
    }
  }

  return energyTR   ;
}

////////////////////////////////////////////////////////////////////////////
//
// Test function for checking of PostStepDoIt random preparation of TR photon
// theta angle relative to particle direction
//
 
       
G4double 
G4ForwardXrayTR::GetThetaTR(G4int iMat, G4int jMat, G4int iTkin) const     
{
  G4double theta = 0.0 ;

  return theta   ;
}



// end of G4ForwardXrayTR implementation file 
//
///////////////////////////////////////////////////////////////////////////
