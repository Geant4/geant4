// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IrregularXrayTRmodel.cc,v 1.2 2001-02-27 07:39:40 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//



#include "G4IrregularXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4IrregularXrayTRmodel::G4IrregularXrayTRmodel(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXrayTRadModel(anEnvelope,a,b)
{
  G4cout<<"Irregular X-ray TR model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

    BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4IrregularXrayTRmodel::~G4IrregularXrayTRmodel()
{
  ;
}



//////////////////////////////////////////////////////////////////////////
//
// User method to code the parameterisation properly
// said.
//
/* *******************************************************

void G4IrregularXrayTRmodel::DoIt( const G4FastTrack& fastTrack , 
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

  G4ParticleMomentum particleDir = fastTrack.GetPrimaryTrack()->
                                   GetMomentumDirection() ;

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
        energyPos = ((*(*fEnergyDistrTable)(iPlace))(0)*W1+
                       (*(*fEnergyDistrTable)(iPlace + 1))(0)*W2)*G4UniformRand() ;
        for(iTransfer=0;iTransfer<fBinTR-1;iTransfer++)
  	{
          if(energyPos >= ((*(*fEnergyDistrTable)(iPlace))(iTransfer)*W1+
                       (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer)*W2)) break ;
   	}
        energyTR = ((*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer))*W1+
               ((*fEnergyDistrTable)(iPlace + 1)->GetLowEdgeEnergy(iTransfer))*W2 ;

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
        directionTR.rotateUz(particleDir) ;
        directionTR.unit() ;

        G4DynamicParticle aPhotonTR(G4Gamma::Gamma(),directionTR,energyTR) ;

        G4ThreeVector positionTR = fastTrack.GetPrimaryTrackLocalPosition() +
	                           G4UniformRand()*distance*direction ;


        G4double distanceTR = fastTrack.GetEnvelopeSolid()->
                              DistanceToOut(positionTR,directionTR) ;

        positionTR = positionTR + distanceTR*directionTR ;

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

***************************************************** */


///////////////////////////////////////////////////////////////////////////
//
// Very rough approximation for radiator interference factor for the case of
// fully irregular radiator. The plate and gas gap thicknesses are distributed 
// according to exponent. The mean values of the plate and gas gap thicknesses 
// are supposed to be much more than XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4IrregularXrayTRmodel::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, plateFactor, gasFactor, factor ;

  plateFactor = 1.0/( 1.0 + fPlateThick*GetPlateLinearPhotoAbs(energy) ) ;

  gasFactor   = 1.0/( 1.0 + fGasThick*GetGasLinearPhotoAbs(energy) ) ;

  factor      = plateFactor*gasFactor ;

  result      = ( 1 - pow(factor,fPlateNumber) )/( 1 - factor ) ;

  result     *= 1 + plateFactor ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








