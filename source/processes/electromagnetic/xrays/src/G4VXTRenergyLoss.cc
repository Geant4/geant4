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
// $Id: G4VXTRenergyLoss.cc,v 1.21 2005/10/11 14:24:34 grichine Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// History:
// 2001-2002 R&D by V.Grichine
// 19.06.03 V. Grichine, modifications in BuildTable for the integration 
//                       in respect of angle: range is increased, accuracy is
//                       improved
// 28.07.05, P.Gumplinger add G4ProcessType to constructor
//

#include "G4Timer.hh"

#include "G4VXTRenergyLoss.hh"
#include "G4Poisson.hh"
#include "G4MaterialTable.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4VSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"

#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"

using namespace std;

// Initialization of local constants

G4double G4XTRenergyLoss::fTheMinEnergyTR   =    1.0*keV  ;
G4double G4XTRenergyLoss::fTheMaxEnergyTR   =  100.0*keV  ;
G4double G4XTRenergyLoss::fTheMaxAngle    =      1.0e-3   ;
G4double G4XTRenergyLoss::fTheMinAngle    =      5.0e-6   ;
G4int    G4XTRenergyLoss::fBinTR            =   50        ;

G4double G4XTRenergyLoss::fMinProtonTkin = 100.0*GeV  ;
G4double G4XTRenergyLoss::fMaxProtonTkin = 100.0*TeV  ;
G4int    G4XTRenergyLoss::fTotBin        =  50        ;
// Proton energy vector initialization

G4PhysicsLogVector* G4XTRenergyLoss::
fProtonEnergyVector = new G4PhysicsLogVector(fMinProtonTkin,
                                             fMaxProtonTkin,
                                                    fTotBin  ) ;

G4double G4XTRenergyLoss::fPlasmaCof = 4.0*pi*fine_structure_const*
                                       hbarc*hbarc*hbarc/electron_mass_c2 ;

G4double G4XTRenergyLoss::fCofTR     = fine_structure_const/pi ;





////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4XTRenergyLoss::G4XTRenergyLoss(G4LogicalVolume *anEnvelope,
				   G4Material* foilMat,G4Material* gasMat,
                                    G4double a, G4double b,
                                    G4int n,const G4String& processName,
                                    G4ProcessType type) :
  G4VDiscreteProcess(processName, type)
  // G4VContinuousProcess(processName, type)
{
  fEnvelope = anEnvelope ;
  //  fPlateNumber = fEnvelope->GetNoDaughters() ;
  fPlateNumber = n ;
  G4cout<<"the number of TR radiator plates = "<<fPlateNumber<<G4endl ;
  if(fPlateNumber == 0)
  {
    G4Exception("No plates in X-ray TR radiator") ;
  }
  // default is XTR dEdx, not flux after radiator
  fExitFlux = false;
  fLambda = DBL_MAX;
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

  pParticleChange = &fParticleChange;

}

///////////////////////////////////////////////////////////////////////////

G4XTRenergyLoss::~G4XTRenergyLoss()
{
   G4int i ;

   if(fEnvelope) delete fEnvelope;

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


G4bool G4XTRenergyLoss::IsApplicable(const G4ParticleDefinition& particle)
{
  return  ( particle.GetPDGCharge() != 0.0 ) ;
}

//////////////////////////////////////////////////////////////////////////////////
//
// GetContinuousStepLimit
//

G4double
G4XTRenergyLoss::GetContinuousStepLimit(const G4Track& ,
				         G4double  ,
				         G4double  ,
                                         G4double& )
{
	G4double StepLimit = DBL_MAX;

	return StepLimit;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calculate step size for XTR process inside raaditor

G4double G4XTRenergyLoss::GetMeanFreePath(const G4Track& aTrack,
					  G4double, // previousStepSize,
                           G4ForceCondition* condition)
{
  G4int iTkin, iPlace;
  G4double lambda, sigma, kinEnergy, mass, gamma;
  G4double charge, chargeSq, massRatio, TkinScaled;
  G4double E1,E2,W,W1,W2;

 *condition = NotForced;
  
  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) lambda = DBL_MAX;
  else
  {
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    kinEnergy = aParticle->GetKineticEnergy();
    mass      = aParticle->GetDefinition()->GetPDGMass();
    gamma     = 1.0 + kinEnergy/mass;
    if(verboseLevel)
    {
      G4cout<<" gamma = "<<gamma<<";   fGamma = "<<fGamma<<G4endl;
    }

    if ( fabs( gamma - fGamma ) < 0.05*gamma ) lambda = fLambda;
    else
    {
      charge = aParticle->GetDefinition()->GetPDGCharge();
      chargeSq  = charge*charge;
      massRatio = proton_mass_c2/mass;
      TkinScaled = kinEnergy*massRatio;

      for(iTkin = 0; iTkin < fTotBin; iTkin++)
      {
        if( TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break ;    
      }
      iPlace = iTkin - 1 ;

      if(iTkin == 0) lambda = DBL_MAX; // Tkin is too small, neglect of TR photon generation
      else          // general case: Tkin between two vectors of the material
      {
        if(iTkin == fTotBin) 
        {
          sigma = (*(*fEnergyDistrTable)(iPlace))(0)*chargeSq;
        }
        else
        {
          E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
          E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
           W = 1.0/(E2 - E1) ;
          W1 = (E2 - TkinScaled)*W ;
          W2 = (TkinScaled - E1)*W ;
          sigma = ( (*(*fEnergyDistrTable)(iPlace  ))(0)*W1 +
                (*(*fEnergyDistrTable)(iPlace+1))(0)*W2   )*chargeSq;
      
        }
        if (sigma < DBL_MIN) lambda = DBL_MAX;
        else                 lambda = 1./sigma; 
        fLambda = lambda;
        fGamma  = gamma;   
        if(verboseLevel)
        {
	  G4cout<<" lambda = "<<lambda/mm<<" mm"<<G4endl;
        }
      }
    }
  }  
  return lambda;
}


//////////////////////////////////////////////////////////////////////////
//
// Build integral energy distribution of XTR photons

void G4XTRenergyLoss::BuildTable()
{
  G4int iTkin, iTR, iPlace;
  G4double radiatorCof = 1.0;           // for tuning of XTR yield

  fEnergyDistrTable = new G4PhysicsTable(fTotBin);

  fGammaTkinCut = 0.0;
  
  // setting of min/max TR energies 
  
  if(fGammaTkinCut > fTheMinEnergyTR)  fMinEnergyTR = fGammaTkinCut ;
  else                                 fMinEnergyTR = fTheMinEnergyTR ;
	
  if(fGammaTkinCut > fTheMaxEnergyTR) fMaxEnergyTR = 2.0*fGammaTkinCut ;  
  else                                fMaxEnergyTR = fTheMaxEnergyTR ;

  G4cout.precision(4) ;
  G4Timer timer ;
  timer.Start() ;
  G4cout<<G4endl;
  G4cout<<"Lorentz Factor"<<"\t"<<"XTR photon number"<<G4endl;
  G4cout<<G4endl;
	
  for( iTkin = 0 ; iTkin < fTotBin ; iTkin++ )      // Lorentz factor loop
  {
     G4PhysicsLogVector* energyVector = new G4PhysicsLogVector( fMinEnergyTR,
                                                                fMaxEnergyTR,
                                                                fBinTR         ) ;

     fGamma = 1.0 + (fProtonEnergyVector->
                            GetLowEdgeEnergy(iTkin)/proton_mass_c2) ;

     fMaxThetaTR = 25.0/(fGamma*fGamma) ;  // theta^2

     fTheMinAngle = 1.0e-3 ; // was 5.e-6, e-6 !!!, e-5, e-4
 
     if( fMaxThetaTR > fTheMaxAngle )    fMaxThetaTR = fTheMaxAngle; 
     else
     {
       if( fMaxThetaTR < fTheMinAngle )  fMaxThetaTR = fTheMinAngle;
     }
G4PhysicsLinearVector* angleVector = new G4PhysicsLinearVector(0.0,
                                                               fMaxThetaTR,
                                                               fBinTR      );

     G4double energySum = 0.0;
     G4double angleSum  = 0.0;

G4Integrator<G4XTRenergyLoss,G4double(G4XTRenergyLoss::*)(G4double)> integral;

     energyVector->PutValue(fBinTR-1,energySum);
     angleVector->PutValue(fBinTR-1,angleSum);

     for( iTR = fBinTR - 2 ; iTR >= 0 ; iTR-- )
     {
        energySum += radiatorCof*fCofTR*integral.Legendre10(
		     this,&G4XTRenergyLoss::SpectralXTRdEdx,
                     energyVector->GetLowEdgeEnergy(iTR),
                     energyVector->GetLowEdgeEnergy(iTR+1) ); 

	//    angleSum  += fCofTR*integral.Legendre96(
	//       this,&G4XTRenergyLoss::AngleXTRdEdx,
	//       angleVector->GetLowEdgeEnergy(iTR),
	//       angleVector->GetLowEdgeEnergy(iTR+1) );

        energyVector->PutValue(iTR,energySum/fTotalDist);
        //  angleVector ->PutValue(iTR,angleSum);
     }
     G4cout
       // <<iTkin<<"\t"
       //   <<"fGamma = "
       <<fGamma<<"\t"  //  <<"  fMaxThetaTR = "<<fMaxThetaTR
       //  <<"sumN = "
       <<energySum      // <<" ; sumA = "<<angleSum
       <<G4endl;
     iPlace = iTkin;
     fEnergyDistrTable->insertAt(iPlace,energyVector);
     //  fAngleDistrTable->insertAt(iPlace,angleVector);
  }     
  timer.Stop();
  G4cout.precision(6);
  G4cout<<G4endl;
  G4cout<<"total time for build X-ray TR energy loss tables = "
        <<timer.GetUserElapsed()<<" s"<<G4endl;
  fGamma = 0.;
  return ;
}

//////////////////////////////////////////////////////////////////////////
//
//

void G4XTRenergyLoss::BuildEnergyTable()
{
  return ;
}

////////////////////////////////////////////////////////////////////////
//
//

void G4XTRenergyLoss::BuildAngleTable()
{
  return ;
} 


//////////////////////////////////////////////////////////////////////////////
//
// The main function which is responsible for the treatment of a particle passage
// trough G4Envelope with discrete generation of G4Gamma

G4VParticleChange* G4XTRenergyLoss::PostStepDoIt( const G4Track& aTrack, 
		                                  const G4Step&  aStep   )
{
  G4int iTkin, iPlace;
  G4double energyTR, theta, phi, dirX, dirY, dirZ;
 

  fParticleChange.Initialize(aTrack);

  if(verboseLevel)
  {
    G4cout<<"Start of G4XTRenergyLoss::PostStepDoIt "<<G4endl ;
    G4cout<<"name of current material =  "
          <<aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<G4endl ;
  }
  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) 
  {
    if(verboseLevel)
    {
      G4cout<<"Go out from G4XTRenergyLoss::PostStepDoIt: wrong volume "<<G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  else
  {
    G4StepPoint* pPostStepPoint        = aStep.GetPostStepPoint();
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   
    // Now we are ready to Generate one TR photon

    G4double kinEnergy = aParticle->GetKineticEnergy() ;
    G4double mass      = aParticle->GetDefinition()->GetPDGMass() ;
    G4double gamma     = 1.0 + kinEnergy/mass ;

    if(verboseLevel > 0 )
    {
      G4cout<<"gamma = "<<gamma<<G4endl ;
    }
    G4double         massRatio   = proton_mass_c2/mass ;
    G4double          TkinScaled = kinEnergy*massRatio ;
    G4ThreeVector      position  = pPostStepPoint->GetPosition();
    G4ParticleMomentum direction = aParticle->GetMomentumDirection();
    G4double           startTime = pPostStepPoint->GetGlobalTime();

    for( iTkin = 0; iTkin < fTotBin; iTkin++ )
    {
      if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break;    
    }
    iPlace = iTkin - 1;

    if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
    {
      if( verboseLevel )
      {
        G4cout<<"Go out from G4XTRenergyLoss::PostStepDoIt:iTkin = "<<iTkin<<G4endl;
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    } 
    else          // general case: Tkin between two vectors of the material
    {
      fParticleChange.SetNumberOfSecondaries(1);

      energyTR = GetXTRrandomEnergy(TkinScaled,iTkin);

      if( verboseLevel )
      {
            G4cout<<"energyTR = "<<energyTR/keV<<" keV"<<G4endl;
      }
      theta = fabs(G4RandGauss::shoot(0.0,pi/gamma));

      if( theta >= 0.1 ) theta = 0.1;

      // G4cout<<" : theta = "<<theta<<endl ;

      phi = twopi*G4UniformRand();

      dirX = sin(theta)*cos(phi);
      dirY = sin(theta)*sin(phi);
      dirZ = cos(theta);

      G4ThreeVector directionTR(dirX,dirY,dirZ);
      directionTR.rotateUz(direction);
      directionTR.unit();

      G4DynamicParticle* aPhotonTR = new G4DynamicParticle(G4Gamma::Gamma(),
                                                           directionTR, energyTR);

      // A XTR photon is set on the particle track inside the radiator 
      // and is moved to the G4Envelope surface for standard X-ray TR models
      // only. The case of fExitFlux=true

      if( fExitFlux )
      {
        const G4RotationMatrix* rotM = pPostStepPoint->GetTouchable()->GetRotation();
        G4ThreeVector transl = pPostStepPoint->GetTouchable()->GetTranslation();
        G4AffineTransform transform = G4AffineTransform(rotM,transl);
        transform.Invert();
        G4ThreeVector localP = transform.TransformPoint(position);
        G4ThreeVector localV = transform.TransformAxis(directionTR);

        G4double distance = fEnvelope->GetSolid()->DistanceToOut(localP, localV);
        if(verboseLevel)
        {
          G4cout<<"distance to exit = "<<distance/mm<<" mm"<<G4endl;
        }
        position         += distance*directionTR;
        startTime        += distance/c_light;
      }
      G4Track* aSecondaryTrack = new G4Track( aPhotonTR, 
		                                startTime, position );
      aSecondaryTrack->SetTouchableHandle(
                         aStep.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID( aTrack.GetTrackID() );

      fParticleChange.AddSecondary(aSecondaryTrack);
      fParticleChange.ProposeEnergy(kinEnergy);     
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}



//////////////////////////////////////////////////////////////////////////////
//
// The main function which is responsible for the treatment of a particle passage
// trough G4Envelope

G4VParticleChange* G4XTRenergyLoss::AlongStepDoIt( const G4Track& aTrack, 
		                                    const G4Step&  aStep   )
{
  G4int iTkin, iPlace,   numOfTR, iTR ;
  G4double energyTR, meanNumOfTR, theta, phi, dirX, dirY, dirZ, rand ;
  G4double W, W1, W2, E1, E2 ;

  aParticleChange.Initialize(aTrack);

  if(verboseLevel)
  {
    G4cout<<"Start of G4XTRenergyLoss::AlongStepDoIt "<<G4endl ;
    G4cout<<"name of current material =  "
          <<aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<G4endl ;
  }
// if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume() != fEnvelope) 

  if( aTrack.GetVolume()->GetLogicalVolume() != fEnvelope ) 
  {
    if(verboseLevel)
    {
      G4cout<<"Go out from G4XTRenergyLoss::AlongStepDoIt: wrong volume "<<G4endl;
    }
    //  return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
  }
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
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break;    
  }
  iPlace = iTkin - 1 ;

  if(iTkin == 0) // Tkin is too small, neglect of TR photon generation
  {
    if(verboseLevel)
    {
      G4cout<<"Go out from G4XTRenergyLoss::AlongStepDoIt:iTkin = "<<iTkin<<G4endl;
    }
    //  return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
  } 
  else          // general case: Tkin between two vectors of the material
  {
    if(iTkin == fTotBin) 
    {
      meanNumOfTR = (*(*fEnergyDistrTable)(iPlace))(0)*chargeSq*distance ;
      numOfTR = G4Poisson(meanNumOfTR) ;
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
              <<G4endl ;
      }
      numOfTR = G4Poisson( meanNumOfTR ) ;
    }
    if( numOfTR == 0 ) // no change, return 
    {
      aParticleChange.SetNumberOfSecondaries(0);
      if(verboseLevel)
      {
      G4cout<<"Go out from G4XTRenergyLoss::AlongStepDoIt: numOfTR = "
            <<numOfTR<<G4endl ;
      }
      // return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep); 
    }
    else
    {
      if(verboseLevel)
      {
        G4cout<<"Number of X-ray TR photons = "<<numOfTR<<G4endl ;
      }
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

      if(verboseLevel)
      {
        G4cout<<"energyTR = "<<energyTR/keV<<"keV"<<G4endl ;
      }
      sumEnergyTR += energyTR ;

      theta = fabs(G4RandGauss::shoot(0.0,pi/gamma)) ;

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

        if( fExitFlux )
        {
          const G4RotationMatrix* rotM = pPostStepPoint->GetTouchable()->GetRotation();
          G4ThreeVector transl = pPostStepPoint->GetTouchable()->GetTranslation();
          G4AffineTransform transform = G4AffineTransform(rotM,transl);
          transform.Invert();
          G4ThreeVector localP = transform.TransformPoint(positionTR);
          G4ThreeVector localV = transform.TransformAxis(directionTR);

          G4double distance = fEnvelope->GetSolid()->DistanceToOut(localP, localV);
          if(verboseLevel)
          {
            G4cout<<"distance to exit = "<<distance/mm<<" mm"<<G4endl;
          }
          positionTR         += distance*directionTR;
          aSecondaryTime        += distance/c_light;
        }
 
        G4Track* aSecondaryTrack = new G4Track( aPhotonTR, 
		                                aSecondaryTime,positionTR ) ;
        aSecondaryTrack->SetTouchableHandle(aStep.GetPostStepPoint()
                                                  ->GetTouchableHandle());
        aSecondaryTrack->SetParentID(aTrack.GetTrackID());

	aParticleChange.AddSecondary(aSecondaryTrack);
      }
      kinEnergy -= sumEnergyTR ;
      aParticleChange.ProposeEnergy(kinEnergy) ;
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

G4complex G4XTRenergyLoss::OneInterfaceXTRdEdx( G4double energy,
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

G4double G4XTRenergyLoss::SpectralAngleXTRdEdx(G4double varAngle)
{
  G4double result =  GetStackFactor(fEnergy,fGamma,varAngle) ;
  if(result < 0.0) result = 0.0 ;
  return result ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double G4XTRenergyLoss::SpectralXTRdEdx(G4double energy)
{
  fEnergy = energy ;
G4Integrator<G4XTRenergyLoss,G4double(G4XTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4XTRenergyLoss::SpectralAngleXTRdEdx,
                             0.0,0.1*fMaxThetaTR) +
         integral.Legendre96(this,&G4XTRenergyLoss::SpectralAngleXTRdEdx,
                             0.1*fMaxThetaTR,0.2*fMaxThetaTR) +         
         integral.Legendre96(this,&G4XTRenergyLoss::SpectralAngleXTRdEdx,
                             0.2*fMaxThetaTR,0.4*fMaxThetaTR) +         
         integral.Legendre96(this,&G4XTRenergyLoss::SpectralAngleXTRdEdx,
                             0.4*fMaxThetaTR,0.7*fMaxThetaTR) +         
         integral.Legendre96(this,&G4XTRenergyLoss::SpectralAngleXTRdEdx,
	                     0.7*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double G4XTRenergyLoss::AngleSpectralXTRdEdx(G4double energy)
{
  G4double result =  GetStackFactor(energy,fGamma,fVarAngle) ;
  if(result < 0) result = 0.0 ;
  return result ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4XTRenergyLoss::AngleXTRdEdx(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator<G4XTRenergyLoss,G4double(G4XTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre10(this,&G4XTRenergyLoss::AngleSpectralXTRdEdx,
			     fMinEnergyTR,fMaxEnergyTR) ;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for plates. Omega is energy !!!

G4double G4XTRenergyLoss::GetPlateFormationZone( G4double omega ,
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

G4complex G4XTRenergyLoss::GetPlateComplexFZ( G4double omega ,
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

void G4XTRenergyLoss::ComputePlatePhotoAbsCof() 
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

G4double G4XTRenergyLoss::GetPlateLinearPhotoAbs(G4double omega) 
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
    G4Exception("Invalid (<I1) energy in G4XTRenergyLoss::GetPlateLinearPhotoAbs");
  }
  else i-- ;
  
  return fPlatePhotoAbsCof[i][1]/omega  + fPlatePhotoAbsCof[i][2]/omega2 + 
         fPlatePhotoAbsCof[i][3]/omega3 + fPlatePhotoAbsCof[i][4]/omega4  ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for gas. Omega is energy !!!

G4double G4XTRenergyLoss::GetGasFormationZone( G4double omega ,
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

G4complex G4XTRenergyLoss::GetGasComplexFZ( G4double omega ,
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

void G4XTRenergyLoss::ComputeGasPhotoAbsCof() 
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

G4double G4XTRenergyLoss::GetGasLinearPhotoAbs(G4double omega) 
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
   G4Exception("Invalid (<I1) energy in G4XTRenergyLoss::GetGasLinearPhotoAbs");
  }
  else i-- ;
  
  return fGasPhotoAbsCof[i][1]/omega  + fGasPhotoAbsCof[i][2]/omega2 + 
         fGasPhotoAbsCof[i][3]/omega3 + fGasPhotoAbsCof[i][4]/omega4  ;

}

//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// Omega is energy !!!

G4double G4XTRenergyLoss::GetPlateZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetPlateFormationZone(omega,gamma,varAngle)*GetPlateLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// G4cout and output in file in some energy range.

void G4XTRenergyLoss::GetPlateZmuProduct() 
{
  ofstream outPlate("plateZmu.dat", ios::out ) ;
  outPlate.setf( ios::scientific, ios::floatfield );

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

G4double G4XTRenergyLoss::GetGasZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetGasFormationZone(omega,gamma,varAngle)*GetGasLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof byformation zone for gas. 
// G4cout and output in file in some energy range.

void G4XTRenergyLoss::GetGasZmuProduct() 
{
  ofstream outGas("gasZmu.dat", ios::out ) ;
  outGas.setf( ios::scientific, ios::floatfield );
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
G4XTRenergyLoss::OneBoundaryXTRNdensity( G4double energy,G4double gamma,
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

G4double G4XTRenergyLoss::GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle )
{
  // return stack factor corresponding to one interface

  return std::real( OneInterfaceXTRdEdx(energy,gamma,varAngle) );
}

//////////////////////////////////////////////////////////////////////////////
//
// For photon energy distribution tables. Integrate first over angle
//

G4double G4XTRenergyLoss::XTRNSpectralAngleDensity(G4double varAngle)
{
  return OneBoundaryXTRNdensity(fEnergy,fGamma,varAngle)*
         GetStackFactor(fEnergy,fGamma,varAngle)             ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double G4XTRenergyLoss::XTRNSpectralDensity(G4double energy)
{
  fEnergy = energy ;
  G4Integrator<G4XTRenergyLoss,G4double(G4XTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4XTRenergyLoss::XTRNSpectralAngleDensity,
                             0.0,0.2*fMaxThetaTR) +
         integral.Legendre10(this,&G4XTRenergyLoss::XTRNSpectralAngleDensity,
	                     0.2*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double G4XTRenergyLoss::XTRNAngleSpectralDensity(G4double energy)
{
  return OneBoundaryXTRNdensity(energy,fGamma,fVarAngle)*
         GetStackFactor(energy,fGamma,fVarAngle)             ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4XTRenergyLoss::XTRNAngleDensity(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator<G4XTRenergyLoss,G4double(G4XTRenergyLoss::*)(G4double)> integral ;
  return integral.Legendre96(this,&G4XTRenergyLoss::XTRNAngleSpectralDensity,
			     fMinEnergyTR,fMaxEnergyTR) ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Check number of photons for a range of Lorentz factors from both energy 
// and angular tables

void G4XTRenergyLoss::GetNumberOfPhotons()
{
  G4int iTkin ;
  G4double gamma, numberE ;

  ofstream outEn("numberE.dat", ios::out ) ;
  outEn.setf( ios::scientific, ios::floatfield );

  ofstream outAng("numberAng.dat", ios::out ) ;
  outAng.setf( ios::scientific, ios::floatfield );

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

G4double G4XTRenergyLoss::GetXTRrandomEnergy( G4double scaledTkin, G4int iTkin )
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

G4double G4XTRenergyLoss::GetXTRenergy( G4int    iPlace, 
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

