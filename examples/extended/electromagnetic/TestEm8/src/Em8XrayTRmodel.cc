// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8XrayTRmodel.cc,v 1.3 2000-04-06 14:40:12 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "Em8XrayTRmodel.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"

// Initialization of local constants

G4double Em8XrayTRmodel::fTheMinEnergyTR   =    1.0*keV  ;
G4double Em8XrayTRmodel::fTheMaxEnergyTR   =  100.0*keV  ;
G4double Em8XrayTRmodel::fTheMaxAngle    =      1.0e-3   ;
G4double Em8XrayTRmodel::fTheMinAngle    =      5.0e-6   ;
G4int    Em8XrayTRmodel::fBinTR            =   50        ;

G4double Em8XrayTRmodel::fMinProtonTkin = 100.0*GeV  ;
G4double Em8XrayTRmodel::fMaxProtonTkin = 100.0*TeV  ;
G4int    Em8XrayTRmodel::fTotBin        =  50        ;
// Proton energy vector initialization

G4PhysicsLogVector* Em8XrayTRmodel::
fProtonEnergyVector = new G4PhysicsLogVector(fMinProtonTkin,
                                             fMaxProtonTkin,
                                                    fTotBin  ) ;

G4double Em8XrayTRmodel::fPlasmaCof = 4.0*pi*fine_structure_const*
                                       hbarc*hbarc*hbarc/electron_mass_c2 ;

G4double Em8XrayTRmodel::fCofTR     = fine_structure_const/pi ;





////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

Em8XrayTRmodel::Em8XrayTRmodel(G4Envelope *anEnvelope, G4double a, G4double b) :
  G4VFastSimulationModel("Em8XrayTRmodel",anEnvelope)
  // ,   G4ForwardXrayTR("Em8XrayTRmodel")
{
  fPlateNumber = anEnvelope->GetNoDaughters() ;
  G4cout<<"the number of TR radiator plates = "<<fPlateNumber<<G4endl ;
  if(fPlateNumber == 0)
  {
    G4Exception("No plates in X-ray TR radiator") ;
  }
  // Mean thicknesses of plates and gas gaps

  fPlateThick = a ;
  fGasThick =   b ;

  // index of plate material
  fMatIndex1 = anEnvelope->GetDaughter(0)->GetLogicalVolume()->
                           GetMaterial()->GetIndex()  ;
  // index of gas material
  fMatIndex2 = anEnvelope->GetMaterial()->GetIndex()  ;

  // plasma energy squared for plate material

  fSigma1 = fPlasmaCof*anEnvelope->GetDaughter(0)->GetLogicalVolume()->
                           GetMaterial()->GetElectronDensity()  ;
  //  fSigma1 = (20.9*eV)*(20.9*eV) ;
  G4cout<<"plate plasma energy = "<<sqrt(fSigma1)/eV<<" eV"<<G4endl ;

  // plasma energy squared for gas material

  fSigma2 = fPlasmaCof*anEnvelope->GetMaterial()->GetElectronDensity()  ;
  G4cout<<"gas plasma energy = "<<sqrt(fSigma2)/eV<<" eV"<<G4endl ;

  // Compute cofs for preparation of linear photo absorption

  ComputePlatePhotoAbsCof() ;
  ComputeGasPhotoAbsCof() ;

}

///////////////////////////////////////////////////////////////////////////

Em8XrayTRmodel::~Em8XrayTRmodel()
{
  ;
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns condition for application of the model depending on particle type


G4bool Em8XrayTRmodel::IsApplicable(const G4ParticleDefinition& particle)
{
  return  ( particle.GetPDGCharge() != 0.0 ) ; 
}

/////////////////////////////////////////////////////////////////////
//
// UserTrigger() method: method which has to decide if
// the parameterisation has to be applied.
// Here ModelTrigger() asks the user (ie you) a 0/1 answer.
//
// Note that quantities like the local/global position/direction etc..
// are available at this level via the fastTrack parameter (allowing 
// to check distance from boundaries, see below  to allow the decision)
//

G4bool Em8XrayTRmodel::ModelTrigger(const G4FastTrack& fastTrack) 
{
  //  G4double mass = fastTrack.GetPrimaryTrack()->GetDefinition()->GetPDGMass() ;
  //  G4double kinEnergy = fastTrack.GetPrimaryTrack()->GetKineticEnergy() ;
  //  G4double gamma = 1.0 + kinEnergy/mass ;  // Lorentz factor
  // G4cout << "gamma = " << gamma << G4endl ;
  //  if (gamma >= 100.0) return true  ;
  //  else                return false ;

  return true ;
}

//////////////////////////////////////////////////////////////////////////////
//
// 

void Em8XrayTRmodel::DoIt( const G4FastTrack& fastTrack , 
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

        theta = abs(RandGauss::shoot(0.0,pi/gamma)) ;

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

//////////////////////////////////////////////////////////////////////////
//
// User method to code the parameterisation properly
// said. This is simple example of creation of one X-ray photon with the 
// energy in the range of around 5 keV produced by relativistic charged 
// particle
//

void Em8XrayTRmodel::ExampleDoIt( const G4FastTrack& fastTrack , 
		                        G4FastStep&  fastStep         )
{

  // The primary track continues along its direction.
  // One secondary (a photon) is added:

  //  G4cout << "      TR `model' applied \n " << endl;

  // Primary: idem as in "DefaultModel":
  //

  G4double distance = fastTrack.GetEnvelopeSolid()->
    DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),
		  fastTrack.GetPrimaryTrackLocalDirection()) ;

  G4ThreeVector position = fastTrack.GetPrimaryTrackLocalPosition() + 
                           distance*fastTrack.GetPrimaryTrackLocalDirection() ;

  // Set final position:

  fastStep.SetPrimaryTrackFinalPosition(position);

  //---------------------------
  // Secondary: Adds one "secondary":
  //
  // First, user has to say how many secondaries will be created:

  fastStep.SetNumberOfSecondaryTracks(1);

  //  Build the secondary direction:

  G4ParticleMomentum direction(fastTrack.GetPrimaryTrackLocalDirection());
  //  direction.setZ(direction.z()*0.5);
  //  direction.setY(direction.y()+direction.z()*0.1);
  direction = direction.unit();                        // necessary !?

  // Dynamics (Note that many constructors exists for G4DynamicParticle

  G4double gammaEnergy = 3.0*keV + G4UniformRand()*2*keV ;
  
  G4DynamicParticle dynamique(G4Gamma::GammaDefinition(),
			      direction,
	//   fastTrack.GetPrimaryTrack()->GetKineticEnergy()/2.
                      gammaEnergy        );
  // -- position:

  G4double Dist = fastTrack.GetEnvelopeSolid()->
               DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),direction) ;

  G4ThreeVector posi = fastTrack.GetPrimaryTrackLocalPosition() + Dist*direction ;
  
  // Creation of the secondary Track:
 
  fastStep.CreateSecondaryTrack( dynamique, 
                                 posi, 
		                 fastTrack.GetPrimaryTrack()->GetGlobalTime());

}

//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for plates. Omega is energy !!!

G4double Em8XrayTRmodel::GetPlateFormationZone( G4double omega ,
                                                G4double gamma ,
                                                G4double varAngle    ) 
{
  G4double cof, lambda ;
  lambda = 1.0/gamma/gamma + varAngle + fSigma1/omega/omega ;
  cof = 2.0*hbarc/omega/lambda ;
  return cof ;
}

////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// plate material

void Em8XrayTRmodel::ComputePlatePhotoAbsCof() 
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

G4double Em8XrayTRmodel::GetPlateLinearPhotoAbs(G4double omega) 
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
    G4Exception("Invalid (<I1) energy in Em8XrayTRmodel::GetPlateLinearPhotoAbs");
  }
  else i-- ;
  
  return fPlatePhotoAbsCof[i][1]/omega  + fPlatePhotoAbsCof[i][2]/omega2 + 
         fPlatePhotoAbsCof[i][3]/omega3 + fPlatePhotoAbsCof[i][4]/omega4  ;
}

//////////////////////////////////////////////////////////////////////
//
// Calculates formation zone for gas. Omega is energy !!!

G4double Em8XrayTRmodel::GetGasFormationZone( G4double omega ,
                                              G4double gamma ,
                                              G4double varAngle   ) 
{
  G4double cof, lambda ;
  lambda = 1.0/gamma/gamma + varAngle + fSigma2/omega/omega ;
  cof = 2.0*hbarc/omega/lambda ;
  return cof ;

}

////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// gas material

void Em8XrayTRmodel::ComputeGasPhotoAbsCof() 
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

G4double Em8XrayTRmodel::GetGasLinearPhotoAbs(G4double omega) 
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
   G4Exception("Invalid (<I1) energy in Em8XrayTRmodel::GetGasLinearPhotoAbs");
  }
  else i-- ;
  
  return fGasPhotoAbsCof[i][1]/omega  + fGasPhotoAbsCof[i][2]/omega2 + 
         fGasPhotoAbsCof[i][3]/omega3 + fGasPhotoAbsCof[i][4]/omega4  ;

}

//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// Omega is energy !!!

G4double Em8XrayTRmodel::GetPlateZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetPlateFormationZone(omega,gamma,varAngle)*GetPlateLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof by formation zone for plate. 
// G4cout and output in file in some energy range.

void Em8XrayTRmodel::GetPlateZmuProduct() 
{
  ofstream outPlate("plateZmu.dat", ios::out ) ;
  outPlate.setf( ios::scientific, ios::floatfield );

  G4int i ;
  G4double omega, varAngle, gamma, result ;
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

G4double Em8XrayTRmodel::GetGasZmuProduct( G4double omega ,
                                             G4double gamma ,
                                             G4double varAngle   ) 
{
  return GetGasFormationZone(omega,gamma,varAngle)*GetGasLinearPhotoAbs(omega) ;
}
//////////////////////////////////////////////////////////////////////
//
// Calculates the product of linear cof byformation zone for gas. 
// G4cout and output in file in some energy range.

void Em8XrayTRmodel::GetGasZmuProduct() 
{
  ofstream outGas("gasZmu.dat", ios::out ) ;
  outGas.setf( ios::scientific, ios::floatfield );
  G4int i ;
  G4double omega, varAngle, gamma, result ;
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
Em8XrayTRmodel::OneBoundaryXTRNdensity( G4double energy,G4double gamma,
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

//////////////////////////////////////////////////////////////////////////
//
//

void Em8XrayTRmodel::BuildTable() 
{
  G4int iMat, jMat, iTkin, iTR, iPlace ;
  G4double radiatorCof = 0.5 ;

  // fAngleDistrTable  = new G4PhysicsTable(fTotBin) ;
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
  for(iTkin=0;iTkin<fTotBin;iTkin++)      // Lorentz factor loop
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
     G4Integrator integral ;
     energyVector->PutValue(fBinTR-1,energySum) ;
     angleVector->PutValue(fBinTR-1,angleSum)   ;

     for(iTR=fBinTR-2;iTR>=0;iTR--)
     {
        energySum += radiatorCof*fCofTR*integral.Legendre10(this,
                     &Em8XrayTRmodel::XTRNSpectralDensity, 
                     energyVector->GetLowEdgeEnergy(iTR),
                     energyVector->GetLowEdgeEnergy(iTR+1) ) ; 

	//    angleSum  += fCofTR*integral.Legendre96(this,
	//        &Em8XrayTRmodel::XTRNAngleDensity,
	//       angleVector->GetLowEdgeEnergy(iTR),
	//      angleVector->GetLowEdgeEnergy(iTR+1) ) ;

        energyVector->PutValue(iTR,energySum) ;
        //  angleVector ->PutValue(iTR,angleSum)   ;
     }
     G4cout<<iTkin<<"\t"
           <<"fGamma = "<<fGamma<<"\t"  //  <<"  fMaxThetaTR = "<<fMaxThetaTR
           <<"sumE = "<<energySum      // <<" ; sumA = "<<angleSum
           <<endl ;
     iPlace = iTkin ;
     fEnergyDistrTable->insertAt(iPlace,energyVector) ;
     //  fAngleDistrTable->insertAt(iPlace,angleVector) ;
  }     
  timer.Stop() ;
  G4cout.precision(6) ;
  G4cout<<G4endl ;
  G4cout<<"total time for build X-ray TR tables = "
        <<timer.GetUserElapsed()<<" s"<<G4endl ;
  return ;
}

//////////////////////////////////////////////////////////////////////////
//
//

void Em8XrayTRmodel::BuildEnergyTable()
{
  return ;
}

////////////////////////////////////////////////////////////////////////
//
//

void Em8XrayTRmodel::BuildAngleTable()
{
  return ;
} 

//////////////////////////////////////////////////////////////////////////////
//
// For photon energy distribution tables. Integrate first over angle
//

G4double Em8XrayTRmodel::XTRNSpectralAngleDensity(G4double varAngle)
{
  return OneBoundaryXTRNdensity(fEnergy,fGamma,varAngle)*
         GetStackFactor(fEnergy,fGamma,varAngle)             ;
}

/////////////////////////////////////////////////////////////////////////
//
// For second integration over energy
 
G4double Em8XrayTRmodel::XTRNSpectralDensity(G4double energy)
{
  fEnergy = energy ;
  G4Integrator integral ;
  return integral.Legendre96(this,&Em8XrayTRmodel::XTRNSpectralAngleDensity,
			     0.0,0.2*fMaxThetaTR) +
         integral.Legendre10(this,&Em8XrayTRmodel::XTRNSpectralAngleDensity,
			     0.2*fMaxThetaTR,fMaxThetaTR) ;
} 
 
//////////////////////////////////////////////////////////////////////////
// 
// for photon angle distribution tables
//

G4double Em8XrayTRmodel::XTRNAngleSpectralDensity(G4double energy)
{
  return OneBoundaryXTRNdensity(energy,fGamma,fVarAngle)*
         GetStackFactor(energy,fGamma,fVarAngle)             ;
} 

///////////////////////////////////////////////////////////////////////////
//
//

G4double Em8XrayTRmodel::XTRNAngleDensity(G4double varAngle) 
{
  fVarAngle = varAngle ;
  G4Integrator integral ;
  return integral.Legendre96(this,&Em8XrayTRmodel::XTRNAngleSpectralDensity,
			     fMinEnergyTR,fMaxEnergyTR) ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Check number of photons for a range of Lorentz factors from both energy 
// and angular tables

void Em8XrayTRmodel::GetNumberOfPhotons()
{
  G4int iTkin ;
  G4double gamma, numberE, numberA ;

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

//
//
///////////////////////////////////////////////////////////////////////

