// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhotoClusterModel.cc,v 1.1 2000-11-14 16:08:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4Timer.hh"

#include "G4PhotoClusterModel.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"


////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4PhotoClusterModel::G4PhotoClusterModel(G4Envelope *anEnvelope) :
  G4VClusterModel("G4PhotoClusterModel",anEnvelope)
{
  fMatIndex = anEnvelope->GetMaterial()->GetIndex()  ;
  ComputePhotoAbsCof() ;
  
}

///////////////////////////////////////////////////////////////////////////

G4PhotoClusterModel::~G4PhotoClusterModel()
{
   for(G4int i=0;i<fIntervalNumber;i++)
   {
     delete[] fPhotoAbsCof[i] ;
   }
   delete[] fPhotoAbsCof ;
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns condition for application of the model depending on particle type


G4bool G4PhotoClusterModel::IsApplicable(const G4ParticleDefinition& particle)
{
  return   &particle == G4Gamma::GammaDefinition() ; 
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

G4bool G4PhotoClusterModel::ModelTrigger(const G4FastTrack& fastTrack) 
{
return ( fastTrack.GetPrimaryTrack()->GetKineticEnergy() > 1*keV  && 
         fastTrack.GetPrimaryTrack()->GetKineticEnergy() < 40*keV     ) ;
}

//////////////////////////////////////////////////////////////////////////////
//
// 

void G4PhotoClusterModel::DoIt( const G4FastTrack& fastTrack , 
		                    G4FastStep&  fastStep         )
{
  G4double  energy ;
  G4double distance, lambda, step ;
  G4ThreeVector clusterPosition ;

  energy  = fastTrack.GetPrimaryTrack()->GetKineticEnergy() ;

  G4ParticleMomentum direction(fastTrack.GetPrimaryTrackLocalDirection());

  distance = fastTrack.GetEnvelopeSolid()->
             DistanceToOut(fastTrack.GetPrimaryTrackLocalPosition(),direction) ;

  G4ThreeVector position = fastTrack.GetPrimaryTrackLocalPosition() + 
                           distance*direction ;

  // Set final position:

  fastStep.SetPrimaryTrackFinalPosition(position);

  // Cluster counting loop

  lambda = 1.0/GetLinearPhotoAbs(energy) ;
  step   = RandExponential::shoot(lambda) ;

  if(step > distance) // no change, return 
  {
    return ;  
  }
  else
  {
    G4ThreeVector globalStartPosition  = fastTrack.GetPrimaryTrack()->
                                         GetPosition() ;
    G4ParticleMomentum globalDirection = fastTrack.GetPrimaryTrack()->
                                         GetMomentumDirection() ; 

    // global (or local ?) cluster coordinates
    //  clusterPosition = fastTrack.GetPrimaryTrackLocalPosition() + 
    //                  stepSum*direction ;  
  
    clusterPosition = globalStartPosition  + step*globalDirection ;   
    fClusterPositionVector.insert(clusterPosition) ;      
    fClusterEnergyVector.insert(energy) ;

    fastStep.KillPrimaryTrack();
    fastStep.SetPrimaryTrackPathLength(step);
    fastStep.SetTotalEnergyDeposited(energy);

    BuildDetectorResponse() ; 
  }
  return ;
}

////////////////////////////////////////////////////////////////////////
//
// Computes matrix of Sandia photo absorption cross section coefficients for
// G4Envelope material

void G4PhotoClusterModel::ComputePhotoAbsCof() 
{
  G4int i, j, numberOfElements ;
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4SandiaTable thisMaterialSandiaTable(fMatIndex) ;
  numberOfElements = (*theMaterialTable)[fMatIndex]->GetNumberOfElements() ;
  G4int* thisMaterialZ = new G4int[numberOfElements] ;

  for(i=0;i<numberOfElements;i++)
  {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[fMatIndex]->
                                      GetElement(i)->GetZ() ;
  }
  fIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
  fIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                           (*theMaterialTable)[fMatIndex]->GetFractionVector() ,
        		     numberOfElements,fIntervalNumber) ;
   
  fPhotoAbsCof = new G4double*[fIntervalNumber] ;
  for(i=0;i<fIntervalNumber;i++) fPhotoAbsCof[i] = new G4double[5] ;
   
  for(i=0;i<fIntervalNumber;i++)
  {
      fPhotoAbsCof[i][0] = thisMaterialSandiaTable.
                                GetPhotoAbsorpCof(i+1,0) ; 
                              
      for(j=1;j<5;j++)
      {
           fPhotoAbsCof[i][j] = thisMaterialSandiaTable.
	                             GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex]->GetDensity() ;
      }
  }
  delete[] thisMaterialZ ;
  return ;
}

//////////////////////////////////////////////////////////////////////
//
// Returns the value of linear photo absorption coefficient (in reciprocal 
// length) for G4Envelope material

G4double G4PhotoClusterModel::GetLinearPhotoAbs(G4double omega) 
{
  G4int i ;
  G4double omega2, omega3, omega4 ; 

  omega2 = omega*omega ;
  omega3 = omega2*omega ;
  omega4 = omega2*omega2 ;

  for(i=0;i<fIntervalNumber;i++)
  {
    if( omega < fPhotoAbsCof[i][0] ) break ;
  }
  if( i == 0 )
  { 
   G4Exception("Invalid (<I1) energy in G4PhotoClusterModel::GetLinearPhotoAbs");
  }
  else i-- ;
  
  return fPhotoAbsCof[i][1]/omega  + fPhotoAbsCof[i][2]/omega2 + 
         fPhotoAbsCof[i][3]/omega3 + fPhotoAbsCof[i][4]/omega4  ;

}

//
//
///////////////////////////////////////////////////////////////////////










