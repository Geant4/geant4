//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Rich advanced example for Geant4
// PadHpdPhotoElectricEffect.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "globals.hh"
#include <cmath>

#include "PadHpdPhotoElectricEffect.hh"
#include "RichTbGeometryParameters.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHandle.hh"
#include "G4GeometryTolerance.hh"
#include "Randomize.hh"
#include "RichTbAnalysisManager.hh"
#include "RichTbRunConfig.hh"
#include "RichTbMaterialParameters.hh"

#include "RichTbAnalysisManager.hh"

PadHpdPhotoElectricEffect::PadHpdPhotoElectricEffect(const G4String& processName ,
   RichTbRunConfig* RConfig)
  :G4VDiscreteProcess(processName),
   DemagnificationFactor(std::vector<G4double>(NumHpdTot)),
   DemagnificationQuadFactor(std::vector<G4double>(NumHpdTot)),
   HpdQE(NumHpdTot, std::vector<G4double>( NumQEbins)),
   HpdWabin(NumHpdTot, std::vector<G4double>( NumQEbins))
{
  rConfig=RConfig;
  PrePhotoElectricVolName="PadHpdWindowQuartz";
  PostPhotoElectricVolName="BiAlkaliPhCathode";
  HpdPhElectronKE=(RConfig-> getHpdPhElectronEnergy())*keV;
  PhCathodeToSilDetDist= HpdPhotoCathodeSiZdist;
  PSFsigma=PadHpdPSFsigma;
  
  for(G4int ihpdq=0; ihpdq<NumHpdTot; ihpdq++ ) {
    
    if( HpdDemagLinearTerm[ihpdq] != 0.0 ) {
    DemagnificationFactor[ihpdq]=1.0 / HpdDemagLinearTerm[ihpdq];
    }else { DemagnificationFactor[ihpdq] =1.0 / 2.3;}

    if( HpdDemagQuadraticTerm[ihpdq] != 0.0 ) {

     DemagnificationQuadFactor[ihpdq]=1.0 / HpdDemagQuadraticTerm[ihpdq];

    }else{DemagnificationQuadFactor[ihpdq]=0.0*(1.0/(1.0*mm)); }


    // Now to apply the error on the HPD demag factor. SE 28-4-02
    // for now  a flat error is applied. the uniform number from -1.0 to
    // 1.0 is obtained and then multiplied with the factor. 
    
    G4double DemagError=  (HpdDemagErrorPercent/100.0)*(2.0*G4UniformRand()-1.0) ;
    
    DemagnificationFactor[ihpdq] = DemagnificationFactor[ihpdq]*(1.0+DemagError);


    std::vector<G4double>qeCurHpd =  InitializeHpdQE(ihpdq);
    std::vector<G4double>waCurHpd =  InitializeHpdWaveL(ihpdq);
    if(qeCurHpd.size() != waCurHpd.size() ) {
      G4cout<<"Wrong size for Hpd QE "<<ihpdq<<" "<<qeCurHpd.size()
	    <<"  "<< waCurHpd.size()<<G4endl;
    } 
    for(size_t iqbin=0; iqbin < qeCurHpd.size(); iqbin++){
      HpdQE[ihpdq][iqbin]=qeCurHpd[iqbin]/100;
      HpdWabin[ihpdq][iqbin]=waCurHpd[iqbin];
    }
  } 
  G4cout<<GetProcessName() <<" is created "<<G4endl;

    
  }
PadHpdPhotoElectricEffect::~PadHpdPhotoElectricEffect() {; }

G4bool PadHpdPhotoElectricEffect::IsApplicable(const G4ParticleDefinition& 
					               aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

G4double PadHpdPhotoElectricEffect::GetMeanFreePath(const G4Track& ,
                                              G4double ,
                                              G4ForceCondition* condition)
{
	*condition = Forced;

	return DBL_MAX;
}

G4VParticleChange* PadHpdPhotoElectricEffect::PostStepDoIt(const G4Track& aTrack,
                                                           const G4Step& aStep)
{

  aParticleChange.Initialize(aTrack);

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();


  
  if (pPostStepPoint->GetStepStatus() != fGeomBoundary){

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

   }


  G4String PrePhName = pPreStepPoint -> GetPhysicalVolume() -> 
    GetLogicalVolume() -> GetMaterial()->GetName();
  G4String PostPhName= pPostStepPoint -> GetPhysicalVolume() -> 
    GetLogicalVolume() -> GetMaterial() ->GetName();
  
 
  if(( PrePhName == PrePhotoElectricVolName &&
      PostPhName == PostPhotoElectricVolName) ||    
     ( PostPhName == PrePhotoElectricVolName &&
       PrePhName == PostPhotoElectricVolName) ) {

  
  }else {
    
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

}
  
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  if (aTrack.GetStepLength()<=kCarTolerance/2){
          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  const G4DynamicParticle* aDynamicPhoton = aTrack.GetDynamicParticle();
  G4double PhotonEnergy = aDynamicPhoton->GetKineticEnergy();

  if(PhotonEnergy <= 0.0 ) {

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }


  //Now use the QE for the current HPD to determine if a
  // photoelectron should be produced or not.

  G4int currentHpdNumber= pPreStepPoint->GetTouchableHandle() 
    -> GetReplicaNumber(1);
  if(currentHpdNumber >= NumHpdTot ){

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
           
  }

  G4double PhotWLength=PhotMomWaveConv/(PhotonEnergy/eV);
  
  G4double PhCathodeQE = getHpdQEff(currentHpdNumber, PhotWLength);
  G4double randomnum = G4UniformRand();


  //the following three lines are copied from few lines later just
  // for histogramming convenience.
  G4ThreeVector GlobalElectronOrigin= pPostStepPoint->GetPosition();

  G4Navigator* theNavigator =
             G4TransportationManager::GetTransportationManager()->
                                       GetNavigatorForTracking();
  G4ThreeVector LocalElectronOrigin = theNavigator->
        			GetGlobalToLocalTransform().
				TransformPoint(GlobalElectronOrigin);



  // For the histogram of the radius of the cherenkov circle.
  // This assumes that the beam is along 001 axis in the global
  // coord system.
    G4double GLx=GlobalElectronOrigin.x();
    G4double GLy=GlobalElectronOrigin.y();
    // G4double PhotCkvRad = std::pow((std::pow(GLx,2)+std::pow(GLy,2)),0.5);
    G4double PhotCkvPhi = std::atan2(GLy,GLx)*180.0/pi;

    if( PhotCkvPhi < - 180.0 )PhotCkvPhi+= 360.0;
  
  if(randomnum <  PhCathodeQE ) {

  
  G4double CurDemagFactor=DemagnificationFactor[currentHpdNumber];
  G4double CurDemagQuadFactor=DemagnificationQuadFactor[currentHpdNumber];

  // now get the Point Spread function.

  G4double PsfRandomAzimuth = twopi*G4UniformRand();
  G4double PsfRandomRad= G4RandGauss::shoot(0.0,PSFsigma);
  G4double PsfX= PsfRandomRad*std::cos( PsfRandomAzimuth);
  G4double PsfY= PsfRandomRad*std::sin( PsfRandomAzimuth);

  
  G4ThreeVector LocalElectronDirection(
         (CurDemagFactor+CurDemagQuadFactor*LocalElectronOrigin.x()-1.0)*LocalElectronOrigin.x()+PsfX,
         (CurDemagFactor+CurDemagQuadFactor*LocalElectronOrigin.y()-1.0)*LocalElectronOrigin.y()+PsfY,
	 -(PhCathodeToSilDetDist-
          (HpdPhCathodeRInner-LocalElectronOrigin.z())));
  //normalize this vector and then transform back to global coord system.
  LocalElectronDirection = LocalElectronDirection.unit();

  const G4ThreeVector GlobalElectronDirection = theNavigator->
        			GetLocalToGlobalTransform().
				TransformAxis(LocalElectronDirection);
  
  G4double ElecKineEnergy=getHpdPhElectronKE();

  //create the electron
  G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(),
                                GlobalElectronDirection, ElecKineEnergy) ;

  aParticleChange.SetNumberOfSecondaries(1) ;
  //  aParticleChange.AddSecondary( aElectron ) ; 
     aParticleChange.AddSecondary( aElectron,GlobalElectronOrigin,true ) ; 
  

  // Kill the incident photon when it has converted to photoelectron.

   aParticleChange.ProposeLocalEnergyDeposit(PhotonEnergy);
   aParticleChange.ProposeEnergy(0.);  
   aParticleChange.ProposeTrackStatus(fStopAndKill); 
  }
  //photon is not killed if it is not converted to photoelectron
  //SE 26-09-01.
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

}

G4double PadHpdPhotoElectricEffect::getHpdQEff(G4int HpdNum, 
  G4double PhotonWLength){

  G4double hq1,hq2, wa1, wa2,aslope,aintc;
  G4double qeff=0.0;
  for (G4int ibinq=0 ; ibinq<NumQEbins-1 ; ibinq++ ){
  wa1 = HpdWabin[HpdNum][ibinq];
  wa2 = HpdWabin[HpdNum][ibinq+1];
  if( PhotonWLength >= wa1 && PhotonWLength <= wa2 ) {
   hq1 =   HpdQE[HpdNum][ibinq];
   hq2 =   HpdQE[HpdNum][ibinq+1];    
   aslope = (hq2-hq1)/(wa2-wa1);
   aintc =  hq1 - (aslope * wa1 );
   qeff= aintc + aslope * PhotonWLength ;
   return qeff;
  }

 }
  return qeff;
}




