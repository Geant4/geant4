// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SynchrotronRadiation.hh,v 1.1 1999-01-07 16:11:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      
//      History: 
//      21-5-98  1 version , V. Grichine
//                   
// 
//
// ------------------------------------------------------------

#ifndef G4SynchrotronRadiation_h
#define G4SynchrotronRadiation_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"

#include "G4Track.hh"
#include "G4Step.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4SynchrotronRadiation : public G4VDiscreteProcess 
{ 
  public:
 
     G4SynchrotronRadiation(const G4String& processName = 
                                                   "SynchrotronRadiation");
 
    ~G4SynchrotronRadiation();

  private:

     G4SynchrotronRadiation & operator=(const G4SynchrotronRadiation &right);

     G4SynchrotronRadiation(const G4SynchrotronRadiation&);

  public:  /////////////////    Post Step functions  //////////////////////////  

     G4double GetMeanFreePath( const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt( const G4Track& track,         
                                      const G4Step& Step    ) ;

     G4double GetPhotonEnergy( const G4Track& trackData,
                               const G4Step&  stepData      ) ;
                 

     G4bool IsApplicable(const G4ParticleDefinition&);

  static G4double GetLambdaConst()  { return fLambdaConst ; } ;

  static G4double GetEnergyConst()  { return fEnergyConst ; } ;

  protected:

  private:

     static const G4double fLambdaConst ;

     static const G4double fEnergyConst ;

     static const G4double fIntegralProbabilityOfSR[200] ;


     const G4double 
     LowestKineticEnergy;   // low  energy limit of the cross-section formula

     const G4double 
     HighestKineticEnergy;  // high energy limit of the cross-section formula
 
     G4int TotBin;          // number of bins in the tables

     G4double CutInRange;

     const G4Gamma*    theGamma; 
     const G4Electron* theElectron;
     const G4Positron* thePositron;

     const G4double* GammaCutInKineticEnergy;
     const G4double* ElectronCutInKineticEnergy;
     const G4double* PositronCutInKineticEnergy;
     const G4double* ParticleCutInKineticEnergy;


     G4double GammaCutInKineticEnergyNow;
     G4double ElectronCutInKineticEnergyNow;
     G4double PositronCutInKineticEnergyNow;
     G4double ParticleCutInKineticEnergyNow;
};

//////////////////////////  INLINE METHODS  /////////////////////////////
//
// gives the MeanFreePath in GEANT4 internal units
//

inline G4double 
G4SynchrotronRadiation::GetMeanFreePath( const G4Track& trackData,
                                               G4double previousStepSize,
                                               G4ForceCondition* condition)
{
   const G4DynamicParticle* aDynamicParticle;
   G4Material* aMaterial;
   G4double MeanFreePath;
   G4bool isOutRange ;
 
   *condition = NotForced ;

   aDynamicParticle = trackData.GetDynamicParticle();
   aMaterial = trackData.GetMaterial();

   G4double gamma = aDynamicParticle->GetTotalEnergy()/
                   (aDynamicParticle->GetMass()              ) ;
 
   G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();

   if (KineticEnergy <  LowestKineticEnergy || gamma<1.0e3)
   {
     MeanFreePath = DBL_MAX ;
   }
   else
   {
     G4TransportationManager* transportMgr;
     G4FieldManager* globalFieldMgr;

     transportMgr = G4TransportationManager::GetTransportationManager() ;
  
     globalFieldMgr = transportMgr->GetFieldManager() ;

     G4bool FieldExists = globalFieldMgr->DoesFieldExist() ;
     G4ThreeVector  FieldValue;
     G4Field*   pField = 0 ;
     if (FieldExists)
     {  
       pField = globalFieldMgr->GetDetectorField() ;
       G4ThreeVector  globPosition = trackData.GetPosition() ;
       G4double  globPosVec[3], FieldValueVec[3] ;
       globPosVec[0] = globPosition.x() ;
       globPosVec[1] = globPosition.y() ;
       globPosVec[2] = globPosition.z() ;

       pField->GetFieldValue( globPosVec, FieldValueVec ) ;
       FieldValue = G4ThreeVector( FieldValueVec[0], 
                                   FieldValueVec[1], 
                                   FieldValueVec[2]   ) ;

       
       
       G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection(); 
       G4ThreeVector unitMcrossB = FieldValue.cross(unitMomentum) ;
       G4double perpB = unitMcrossB.mag() ;
       G4double beta = aDynamicParticle->GetTotalMomentum()/
                      (aDynamicParticle->GetTotalEnergy()    ) ;
       if(perpB > 0.0)
       {
         MeanFreePath = fLambdaConst*beta/perpB ;
       }       
       else
       {
         MeanFreePath = DBL_MAX ;
       }
     }
     else
     {
       MeanFreePath = DBL_MAX ;
     }
   }
   return MeanFreePath; 
} 



inline G4bool 
G4SynchrotronRadiation::IsApplicable( const G4ParticleDefinition& particle )
{
   return(   (&particle == (const G4ParticleDefinition *)theElectron)
           ||(&particle == (const G4ParticleDefinition *)thePositron)
         ) ;
}
  
#endif  // end of G4SynchrotronRadiation.hh
 
