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
// $Id: G4MultipleScatteringx.hh,v 1.5 2001-08-23 08:31:04 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//    --------- G4MultipleScatteringx physics process --------
//               by Laszlo Urban, March 2001   
//
// 07-08-01, new methods Store/Retrieve PhysicsTable 
// 23-08-01, new angle and z distribution,energy dependence reduced,
//           Store,Retrieve methods commented out temporarily, L.Urban
//            
// --------------------------------------------------------------

// class description
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MultipleScatteringx_h
#define G4MultipleScatteringx_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4PhysicsLogVector.hh"
#include "G4GPILSelection.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MultipleScatteringx : public G4VContinuousDiscreteProcess

{
 public:

   G4MultipleScatteringx(const G4String& processName="mulscat");

  ~G4MultipleScatteringx();
          
   G4bool IsApplicable ( const G4ParticleDefinition& );

   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

   void PrintInfoDefinition();
   
 //  G4bool StorePhysicsTable(G4ParticleDefinition* ,
 //			    const G4String& directory, G4bool);
       // store TransportMeanFreePath tables into an external file
       // specified by 'directory' (must exist before invokation)

 //  G4bool RetrievePhysicsTable(G4ParticleDefinition* ,
 //		               const G4String& directory, G4bool);
       // retrieve TransportMeanFreePath tables from an external file
       // specified by 'directory' 

   G4double GetContinuousStepLimit(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& currentSafety) ; 

   G4double GetMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition) ;

   G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) ;


   G4double GetLambda(G4double KineticEnergy,G4Material* material);

   void Setpalfa(G4double value) { palfa = value ; } ;
   void Setpbeta(G4double value) { pbeta = value ; } ;
   void Setpgamma(G4double value) { pgamma = value ; } ;
   void Setpq0(G4double value) { pq0 = value ; } ;
   void Setpq1(G4double value) { pq1 = value ; } ;
   void Setpc0(G4double value) { pc0 = value ; } ;
   void Setpcz(G4double value) { pcz = value ; } ;
   void Setdtrl(G4double value) { dtrl = value ; } ;

   void SetBoundary(G4bool value) { boundary = value ;} ;
   void SetFactlim(G4double val) { factlim=val;};

   G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                  G4double  previousStepSize,
                                                  G4double  currentMinimumStep,
                                                  G4double& currentSafety,
                                                  G4GPILSelection* selection);

   void SetTuning(G4double value)               {tuning = value;};
   void SetCparm (G4double value)               {cparm  = value;};
   void SetLateralDisplacementFlag(G4bool flag) {fLatDisplFlag = flag;};
   
   void SetNuclCorrPar(G4double val)            {NuclCorrPar = val;};
   void SetFactPar(G4double val)                {FactPar = val;};

 protected:

   G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition& aParticleType,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight);

   G4double TrueToGeomTransformation(const G4DynamicParticle* aParticle,
                                           G4Material* aMaterial,
                                           G4double truePathLength);
 private:

 //  hide assignment operator as  private
   G4MultipleScatteringx & operator = (const G4MultipleScatteringx &right);
   G4MultipleScatteringx ( const G4MultipleScatteringx &);

 private:        // data members ...............................................

   G4PhysicsTable* theTransportMeanFreePathTable;

   G4double fTransportMeanFreePath;

   G4double biglambda,taubig,tausmall,taulim ;

   G4double LowestKineticEnergy;
   G4double HighestKineticEnergy;
   G4int TotBin;

   const G4Electron* theElectron;
   const G4Positron* thePositron;

   G4Material* lastMaterial;
   G4double lastKineticEnergy;
   G4int materialIndex ;
  
   G4double tLast;
   G4double zLast;

   // model parameters
   G4bool boundary ;               // spec. handling near boundaries
   G4double factlim ;
   G4GPILSelection  valueGPILSelectionMSC ;

   G4double pcz,zmean ;                  // z distribution 

   G4double palfa,pbeta,pgamma,pq0,pq1,pc0 ; // theta distr.

   G4double range,T1,lambda1,cth1,z1,t1,dtrl ;

   G4double tuning;        //  param. for lambda tuning
   G4double cparm;         //          "

   // with/without lateral displacement
   G4bool fLatDisplFlag;

   // nuclear size effect correction
   G4double NuclCorrPar;
   G4double FactPar;

   G4ParticleChangeForMSC fParticleChange;
   
};

#include "G4MultipleScatteringx.icc"

#endif
 

 
