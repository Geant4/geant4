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
//
//------------- G4MultipleScattering physics process --------------------------
//               by Laszlo Urban, March 2001
//
// 07-08-01 new methods Store/Retrieve PhysicsTable
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 11-09-01 G4MultipleScatteringx put as default: G4MultipleScattering
//          Store,Retrieve methods reactived (mma)
// 13-09-01 Unused TrueToGeomTransformation method deleted,
//          class description (L.Urban)
// 19-09-01 come back to previous process name msc
// 17-04-02 NEW angle distribution + boundary algorithm modified, L.Urban
// 22-04-02 boundary algorithm modified -> important improvement in timing !!!!
//          (L.Urban)
// 24-05-02 changes in data members, L.Urban
// 30-10-02 changes in data members, L.Urban
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 05-02-03 changes in data members, L.Urban
//
//------------------------------------------------------------------------------

// class description
//
//  The class simulates the multiple scattering for any kind
//  of charged particle.
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MultipleScatteringSTD_h
#define G4MultipleScatteringSTD_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4GPILSelection.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MultipleScatteringSTD : public G4VContinuousDiscreteProcess

{
 public:    // with description

   G4MultipleScatteringSTD(const G4String& processName="msc");

  ~G4MultipleScatteringSTD();
          
   G4bool IsApplicable ( const G4ParticleDefinition& );
     // returns true for charged particles, false otherwise

   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
     // This function overloads the corresponding virtual function
     // of the base class G4VContinuousDiscreteProcess.
     // It is invoked by the G4ParticleWithCuts()::SetCut() method.
     // It prepares the table of the transport mean free paths 
     // for every material.

   void PrintInfoDefinition();
     // Print few lines of informations about the process: validity range,
     // origine ..etc..
     // Invoked by BuildPhysicsTable().

   
   G4bool StorePhysicsTable(G4ParticleDefinition* ,
 			    const G4String& directory, G4bool);
     // store TransportMeanFreePath tables into an external file
     // specified by 'directory' (must exist before invokation)

   G4bool RetrievePhysicsTable(G4ParticleDefinition* ,
 		               const G4String& directory, G4bool);
     // retrieve TransportMeanFreePath tables from an external file
     // specified by 'directory'


   G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                  G4double  previousStepSize,
                                                  G4double  currentMinimumStep,
                                                  G4double& currentSafety,
                                                  G4GPILSelection* selection);
     // The function overloads the corresponding function of the base
     // class.It limits the step near to boundaries only
     // and invokes the method GetContinuousStepLimit at every step.

   G4double GetContinuousStepLimit(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& currentSafety);
     // It performs the true step length --> geometrical step length
     // transformation. It is invoked by the
     // AlongStepGetPhysicalInteractionLength method.

   G4double GetMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition);
     // It sets the force condition to true only
     // in order to have the PostStepDoIt called at every step.
     // This function overloads a virtual function of the base class.
     // It is invoked by the ProcessManager of the Particle.


   G4double GetTransportMeanFreePath(
                          G4double KineticEnergy,G4Material* material);
     // Just a utility method to get the values of the transport
     //  mean free path . (It is not used inside the class.)

   G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);
     // The geometrical step length --> true path length transformation
     // is performed here (the inverse of the transformation done
     // by GetContinuousStepLimit).

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep);
     // It computes the final state of the particle: samples the
     // ScatteringSTD angle and computes the lateral displacement.
     // The final state is returned as a ParticleChange object.
     // This function overloads a virtual function of the base class.
     // It is invoked by the ProcessManager of the Particle.

   void Setsamplez(G4bool value)               {samplez = value;};
     // geom. step length distribution should be sampled or not

   void Setdtrl(G4double value)                 {dtrl = value;};
     // to reduce the energy/step dependence

   void SetBoundary(G4bool value)      {boundary = value;};
   void Setfacxsi(G4double value)      {facxsi = value;
                                    G4cout << " facxsi=" << facxsi << G4endl;};
   void SetFacrange(G4double val)  {facrange=val;
                                    nsmallstep = G4int(log((cf+facrange-1.)/
                                                 facrange)/log(cf))+1 ;
                                    G4cout << " fr=" << facrange
                                           << "  nsmall=" << nsmallstep << G4endl ;};
      // Steplimit after boundary crossing = facrange*range
      // estimated nb of steps at boundary nsmallstep = 1/facrange

   void SetLateralDisplacementFlag(G4bool flag) {fLatDisplFlag = flag;};
     // lateral displacement to be/not to be computed

   void SetNuclCorrPar(G4double val)            {NuclCorrPar = val;};
   void SetFactPar(G4double val)                {FactPar = val;};
     // corrs to transport cross section for high energy

 protected:    // with description

   virtual G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition& aParticleType,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight);
     // It computes the transport cross section.
     // The transport mean free path is 1/(transport cross section).

 private:

 //  hide assignment operator as  private
   G4MultipleScatteringSTD & operator = (const G4MultipleScatteringSTD &right);
   G4MultipleScatteringSTD ( const G4MultipleScatteringSTD &);

 private:        // data members

   G4PhysicsTable* theTransportMeanFreePathTable;

   G4double taubig,tausmall,taulim;

   G4double LowestKineticEnergy;
   G4double HighestKineticEnergy;
   G4int    TotBin;

   G4int    materialIndex;

   G4double tLast;
   G4double zLast;

   // model parameters
   G4bool   boundary;                         // spec. handling near boundaries
   G4double facrange,tlimit,tlimitmin,cf;
   G4int    stepno, stepnolastmsc, nsmallstep;
   G4double laststep ;
   G4GPILSelection  valueGPILSelectionMSC;

   G4double zmean;                            // z(geom.step length)
   G4bool samplez ;                           //  distribution

   G4double range,T0,T1,lambda0,lambda1,      // used to reduce the energy
            Tlow,alam,blam,dtrl,lambdam,      // (or step length) dependence
            clam,zm,cthm;

   // with/without lateral displacement
   G4bool fLatDisplFlag;

   // nuclear size effect correction
   G4double NuclCorrPar;
   G4double FactPar;

   G4ParticleChangeForMSC fParticleChange; 

   G4double alfa1,alfa2,alfa3,b,xsi,c0,facxsi ;  // angle distr. parameters
                                                 // facxsi : some tuning 
                                                 // possibility in the tail 

};

#include "G4MultipleScatteringSTD.icc"

#endif
 

 
