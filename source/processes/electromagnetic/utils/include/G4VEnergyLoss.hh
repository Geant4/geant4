// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VEnergyLoss.hh,v 1.3 2000-02-17 09:08:42 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// Class Description 
//
//  General service class for the energy loss classes
//  
//  It contains code needed to compute the range tables,
//  time tables, the inverse range tables and some auxiliary
//  tables.
//  The energy loss fluctuation code is here,too.
//
//  All the EnergyLoss classes are inherited from G4VEnergyLoss
//  class.
//
//  -----------------------------------------------------------
//  created  on 28 January 2000  by L. Urban               
//  -----------------------------------------------------------

#ifndef G4VEnergyLoss_h
#define G4VEnergyLoss_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Electron.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"


class G4VEnergyLoss : public G4VContinuousDiscreteProcess 
{
  public:     

      G4VEnergyLoss(const G4String& ,
				   G4ProcessType   aType = fNotDefined );
      G4VEnergyLoss(G4VEnergyLoss &);

      virtual ~G4VEnergyLoss();

      virtual G4double GetContinuousStepLimit(const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety) = 0 ;

      virtual G4VParticleChange* AlongStepDoIt(const G4Track& track,
                                     const G4Step& Step) = 0 ;

      virtual G4double GetMeanFreePath(const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) = 0;

      virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                            const G4Step& Step) = 0;



  protected:// with description

    // code for the energy loss fluctuation

    G4double GetLossWithFluct(const G4DynamicParticle* aParticle,
                              G4Material* aMaterial,
                              G4double	 threshold);


   private:

  // hide default constructor and assignment operator as private 
      G4VEnergyLoss();
      G4VEnergyLoss & operator=(const G4VEnergyLoss &right);

  protected:

    // data members to speed up the fluctuation calculation
    G4Material* lastMaterial;
    G4int imat;
    G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
    G4double e1LogFluct,e2LogFluct,ipotLogFluct;

    const G4double MaxExcitationNumber ;
    const G4int nmaxCont1,nmaxCont2 ;

  // static part of the class 

  public:  // With description

    static void SetRndmStep     (G4bool   value) {rndmStepFlag   = value;}
    // use / do not use randomisation in energy loss steplimit
    // ( default = no randomisation)

    static void SetEnlossFluc   (G4bool   value) {EnlossFlucFlag = value;}
    // compute energy loss with/without fluctuation
    // ( default : with fluctuation)

    static void SetStepFunction (G4double c1, G4double c2)
                               {dRoverRange = c1; finalRange = c2;
                                c1lim=dRoverRange ;
                                c2lim=2.*(1-dRoverRange)*finalRange;
                                c3lim=-(1.-dRoverRange)*finalRange*finalRange;
                               }
    // sets values for data members used to compute the step limit:
    //   dRoverRange : max. relative range change in one step,
    //   finalRange  : if range <= finalRange --> last step for the particle.


  protected: // With description

    // Build range table starting from the DEDXtable
    static G4PhysicsTable*
     BuildRangeTable(G4PhysicsTable* theDEDXTable,
                     G4PhysicsTable* theRangeTable,
                     G4double Tmin,G4double Tmax,G4int nbin);

    // Build time tables starting from the DEDXtable
    static G4PhysicsTable*
     BuildLabTimeTable(G4PhysicsTable* theDEDXTable,
                       G4PhysicsTable* theLabTimeTable,
                       G4double Tmin,G4double Tmax,G4int nbin);

    static G4PhysicsTable*
     BuildProperTimeTable(G4PhysicsTable* theDEDXTable,
                       G4PhysicsTable* ProperTimeTable,
                       G4double Tmin,G4double Tmax,G4int nbin);

    // Build tables of coefficients needed for inverting the range table 
    static G4PhysicsTable*
     BuildRangeCoeffATable(G4PhysicsTable* theRangeTable,
                           G4PhysicsTable* theCoeffATable,
                           G4double Tmin,G4double Tmax,G4int nbin);
    static G4PhysicsTable*
     BuildRangeCoeffBTable(G4PhysicsTable* theRangeTable,
                           G4PhysicsTable* theCoeffBTable,
                           G4double Tmin,G4double Tmax,G4int nbin);
    static G4PhysicsTable*
     BuildRangeCoeffCTable(G4PhysicsTable* theRangeTable,
                           G4PhysicsTable* theCoeffCTable,
                           G4double Tmin,G4double Tmax,G4int nbin);

    // Invert range table
    static G4PhysicsTable*
     BuildInverseRangeTable(G4PhysicsTable* theRangeTable,
                            G4PhysicsTable* theRangeCoeffATable,
                            G4PhysicsTable* theRangeCoeffBTable,
                            G4PhysicsTable* theRangeCoeffCTable,
                            G4PhysicsTable* theInverseRangeTable,
                            G4double Tmin,G4double Tmax,G4int nbin);

  private:

    static void BuildRangeVector(G4PhysicsTable* theDEDXTable,
                        G4double Tmin,G4double Tmax,G4int nbin,
                        G4int materialIndex,G4PhysicsLogVector* rangeVector);

    static G4double RangeIntLin(G4PhysicsVector* physicsVector
                                                        ,G4int nbin);

    static G4double RangeIntLog(G4PhysicsVector* physicsVector
                                                        ,G4int nbin);

    static void BuildLabTimeVector(G4PhysicsTable* theDEDXTable,
                        G4double Tmin,G4double Tmax,G4int nbin,
                        G4int materialIndex,G4PhysicsLogVector* rangeVector);

    static void BuildProperTimeVector(G4PhysicsTable* theDEDXTable,
                        G4double Tmin,G4double Tmax,G4int nbin,
                        G4int materialIndex,G4PhysicsLogVector* rangeVector);

    static G4double LabTimeIntLog(G4PhysicsVector* physicsVector
                                                        ,G4int nbin);

    static G4double ProperTimeIntLog(G4PhysicsVector* physicsVector,
                                                         G4int nbin);

    static void InvertRangeVector(G4PhysicsTable* theRangeTable,
                                  G4PhysicsTable* theRangeCoeffATable,
                                  G4PhysicsTable* theRangeCoeffBTable,
                                  G4PhysicsTable* theRangeCoeffCTable,
                                  G4double Tmin,G4double Tmax,G4int nbin,
                       G4int materialIndex,G4PhysicsLogVector* rangeVector);


  // data members
  protected:

   // variables for the integration routines
   static G4double ParticleMass,taulow,tauhigh,ltaulow,ltauhigh;


   static G4double dRoverRange;     // dRoverRange is the maximum allowed
                                     // deltarange/range in one Step
   static G4double finalRange;      // final step before stopping
   static G4double c1lim,c2lim,c3lim ; // coeffs for computing steplimit

   static G4bool   rndmStepFlag;    // control the randomization of the step
   static G4bool   EnlossFlucFlag;  // control the energy loss fluctuation


};

#endif



