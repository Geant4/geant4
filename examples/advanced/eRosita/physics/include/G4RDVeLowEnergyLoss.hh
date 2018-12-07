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
//
//
// 3.4.2000 Veronique Lefebure:
//          Move utils/include/G4VEnergyLoss.hh to 
//               lowenergy/include/G4RDVeLowEnergyLoss.hh
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
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
//  All the EnergyLoss classes are inherited from G4RDVeLowEnergyLoss
//  class.
//
//  -----------------------------------------------------------
//  created  on 28 January 2000  by L. Urban
//  -----------------------------------------------------------
//
//  Modifications:
// 20/09/00 V.Ivanchenko update fluctuations
// 23/11/01 V.Ivanchenko Move static member-functions from header to source
// 22/01/03 V.Ivanchenko Cut per range
//
// Class description:
// Abstract class for Low Energy Electromagnetic electron energy loss
// Further documentation available from http://www.ge.infn.it/geant4/lowE

//  -----------------------------------------------------------

#ifndef G4RDVeLowEnergyLoss_h
#define G4RDVeLowEnergyLoss_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Electron.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4MaterialCutsCouple.hh"

class G4RDVeLowEnergyLoss : public G4VContinuousDiscreteProcess
{
  public:

      G4RDVeLowEnergyLoss(const G4String& ,
				   G4ProcessType   aType = fNotDefined );
      G4RDVeLowEnergyLoss(G4RDVeLowEnergyLoss &);

      virtual ~G4RDVeLowEnergyLoss();

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
                              const G4MaterialCutsCouple* couple,
                              G4double	MeanLoss,
                              G4double  step);


   private:

  // hide default constructor and assignment operator as private
      G4RDVeLowEnergyLoss();
      G4RDVeLowEnergyLoss & operator=(const G4RDVeLowEnergyLoss &right);

  protected:

    // data members to speed up the fluctuation calculation
    const G4Material* lastMaterial;
    G4int imat;
    G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
    G4double e1LogFluct,e2LogFluct,ipotLogFluct;

    const G4int nmaxCont1,nmaxCont2 ;

  // static part of the class

  public:  // With description

    static void SetRndmStep     (G4bool   value);
    // use / do not use randomisation in energy loss steplimit
    // ( default = no randomisation)

    static void SetEnlossFluc   (G4bool   value);
    // compute energy loss with/without fluctuation
    // ( default : with fluctuation)

    static void SetStepFunction (G4double c1, G4double c2);
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

    static void BuildRangeVectorNew(const G4PhysicsTable*,G4int,
                                          G4int,G4PhysicsLogVector*);

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



