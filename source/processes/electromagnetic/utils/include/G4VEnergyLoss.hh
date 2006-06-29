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
// $Id: G4VEnergyLoss.hh,v 1.18 2006-06-29 19:54:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
// 26.10.01 static inline functions moved to .cc file (mma)
// 08.11.01 some static methods,data members are not static L.Urban
// 15.01.03 Migrade to cut per region (V.Ivanchenko)
// ------------------------------------------------------------
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
#include "G4MaterialCutsCouple.hh"

class G4EnergyLossMessenger;

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
                              const G4MaterialCutsCouple* couple,
                              G4double ChargeSquare,
                              G4double	 MeanLoss,
                              G4double step);

    // Build range table starting from the DEDXtable
    G4PhysicsTable*
    BuildRangeTable(G4PhysicsTable* theDEDXTable,
                    G4PhysicsTable* theRangeTable,
                    G4double Tmin,G4double Tmax,G4int nbin);

    // Build time tables starting from the DEDXtable
    G4PhysicsTable*
    BuildLabTimeTable(G4PhysicsTable* theDEDXTable,
                      G4PhysicsTable* theLabTimeTable,
                      G4double Tmin,G4double Tmax,G4int nbin);

    G4PhysicsTable*
    BuildProperTimeTable(G4PhysicsTable* theDEDXTable,
                      G4PhysicsTable* ProperTimeTable,
                      G4double Tmin,G4double Tmax,G4int nbin);

    // Build tables of coefficients needed for inverting the range table
    G4PhysicsTable*
    BuildRangeCoeffATable(G4PhysicsTable* theRangeTable,
                          G4PhysicsTable* theCoeffATable,
                          G4double Tmin,G4double Tmax,G4int nbin);
    G4PhysicsTable*
    BuildRangeCoeffBTable(G4PhysicsTable* theRangeTable,
                          G4PhysicsTable* theCoeffBTable,
                          G4double Tmin,G4double Tmax,G4int nbin);
    G4PhysicsTable*
    BuildRangeCoeffCTable(G4PhysicsTable* theRangeTable,
                          G4PhysicsTable* theCoeffCTable,
                          G4double Tmin,G4double Tmax,G4int nbin);

    // Invert range table
    G4PhysicsTable*
    BuildInverseRangeTable(G4PhysicsTable* theRangeTable,
                           G4PhysicsTable* theRangeCoeffATable,
                           G4PhysicsTable* theRangeCoeffBTable,
                           G4PhysicsTable* theRangeCoeffCTable,
                           G4PhysicsTable* theInverseRangeTable,
                           G4double Tmin,G4double Tmax,G4int nbin);

   private:

  // hide default constructor and assignment operator as private
      G4VEnergyLoss();
      G4VEnergyLoss & operator=(const G4VEnergyLoss &right);

      void BuildRangeVector(G4PhysicsTable* theDEDXTable,
                     G4double Tmin,G4double Tmax,G4int nbin,
                     G4int materialIndex,G4PhysicsLogVector* rangeVector);

      G4double RangeIntLin(G4PhysicsVector* physicsVector,G4int nbin);

      G4double RangeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

      void BuildLabTimeVector(G4PhysicsTable* theDEDXTable,
                     G4double Tmin,G4double Tmax,G4int nbin,
                     G4int materialIndex,G4PhysicsLogVector* rangeVector);

      void BuildProperTimeVector(G4PhysicsTable* theDEDXTable,
                     G4double Tmin,G4double Tmax,G4int nbin,
                     G4int materialIndex,G4PhysicsLogVector* rangeVector);

      G4double LabTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

      G4double ProperTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

      void InvertRangeVector(G4PhysicsTable* theRangeTable,
                             G4PhysicsTable* theRangeCoeffATable,
                             G4PhysicsTable* theRangeCoeffBTable,
                             G4PhysicsTable* theRangeCoeffCTable,
                             G4double Tmin,G4double Tmax,G4int nbin,
                       G4int materialIndex,G4PhysicsLogVector* rangeVector);


  protected:

    G4double ParticleMass;

  private:

    // data members to speed up the fluctuation calculation
    const G4Material* lastMaterial;
    G4int imat;
    G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
    G4double e1LogFluct,e2LogFluct,ipotLogFluct;

    const G4int nmaxCont1,nmaxCont2 ;

    // for some integration routines
    G4double taulow,tauhigh,ltaulow,ltauhigh;

  // static part of the class

  public:  // With description

    static void SetRndmStep(G4bool value);
    // use / do not use randomisation in energy loss steplimit
    // ( default = no randomisation)

    static void SetEnlossFluc(G4bool value);
    // compute energy loss with/without fluctuation
    // ( default : with fluctuation)

    static void SetSubSec(G4bool value);
    // switch on/off the generation of the subcutoff secondaries
    // ( default = no subcutoff secondary generation )

    static void SetMinDeltaCutInRange(G4double value);
    // sets minimal cut value for the subcutoff secondaries
    // (i.e. the kinetic energy of these secondaries can not be
    //	smaller than the energy corresponds to MinDeltaCutInRange).


    static void SetStepFunction (G4double c1, G4double c2);
    // sets values for data members used to compute the step limit:
    //   dRoverRange : max. relative range change in one step,
    //   finalRange  : if range <= finalRange --> last step for the particle.


  protected: // With description

     static G4bool EqualCutVectors( G4double* vec1, G4double* vec2 );
     static G4double* CopyCutVectors( G4double* dest, G4double* source );
     G4bool CutsWhereModified();

  // data members
  protected:

   static G4double dRoverRange;     // dRoverRange is the maximum allowed
                                     // deltarange/range in one Step
   static G4double finalRange;      // final step before stopping
   static G4double finalRangeRequested; //from UI command
   static G4double c1lim,c2lim,c3lim ; // coeffs for computing steplimit

   static G4bool   rndmStepFlag;    // control the randomization of the step
   static G4bool   EnlossFlucFlag;  // control the energy loss fluctuation
   static G4bool       subSecFlag;  // control the generation of subcutoff delta

   static G4double MinDeltaCutInRange; // minimum cut for delta rays
   static G4double* MinDeltaEnergy ;
   static G4bool* LowerLimitForced ;

   static G4bool setMinDeltaCutInRange ;

   static G4EnergyLossMessenger* ELossMessenger;
};

#endif



