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
// $Id: G4Step.hh 93028 2015-09-30 16:09:00Z gcosmo $
//
//
//---------------------------------------------------------------
//
// G4Step.hh
//
// Class Description:
//   This class represents the Step of a particle tracked.
//   It includes information of 
//     1) List of Step points which compose the Step,
//     2) static information of particle which generated the 
//        Step, 
//     3) trackID and parent particle ID of the Step,
//     4) termination condition of the Step,
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------
//   Modified for the new G4ParticleChange          12 Mar. 1998  H.Kurahige
//   Correct treatment of touchable in G4Step::UpdateTrack
//                                                  12 May. 1998 H.Kurashige
// ---------------------------------------------------------------
//   Separate implementation of inline functions inti G4Step.icc
//   Add updating mass/charge                        6 Oct. 1999 H.Kurashige
//   add nonIonizingEnergyLoss                      26 Mar. 2007 H.Kurashige 
//
//   Repository test - Dennis Wright
//
#ifndef G4Step_h
#define G4Step_h 1

#include <stdlib.h>                 // Include from 'system'
#include <cmath>                    // Include from 'system'
#include "G4ios.hh"                 // Include from 'system'
#include <iomanip>                  // Include from 'system'
#include "globals.hh"               // Include from 'global'
#include "G4ThreeVector.hh"         // Include from 'global'
#include "G4VPhysicalVolume.hh"     // Include from 'geometry'
#include "G4StepPoint.hh"           // Include from 'track'
#include "G4StepStatus.hh"          // Include from 'track'
class G4Polyline;                   // Forward declaration.
class G4Track;                      // Forward declaration.
#include "G4TrackVector.hh"           // Include from 'tracking'

////////////
class G4Step
////////////
{

//--------
   public:

// Constructor/Destrcutor
   G4Step();
   ~G4Step();

// Copy Counstructor and assignment operator
   G4Step(const G4Step& );
   G4Step & operator=(const G4Step &);   

//--------
   public: // WIth description

// Get/Set functions 
   // currnet track
   G4Track* GetTrack() const;
   void SetTrack(G4Track* value);

   // step points 
   G4StepPoint* GetPreStepPoint() const;
   void SetPreStepPoint(G4StepPoint* value);

   G4StepPoint* GetPostStepPoint() const;
   void SetPostStepPoint(G4StepPoint* value);

   // step length
   G4double GetStepLength() const;
   void SetStepLength(G4double value);
    // Before the end of the AlongStepDoIt loop,StepLength keeps
    // the initial value which is determined by the shortest geometrical Step
    // proposed by a physics process. After finishing the AlongStepDoIt,
    // it will be set equal to 'StepLength' in G4Step. 

   // total energy deposit 
   G4double GetTotalEnergyDeposit() const;
   void SetTotalEnergyDeposit(G4double value);

   // total non-ionizing energy deposit 
   G4double GetNonIonizingEnergyDeposit() const;
   void SetNonIonizingEnergyDeposit(G4double value);

   // cotrole flag for stepping
   G4SteppingControl GetControlFlag() const;
   void SetControlFlag(G4SteppingControl StepControlFlag);

    // manipulation of total energy deposit 
   void AddTotalEnergyDeposit(G4double value);
   void ResetTotalEnergyDeposit();

   // manipulation of non-ionizng energy deposit 
   void AddNonIonizingEnergyDeposit(G4double value);
   void ResetNonIonizingEnergyDeposit();


  // Get/Set/Clear flag for initial/last step
   // NOTE:  following flags are not used 
   //        will be ready in later release
   G4bool IsFirstStepInVolume() const;
   G4bool IsLastStepInVolume() const;

   void SetFirstStepFlag();
   void ClearFirstStepFlag();
   void SetLastStepFlag();
   void ClearLastStepFlag();

  // difference of position, time, momentum and energy
   G4ThreeVector GetDeltaPosition() const;
   G4double GetDeltaTime() const;

  // These methods will be deleted 
  // NOTE: use  GetTotalEnergyDeposit() to obtain 
  //       energy loss in the material 
  // 
   G4ThreeVector GetDeltaMomentum() const;
   G4double GetDeltaEnergy() const;


// Other member functions
   void InitializeStep( G4Track* aValue );
   // initiaize contents of G4Step

   void UpdateTrack( );
   // update track by using G4Step information

   void CopyPostToPreStepPoint( );
   // copy PostStepPoint to PreStepPoint 
  
   G4Polyline* CreatePolyline () const;
   // for visualization

//-----------
   protected:
//-----------

// Member data
   G4double fTotalEnergyDeposit;
     // Accummulated total energy desposit in the current Step

   G4double fNonIonizingEnergyDeposit;
     // Accummulated non-ionizing energy desposit in the current Step

//---------
   private:
//---------

// Member data
   G4StepPoint* fpPreStepPoint;
   G4StepPoint* fpPostStepPoint;
   G4double fStepLength;
     // Step length which may be updated at each invocation of 
     // AlongStepDoIt and PostStepDoIt
   G4Track* fpTrack;
     //
   G4SteppingControl fpSteppingControlFlag;     
    // A flag to control SteppingManager behavier from process

  // flag for initial/last step
   G4bool fFirstStepInVolume;
   G4bool fLastStepInVolume;

// Secondary buckets
public:
  // secodaries in the current step
   G4int GetNumberOfSecondariesInCurrentStep() const;

   const std::vector<const G4Track*>* GetSecondaryInCurrentStep() const; 

   // NOTE: Secondary bucket of the Step contains  
   //       all secondaries during tracking the current track 
   //       (i.e. NOT secondaries produced in the current step)
   // all following methods give same object (i.e. G4TrackVector  )
   // but 2nd one will create bucket in addition  
   const G4TrackVector* GetSecondary() const ;
   G4TrackVector* GetfSecondary();
   G4TrackVector* NewSecondaryVector();

   // just delete secondary bucket
   //  NOTE: G4Track objects inside the bucket are not deleted 
   void DeleteSecondaryVector();

   // Add secondary tracks to the bucket 
   void SetSecondary( G4TrackVector* value);

private: 
   // Secondaty bucket implemented by using  std::vector of G4Track*   
   G4TrackVector* fSecondary;

   // number of secondaries which have been created by the last step
   G4int  nSecondaryByLastStep;

   typedef const G4Track* CT;
   std::vector<CT>* secondaryInCurrentStep;

  // Prototyping implementation of smooth representation of curved
  // trajectories. (jacek 30/10/2002)
public:
  // Auxiliary points are ThreeVectors for now; change to
  // G4VAuxiliaryPoints or some such (jacek 30/10/2002)
  void SetPointerToVectorOfAuxiliaryPoints( std::vector<G4ThreeVector>* theNewVectorPointer ) {
    fpVectorOfAuxiliaryPointsPointer = theNewVectorPointer;
  }
  std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const {
    return fpVectorOfAuxiliaryPointsPointer;
  }
private:
  // Explicity including the word "Pointer" in the name as I keep
  // forgetting the * (jacek 30/10/2002)
  std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer;

};

#include "G4Step.icc"


#endif
