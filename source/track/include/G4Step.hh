// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Step.hh,v 1.5 1999-10-06 01:21:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4Step.hh
//
// Description:
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
//   Add updating mass/charge  
//                                                  6 Oct. 1999 H.Kurashige
//
#ifndef G4Step_h
#define G4Step_h 1

#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"               // Include from 'system'
#include <iomanip.h>                // Include from 'system'
#include "globals.hh"               // Include from 'global'
#include "G4ThreeVector.hh"         // Include from 'global'
#include "G4VPhysicalVolume.hh"     // Include from 'geometry'
#include "G4StepPoint.hh"           // Include from 'track'
#include "G4StepStatus.hh"          // Include from 'track'
class G4Polyline;                   // Forward declaration.
class G4Track;                      // Forward declaration.

////////////
class G4Step
////////////
{

//--------
   public:
//--------

// Constructor/Destrcutor
   G4Step();
   ~G4Step();

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

  // total energy deposit 
   G4double GetTotalEnergyDeposit() const;
   void SetTotalEnergyDeposit(G4double value);

   // cotrole flag for stepping
   G4SteppingControl GetControlFlag() const;
   void SetControlFlag(G4SteppingControl StepControlFlag);

   // difference of position, time, momentum and energy
   G4ThreeVector GetDeltaPosition() const;
   G4double GetDeltaTime() const;

   G4ThreeVector GetDeltaMomentum() const;
   G4double GetDeltaEnergy() const;

   // manipulation of total energy deposit 
   void AddTotalEnergyDeposit(G4double value);
   void ResetTotalEnergyDeposit();


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
};

#include "G4Step.icc"


#endif










