// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Step.hh,v 1.3 1999-04-13 09:43:26 kurasige Exp $
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
   G4StepPoint* GetPreStepPoint() const;
   void SetPreStepPoint(G4StepPoint* value);

   G4StepPoint* GetPostStepPoint() const;
   void SetPostStepPoint(G4StepPoint* value);

   G4double GetStepLength() const;
   void SetStepLength(G4double value);

   G4Track* GetTrack() const;
   void SetTrack(G4Track* value);

   G4ThreeVector GetDeltaPosition() const;
   G4double GetDeltaTime() const;

   G4ThreeVector GetDeltaMomentum() const;

   G4double GetDeltaEnergy() const;

   G4double GetTotalEnergyDeposit() const;

   void SetTotalEnergyDeposit(G4double value);
   void AddTotalEnergyDeposit(G4double value);
   void ResetTotalEnergyDeposit();

   G4SteppingControl GetControlFlag() const;
   void SetControlFlag(G4SteppingControl StepControlFlag);

// Other member functions
   void InitializeStep( G4Track* aValue );

   void UpdateTrack( );

   void CopyPostToPreStepPoint( );

//TS moved to SteppingVerbose
//   void ShowStep() const ;
     // Print all information of the Step to stdout

   G4Polyline* CreatePolyline () const;


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


//-----------------------------------------------------------------
//  In-line definitions
//-----------------------------------------------------------------

// Get/Set functions 
   inline G4StepPoint* G4Step::GetPreStepPoint() const 
   { return fpPreStepPoint; }
   inline void G4Step::SetPreStepPoint(G4StepPoint* value)
   { fpPreStepPoint = value; }

   inline G4StepPoint* G4Step::GetPostStepPoint() const
   { return fpPostStepPoint; }
   inline void G4Step::SetPostStepPoint(G4StepPoint* value)
   { fpPostStepPoint = value; }

   inline G4double G4Step::GetStepLength() const
   { return fStepLength; }
   inline void G4Step::SetStepLength(G4double value)
   { fStepLength = value; }

   inline G4ThreeVector G4Step::GetDeltaPosition() const
   { return fpPostStepPoint->GetPosition()
            - fpPreStepPoint->GetPosition(); }

   inline G4double G4Step::GetDeltaTime() const
   { return fpPostStepPoint->GetLocalTime()
            - fpPreStepPoint->GetLocalTime(); }

   inline G4ThreeVector G4Step::GetDeltaMomentum() const
   { return fpPostStepPoint->GetMomentum()
            - fpPreStepPoint->GetMomentum(); }

   inline G4double G4Step::GetDeltaEnergy() const
   { return fpPostStepPoint->GetKineticEnergy()
            - fpPreStepPoint->GetKineticEnergy(); }

   inline G4double G4Step::GetTotalEnergyDeposit() const
   { return fTotalEnergyDeposit; }

   inline void G4Step::SetTotalEnergyDeposit(G4double value)
   { fTotalEnergyDeposit = value;   }

   inline void G4Step::AddTotalEnergyDeposit(G4double value)
   { fTotalEnergyDeposit += value;   }

   inline void G4Step::ResetTotalEnergyDeposit()
   { fTotalEnergyDeposit = 0.; }

   inline void G4Step::SetControlFlag(G4SteppingControl value)
   {
       fpSteppingControlFlag = value;     
   }

   inline G4SteppingControl G4Step::GetControlFlag() const
   {
       return fpSteppingControlFlag;     
   }

   inline void G4Step::CopyPostToPreStepPoint( )
   { 
   // Default equal operator is used for copy
      *(fpPreStepPoint) = *(fpPostStepPoint);
   }


//-------------------------------------------------------------
// To implement bi-directional association between G4Step and
// and G4Track, a combined usage of 'forward declaration' and
// 'include' is necessary.
//-------------------------------------------------------------
#include "G4Track.hh"

   inline G4Track* G4Step::GetTrack() const
   { return fpTrack; }
   inline void G4Step::SetTrack(G4Track* value)
   { fpTrack = value; }


// Other member functions
   inline void G4Step::InitializeStep( G4Track* aValue )
   {
     // Initialize G4Step attributes
     fStepLength = 0.;
     fTotalEnergyDeposit = 0.;
     fpTrack = aValue;
     fpTrack->SetStepLength(0.);

     // Initialize G4StepPoint attributes.
     // To avoid the circular dependency between G4Track, G4Step
     // and G4StepPoint, G4Step has to manage the copy actions.
     fpPreStepPoint->SetPosition(fpTrack->GetPosition());
     fpPreStepPoint->SetGlobalTime(fpTrack->GetGlobalTime());
     fpPreStepPoint->SetLocalTime(fpTrack->GetLocalTime());
     fpPreStepPoint->SetProperTime(fpTrack->GetProperTime());
     fpPreStepPoint->SetMomentumDirection(fpTrack->GetMomentumDirection());
     fpPreStepPoint->SetKineticEnergy(fpTrack->GetKineticEnergy());
     fpPreStepPoint->SetTouchable(fpTrack->GetTouchable());
     fpPreStepPoint->SetPolarization(fpTrack->GetPolarization());
     fpPreStepPoint->SetSafety(0.);
     fpPreStepPoint->SetStepStatus(fUndefined);
     fpPreStepPoint->SetProcessDefinedStep(0);
     fpPreStepPoint->SetMass(fpTrack->GetDynamicParticle()->GetMass());    
     fpPreStepPoint->SetWeight(fpTrack->GetWeight());

     (*fpPostStepPoint) = (*fpPreStepPoint);
   }

   inline void G4Step::UpdateTrack( )
   { 
   // To avoid the circular dependency between G4Track, G4Step
   // and G4StepPoint, G4Step has to manage the update actions.
   fpTrack->SetPosition(fpPostStepPoint->GetPosition());
   fpTrack->SetGlobalTime(fpPostStepPoint->GetGlobalTime());
   fpTrack->SetLocalTime(fpPostStepPoint->GetLocalTime());
   fpTrack->SetProperTime(fpPostStepPoint->GetProperTime());
   fpTrack->SetMomentumDirection(fpPostStepPoint->GetMomentumDirection());
   fpTrack->SetKineticEnergy(fpPostStepPoint->GetKineticEnergy());
   fpTrack->SetPolarization(fpPostStepPoint->GetPolarization());
   fpTrack->SetStepLength(fStepLength);
   // NextTouchable is updated 
   // (G4Track::Touchable points touchable of Pre-StepPoint)
   fpTrack->SetNextTouchable(fpPostStepPoint->GetTouchable());
   fpTrack->SetWeight(fpPostStepPoint->GetWeight());
   }

#endif










