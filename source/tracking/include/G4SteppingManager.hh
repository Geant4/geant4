// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingManager.hh,v 1.15 2001-02-08 07:48:39 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//---------------------------------------------------------------
//
// G4SteppingManager.hh
//
// class description:
//  This is the class which plays an essential role in tracking 
//  the particle. It takes cares all message passing between
//  objects in the different categories (for example, 
//  geometry(including transportation), interactions in 
//  matter, etc). It's public method 'stepping' steers to step 
//  the particle.
//  Geant4 kernel use only
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige


class G4SteppingManager;

#ifndef G4SteppingManager_h
#define G4SteppingManager_h 1

#include "G4ios.hh"                 // Include from 'system'
#include "g4std/iomanip"                  // Include from 'system'
#include "g4std/vector"                  // Include from 'system'
#include "globals.hh"                 // Include from 'global'
#include "Randomize.hh"               // Include from 'global'
#include "G4Navigator.hh"             // Include from 'geometry'
#include "G4LogicalVolume.hh"         // Include from 'geometry'
#include "G4VPhysicalVolume.hh"       // Include from 'geometry'
//#include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'
class G4VSensitiveDetector;
#include "G4ProcessManager.hh"        // Include from 'piim'
//#include "G4DynamicParticleVector.hh" // Include from 'piim'
//#include "G4Hit.hh"                   // Include from 'Hit/dig'

#include "G4Track.hh"                 // Include from 'tracking'
#include "G4TrackVector.hh"           // Include from 'tracking'
#include "G4TrackStatus.hh"           // Include from 'tracking'
#include "G4StepStatus.hh"            // Include from 'tracking'
#include "G4UserSteppingAction.hh"    // Include from 'tracking'
#include "G4Step.hh"                  // Include from 'tracking'
#include "G4StepPoint.hh"             // Include from 'tracking'
#include "G4VSteppingVerbose.hh"       // Include from 'tracking'
#include "G4VTouchable.hh"            // Include from 'geometry'
#include "G4TouchableHistory.hh"      // Include from 'geometry'

//  must be changed in elegant way
// RogueWave Tools.h++
   typedef G4std::vector<G4int> 
             G4SelectedAtRestDoItVector;
   typedef G4std::vector<G4int> 
             G4SelectedAlongStepDoItVector;
   typedef G4std::vector<G4int>
             G4SelectedPostStepDoItVector;


///////////////////////
class G4SteppingManager 
///////////////////////
{
//--------
public: //without description
//--------

// Constructor/Destructor

   G4SteppingManager();
      // SteppingManger should be dynamically persistent, therefore 
      // you need to invoke new() when you call this constructor.
      // "Secodary track vector" will be dynamically created by this 
      // cosntructor. "G4UserSteppingAction" will be also constructed 
      // in this constructor, and "this" pointer will be passed to 
      // "G4UserSteppingAction". 

   ~G4SteppingManager();

// Get/Set functions

   G4TrackVector* GetSecondary() const;
   void SetUserAction(G4UserSteppingAction* apAction);
   G4Track* GetTrack() const;
   void SetVerboseLevel(G4int vLevel);
   void SetVerbose(G4VSteppingVerbose*);
   G4Step* GetStep() const;


// Other member functions

   G4StepStatus Stepping();
      // Steers to move the give particle from the TrackingManger 
      // by one Step.

  void SetInitialStep(G4Track* valueTrack);
     // Sets up initial track information (enegry, position, etc) to 
     // the PreStepPoint of the G4Step. This funciton has to be called 
     // just once before the stepping loop in the "TrackingManager".

  void GetProcessNumber();

// Get methods
   G4double GetPhysicalStep();
   G4double GetGeometricalStep();
   G4double GetCorrectedStep();
   G4bool GetPreStepPointIsGeom();
   G4bool GetFirstStep();
   G4StepStatus GetfStepStatus();
   G4double GetTempInitVelocity();
   G4double GetTempVelocity();
   G4double GetMass();
   G4double GetsumEnergyChange();
   G4VParticleChange* GetfParticleChange();
   G4Track* GetfTrack();
   G4TrackVector* GetfSecondary();
   G4Step* GetfStep();
   G4StepPoint* GetfPreStepPoint();
   G4StepPoint* GetfPostStepPoint();
   G4VPhysicalVolume* GetfCurrentVolume();
   G4VSensitiveDetector* GetfSensitive();
   G4VProcess* GetfCurrentProcess();
   G4ProcessVector* GetfAtRestDoItVector();
   G4ProcessVector* GetfAlongStepDoItVector();
   G4ProcessVector* GetfPostStepDoItVector();
   G4ProcessVector* GetfAlongStepGetPhysIntVector();
   G4ProcessVector* GetfPostStepGetPhysIntVector();
   G4ProcessVector* GetfAtRestGetPhysIntVector();
   G4double GetcurrentMinimumStep();
   G4double GetnumberOfInteractionLengthLeft();
   size_t GetfAtRestDoItProcTriggered();
   size_t GetfAlongStepDoItProcTriggered();
   size_t GetfPostStepDoItProcTriggered();
   G4int GetfN2ndariesAtRestDoIt();
   G4int GetfN2ndariesAlongStepDoIt();
   G4int GetfN2ndariesPostStepDoIt();
   G4Navigator* GetfNavigator();
   G4int GetverboseLevel();
   size_t GetMAXofAtRestLoops();
   size_t GetMAXofAlongStepLoops();
   size_t GetMAXofPostStepLoops();
   G4SelectedAtRestDoItVector* GetfSelectedAtRestDoItVector();
   G4SelectedAlongStepDoItVector* GetfSelectedAlongStepDoItVector();
   G4SelectedPostStepDoItVector* GetfSelectedPostStepDoItVector();
   G4double   GetfPreviousStepSize();
   const G4VTouchable* GetfTouchable1();
   const G4VTouchable* GetfTouchable2();
   G4bool GetfIsTouchable1Free();
   G4bool GetfIsTouchable2Free();
   G4SteppingControl GetStepControlFlag();
   G4Navigator GetNavigator();
   G4UserSteppingAction* GetUserAction();
   G4double GetphysIntLength();
   G4ForceCondition GetfCondition();
   G4GPILSelection  GetfGPILSelection();
  //
   G4bool KillVerbose;
//---------   
   private:
//---------   

// Member functions

   void DefinePhysicalStepLength();
      // Calculate corresponding physical length from the mean free path 
      // left for each discrete phyiscs process. The minimum allowable
      // Step for each continious process will be also calculated.
   void InvokeAtRestDoItProcs();
   void InvokeAlongStepDoItProcs();
   void InvokePostStepDoItProcs();
   void SetNavigator(G4Navigator* value);
   G4double CalculateSafety();
      // Return the estimated safety value at the PostStepPoint
   const G4VTouchable* GetFreeTouchable();
      // Get Touchable which is free, i.e. not assigined to Track/StepPoint
      // If no free touchable is availabe, the NULL will be returned
      // Once you get a Touchable, it will be set to NotFree. 
   void SetAnotherTouchableFree(const G4VTouchable* value);
      // Set the partner of the given Touchable to be free. For example,
      // the argument has fTouchable1, then fTouchable2 will be set to 
      // free. The state of fTouchable1 is intact.


// Member data
   static const size_t SizeOfSelectedDoItVector=100;

   G4UserSteppingAction* fUserSteppingAction;

   G4VSteppingVerbose* fVerbose;

   G4double PhysicalStep;
   G4double GeometricalStep;
   G4double CorrectedStep;
   G4bool PreStepPointIsGeom;
   G4bool FirstStep;
   G4StepStatus fStepStatus;

   G4double TempInitVelocity;
   G4double TempVelocity;
   G4double Mass;

   G4double sumEnergyChange;

   G4VParticleChange* fParticleChange;
   G4Track* fTrack;
   G4TrackVector* fSecondary;
   G4Step* fStep;
   G4StepPoint* fPreStepPoint;
   G4StepPoint* fPostStepPoint;

   G4VPhysicalVolume* fCurrentVolume;
   G4VSensitiveDetector* fSensitive;
   G4VProcess* fCurrentProcess;
      // The pointer to the process of which DoIt or
      // GetPhysicalInteractionLength has been just executed.


   G4ProcessVector* fAtRestDoItVector;
   G4ProcessVector* fAlongStepDoItVector;
   G4ProcessVector* fPostStepDoItVector;

   G4ProcessVector* fAtRestGetPhysIntVector;
   G4ProcessVector* fAlongStepGetPhysIntVector;
   G4ProcessVector* fPostStepGetPhysIntVector;

   size_t MAXofAtRestLoops;
   size_t MAXofAlongStepLoops;
   size_t MAXofPostStepLoops;

   G4double currentMinimumStep;
   G4double numberOfInteractionLengthLeft;

   size_t fAtRestDoItProcTriggered;
   size_t fAlongStepDoItProcTriggered;
   size_t fPostStepDoItProcTriggered;

   G4int fN2ndariesAtRestDoIt;
   G4int fN2ndariesAlongStepDoIt;
   G4int fN2ndariesPostStepDoIt;
      // These are the numbers of secondaries generated by the process
      // just executed.

   G4Navigator *fNavigator;

   G4int verboseLevel;

   G4SelectedAtRestDoItVector* fSelectedAtRestDoItVector;
   G4SelectedAlongStepDoItVector* fSelectedAlongStepDoItVector;
   G4SelectedPostStepDoItVector* fSelectedPostStepDoItVector;

   G4double   fPreviousStepSize;

   const G4VTouchable* fTouchable1;
   const G4VTouchable* fTouchable2;
   G4bool fIsTouchable1Free;
   G4bool fIsTouchable2Free;

   G4SteppingControl StepControlFlag;

   G4double proposedSafety;
      // This keeps the minimum safety value proposed by AlongStepGPILs.
   G4ThreeVector endpointSafOrigin;
   G4double endpointSafety;
      // To get the true safety value at the PostStepPoint, you have
      // to subtract the distance to 'endpointSafOrigin' from this value.
   G4double physIntLength;
   G4ForceCondition fCondition;
   G4GPILSelection  fGPILSelection;
      // Above three variables are for the method 
      // DefinePhysicalStepLength(). To pass these information to
      // the method Verbose, they are kept at here. Need a more 
      // elegant mechanism.

};


//*******************************************************************
//
//  Inline function 
//
//*******************************************************************

  inline G4double G4SteppingManager::GetPhysicalStep(){
    return PhysicalStep;
  }

  inline  G4double G4SteppingManager::GetGeometricalStep(){
    return GeometricalStep;
  }

  inline  G4double G4SteppingManager::GetCorrectedStep(){
    return CorrectedStep;
  }
   inline G4bool G4SteppingManager::GetPreStepPointIsGeom(){
     return PreStepPointIsGeom;
   }

   inline G4bool G4SteppingManager::GetFirstStep(){
     return FirstStep;
   }

   inline G4StepStatus G4SteppingManager::GetfStepStatus(){
     return fStepStatus;
   }

   inline G4double G4SteppingManager::GetTempInitVelocity(){
     return TempInitVelocity;
   }
   inline G4double G4SteppingManager::GetTempVelocity(){
     return TempVelocity;
   }
   inline G4double G4SteppingManager::GetMass(){
     return Mass;
   }

   inline G4double G4SteppingManager::GetsumEnergyChange(){
     return sumEnergyChange;
   }

   inline G4VParticleChange* G4SteppingManager::GetfParticleChange(){
     return fParticleChange;
   }

   inline G4Track* G4SteppingManager::GetfTrack(){
     return fTrack;
   }

   inline G4TrackVector* G4SteppingManager::GetfSecondary(){
     return fSecondary;
   }
   inline G4Step* G4SteppingManager::GetfStep(){
     return fStep;
   }
   inline G4StepPoint* G4SteppingManager::GetfPreStepPoint(){
     return fPreStepPoint;
   }
   inline G4StepPoint* G4SteppingManager::GetfPostStepPoint(){
     return fPostStepPoint;
   }

   inline G4VPhysicalVolume* G4SteppingManager::GetfCurrentVolume(){
     return fCurrentVolume;
   }
   inline G4VSensitiveDetector* G4SteppingManager::GetfSensitive(){
     return fSensitive;
   }
   inline G4VProcess* G4SteppingManager::GetfCurrentProcess(){
     return fCurrentProcess;
   }

   inline G4ProcessVector* G4SteppingManager::GetfAtRestDoItVector(){
     return fAtRestDoItVector;
   }
   inline G4ProcessVector* G4SteppingManager::GetfAlongStepDoItVector(){
     return fAlongStepDoItVector;
   }
   inline G4ProcessVector* G4SteppingManager::GetfPostStepDoItVector(){
     return fPostStepDoItVector;
   }

   inline G4ProcessVector* G4SteppingManager::GetfAtRestGetPhysIntVector(){
     return fAtRestGetPhysIntVector;
   }

   inline G4ProcessVector* G4SteppingManager::GetfAlongStepGetPhysIntVector(){
     return fAlongStepGetPhysIntVector;
   }

   inline G4ProcessVector* G4SteppingManager::GetfPostStepGetPhysIntVector(){
     return fPostStepGetPhysIntVector;
   }

   inline size_t G4SteppingManager::GetMAXofAtRestLoops(){
     return MAXofAtRestLoops;
   }
   inline size_t G4SteppingManager::GetMAXofAlongStepLoops(){
     return MAXofAlongStepLoops;
   }
   inline size_t G4SteppingManager::GetMAXofPostStepLoops(){
     return MAXofPostStepLoops;
   }

   inline G4double G4SteppingManager::GetcurrentMinimumStep(){
     return currentMinimumStep;
   }
   inline G4double G4SteppingManager::GetnumberOfInteractionLengthLeft(){
     return numberOfInteractionLengthLeft;
   }

   inline size_t G4SteppingManager::GetfAtRestDoItProcTriggered(){
     return fAtRestDoItProcTriggered;
   }
   inline size_t G4SteppingManager::GetfAlongStepDoItProcTriggered(){
     return fAtRestDoItProcTriggered;
   }
   inline size_t G4SteppingManager::GetfPostStepDoItProcTriggered(){
     return fPostStepDoItProcTriggered;
   }
   inline G4int G4SteppingManager::GetfN2ndariesAtRestDoIt(){
     return fN2ndariesAtRestDoIt;
   }
   inline G4int G4SteppingManager::GetfN2ndariesAlongStepDoIt(){
     return fN2ndariesAlongStepDoIt;
   }
   inline G4int G4SteppingManager::GetfN2ndariesPostStepDoIt(){
     return fN2ndariesPostStepDoIt;
   }

   inline G4Navigator* G4SteppingManager::GetfNavigator(){
     return fNavigator;
   }
   inline G4int G4SteppingManager::GetverboseLevel(){
     return verboseLevel;
   }

   inline G4SelectedAtRestDoItVector* G4SteppingManager::GetfSelectedAtRestDoItVector(){
     return fSelectedAtRestDoItVector;
   }

   inline G4SelectedAlongStepDoItVector* G4SteppingManager::GetfSelectedAlongStepDoItVector(){
     return fSelectedAlongStepDoItVector;
   }

   inline G4SelectedPostStepDoItVector* G4SteppingManager::GetfSelectedPostStepDoItVector(){
     return fSelectedPostStepDoItVector;
   }

   inline G4double   G4SteppingManager::GetfPreviousStepSize(){
     return fPreviousStepSize;
   }

   inline const G4VTouchable* G4SteppingManager::GetfTouchable1(){
     return fTouchable1;
   }
   inline const G4VTouchable* G4SteppingManager::GetfTouchable2(){
     return fTouchable2;
   }
   inline G4bool G4SteppingManager::GetfIsTouchable1Free(){
     return fIsTouchable1Free;
   }
   inline G4bool G4SteppingManager::GetfIsTouchable2Free(){
     return fIsTouchable2Free;
   }

   inline G4SteppingControl G4SteppingManager::GetStepControlFlag(){
     return StepControlFlag;
   }
   inline G4double G4SteppingManager::GetphysIntLength(){
     return physIntLength;
   }
   inline G4ForceCondition G4SteppingManager::GetfCondition(){
     return fCondition;
   }
   inline G4GPILSelection  G4SteppingManager::GetfGPILSelection(){
     return fGPILSelection;
   }


  inline G4TrackVector* G4SteppingManager::GetSecondary() const {
    return fSecondary; 
  }

  inline void G4SteppingManager::SetNavigator(G4Navigator* value){
    fNavigator = value; 
  }

  inline void G4SteppingManager::SetUserAction(G4UserSteppingAction* apAction){
    if (apAction != NULL) fUserSteppingAction = apAction;
  }
  inline G4UserSteppingAction* G4SteppingManager::GetUserAction(){
    return fUserSteppingAction;
  }

  inline G4Track* G4SteppingManager::GetTrack() const {
    return fTrack; 
  }

  inline void G4SteppingManager::SetVerboseLevel(G4int vLevel){
    verboseLevel = vLevel; 
  }

  inline void G4SteppingManager::SetVerbose(G4VSteppingVerbose* yourVerbose){
     fVerbose = yourVerbose;
  }

  inline G4Step* G4SteppingManager::GetStep() const {
    return fStep;
  }

  inline G4double G4SteppingManager::CalculateSafety(){
    return G4std::max( endpointSafety -
		(endpointSafOrigin - fPostStepPoint->GetPosition()).mag(),
	        0.);
  }

  inline const G4VTouchable* G4SteppingManager::GetFreeTouchable(){
    if(fIsTouchable1Free){
      fIsTouchable1Free = false;
      return fTouchable1;
    } 
    else if(fIsTouchable2Free){
      fIsTouchable2Free = false;
      return fTouchable2;
    }
    else {
      return NULL;
    }
  }

  inline void G4SteppingManager::SetAnotherTouchableFree(const G4VTouchable* value){
    if(value == fTouchable1){
      fIsTouchable2Free = true;
    }
    else if(value == fTouchable2){
      fIsTouchable1Free = true;
    }
  }

//#include "G4SteppingManager.icc"

#endif





