#ifndef G4VGFlashSensitiveDetector_h
#define G4VGFlashSensitiveDetector_h 1

#include "G4Step.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "GFlashEnergySpot.hh"
#include "G4GFlashSpot.hh"

// class description:
//
//  This is the abstract base class of the sensitive detector for use with GFlash. The user's
// sensitive detector which generates hits must be derived from this
// class, and G4VSensitiveDetector.

class G4VGFlashSensitiveDetector 
{

  public: // with description
      G4VGFlashSensitiveDetector() {}
      G4VGFlashSensitiveDetector(const G4VGFlashSensitiveDetector &right) {}
      // Constructors. The user's concrete class must use one of these constructors
      // by the constructor initializer of the derived class. The name of
      // the sensitive detector must be the same as for the corresponding GG4VSensitiveDetector.

  public:
      virtual ~G4VGFlashSensitiveDetector() {}

      G4int operator==(const G4VGFlashSensitiveDetector &right) const {return this == &right;}
      G4int operator!=(const G4VGFlashSensitiveDetector &right) const {return this != &right;};

  protected: // with description
      virtual G4bool ProcessHits(G4GFlashSpot*aSpot,G4TouchableHistory*ROhist) = 0;
      //  The user MUST implement this method for generating hit(s) from the GFlashSpots 
      //  Be aware that this method is a protected method and it will be invoked 
      // by Hit() method of the Base class once the Readout geometry that may be associated to the
      // corresponding G4VSensitiveDetector was taken into account.


  public: // with description
      inline G4bool Hit(G4GFlashSpot * aSpot)
      {
        G4bool result = true; 
	  	  G4VSensitiveDetector * This = dynamic_cast<G4VSensitiveDetector *>(this);
	  		if(!This)
	  		{
	    		G4Exception("Sensitive detector using GFlash needs to also inherit from G4VSensitiveDetector");
	  		}
        if(This->isActive())
        { 
 	        G4VReadOutGeometry * ROgeometry = 0;
          G4TouchableHistory* ROhis = 0;

	  			if(This) ROgeometry = This->GetROgeometry();
	  			if(ROgeometry)
	  			{
	    			// fake pre-step point for touchable from read-out geometry.
	    			G4Step fakeStep;
	    			G4StepPoint * tmpPoint = fakeStep.GetPreStepPoint();
	    			tmpPoint->SetTouchableHandle(aSpot->GetTouchableHandle());
	    			tmpPoint->SetPosition(aSpot->GetPosition());
	    			tmpPoint->SetMomentumDirection(aSpot->GetOriginatorTrack()->GetPrimaryTrack()->GetMomentumDirection());
  	    		result = ROgeometry->CheckROVolume(&fakeStep, ROhis); 
	  			}
	  			if(result) result = ProcessHits(aSpot, ROhis); 
				}
				else 
				{
	  			result = false;
				}
				return result;
      }
      //  This is the public method invoked by GFlashHitMaker for generating
      // hit(s). The actual user's implementation for generating hit(s) must be
      // implemented in GenerateHits() virtual protected method. 
};

#endif

