// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserSteppingAction.hh,v 1.1 1999-01-07 16:14:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4UserSteppingAction.hh
//
//  Description:
//    This class represents actions taken place by the user at each
//    end of stepping. 
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

class G4UserSteppingAction;

#ifndef G4UserSteppingAction_h
#define G4UserSteppingAction_h 1

#include "globals.hh"                  // Include from 'global'
#include "G4Track.hh"                  // Include from 'tracking'
class G4SteppingManager;               // Forward declaration

///////////////////////////
class G4UserSteppingAction 
///////////////////////////
{

//--------
   public:
//--------

// Constructor and destructors
   G4UserSteppingAction(){}
   virtual ~G4UserSteppingAction(){}

// Member functions
   void SetSteppingManagerPointer(G4SteppingManager* pValue);
   virtual void UserSteppingAction(){}

//----------- 
   protected:
//----------- 

// Member functions
   inline const G4SteppingManager* GetSteppingManager() {
      return fpSteppingManager;
   }
   inline G4SteppingManager* GetOmnipotentSteppingManager() {
      return fpSteppingManager;
   }

//---------
   private:
//---------

// Member data
   G4SteppingManager* fpSteppingManager;

};

#endif


