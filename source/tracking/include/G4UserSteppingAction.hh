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
// $Id: G4UserSteppingAction.hh,v 1.7 2001-07-11 10:08:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4UserSteppingAction.hh
//
// class description:
//    This class represents actions taken place by the user at each
//    end of stepping. 
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4UserSteppingAction_h
#define G4UserSteppingAction_h 1

class G4Step;
class G4SteppingManager;               // Forward declaration

///////////////////////////
class G4UserSteppingAction 
///////////////////////////
{

//--------
public: // with description
//--------

// Constructor and destructors
   G4UserSteppingAction(){;}
   virtual ~G4UserSteppingAction(){;}

// Member functions
   void SetSteppingManagerPointer(G4SteppingManager* pValue);
   virtual void UserSteppingAction(const G4Step* aStep){;}

//-----------
   protected:
//-----------

// Member data
   G4SteppingManager* fpSteppingManager;

};

#endif


