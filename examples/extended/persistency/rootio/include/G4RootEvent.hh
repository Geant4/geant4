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
// File: G4RootEvent.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef ROOT_FADS_EVENT_HH
#define ROOT_FADS_EVENT_HH 1

class G4Event;
class GenEvent;
class G4MCTEvent;
class G4Pevent;

// Class inherited:
#include "TObject.h"

// Class Description:
//   FADS persistent event object.
// 
//   This class has pointers to HepMCGenEvent, MCTruth and G4Event.

class G4RootEvent
 : public TObject
{
    public: // With description
      G4RootEvent() {};
      // Constructor

      virtual ~G4RootEvent();
      // Destructor

    public: // With description
      int GetEventID() { return m_id; };
      // returns the event ID.

      G4Pevent* MakeTransientObject();
      // Make transient G4Pevent from G4RootEvent

    private:
      int m_id;

    ClassDef (G4RootEvent,1)

}; // End of class G4RootEvent

#endif

