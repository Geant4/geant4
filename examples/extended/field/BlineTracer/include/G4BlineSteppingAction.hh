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
// $Id: G4BlineSteppingAction.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 
// --------------------------------------------------------------------
//
// G4BlineSteppingAction
//
// Class description:
//
// Defines the stepping action used when tracing magnetic field lines.

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlineSteppingAction_h
#define G4BlineSteppingAction_h 1

#include "G4Types.hh"
#include "G4UserSteppingAction.hh"

class G4BlineTracer;

class G4BlineSteppingAction : public G4UserSteppingAction
{
  public:  // with description

    G4BlineSteppingAction(G4BlineTracer* aBlineTool);
    virtual ~G4BlineSteppingAction();
    virtual void UserSteppingAction(const G4Step*);
    
  private:

    G4BlineTracer* fBlineTool;
};

#endif
