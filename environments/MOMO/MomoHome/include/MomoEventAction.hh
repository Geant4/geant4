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
// Taken from
// $Id: MomoEventAction.hh,v 1.1 2003/11/04 09:02:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
 
#ifndef MomoEventAction_h
#define MomoEventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;


class MomoEventAction : public G4UserEventAction
{
  public:
    MomoEventAction();
   ~MomoEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
};


#endif

    
