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
// $Id: RE01RunAction.hh,v 1.1 2004/11/26 07:37:40 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 

#ifndef RE01RunAction_h
#define RE01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;

class RE01RunAction : public G4UserRunAction
{
  public:
    RE01RunAction();
   ~RE01RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
};


#endif





