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
// $Id: MyRunAction.hh,v 1.4 2001-07-11 10:10:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyRunAction_h
#define MyRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"

#ifdef G4VIS_USE_OPACS
#include <OHistogram.h>
#endif

class G4Run;

class MyRunAction : public G4UserRunAction
{
  public:
    MyRunAction();
    ~MyRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);
#ifdef G4VIS_USE_OPACS
    static OHistogram get_1d(){return myHisto1D; }
    static OHistogram Get2d(){return myHisto2D; }
#endif

  private:
    G4Timer* timer;
#ifdef G4VIS_USE_OPACS
    static OHistogram myHisto1D;
    static OHistogram myHisto2D;
#endif
};

#endif

