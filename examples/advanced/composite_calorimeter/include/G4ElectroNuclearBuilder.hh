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
#ifndef G4ElectroNuclearBuilder_h
#define G4ElectroNuclearBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4GammaNuclearReaction.hh"
#include "G4ElectroNuclearReaction.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

class G4ElectroNuclearBuilder 
{
  public: 
    G4ElectroNuclearBuilder();
    virtual ~G4ElectroNuclearBuilder();

  public: 
    virtual void Build();

  protected:
    G4PhotoNuclearProcess thePhotoNuclearProcess;
    G4ElectronNuclearProcess theElectronNuclearProcess;
    G4PositronNuclearProcess thePositronNuclearProcess;
    G4ElectroNuclearReaction * theElectroReaction;
    G4GammaNuclearReaction * theGammaReaction;  
    
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4GammaParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
};

// 2002 by J.P. Wellisch

#endif





