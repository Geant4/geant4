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
// $Id: G4MuonDecayChannel.hh,v 1.5 2001-07-11 10:01:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4MuonDecayChannel_h
#define G4MuonDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4MuonDecayChannel :public G4VDecayChannel
{
  // Class Decription
  //  This class describes muon decay kinemtics.
  //  This version neglects muon polarization  
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
  //

  public:  // With Description
    //Constructors 
      G4MuonDecayChannel(const G4String& theParentName,
			 G4double        theBR);
    //  Destructor
      virtual ~G4MuonDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);     
  
};  


#endif
