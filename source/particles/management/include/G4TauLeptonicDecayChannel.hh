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
// $Id: G4TauLeptonicDecayChannel.hh,v 1.1 2002-03-08 08:47:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4TauLeptonicDecayChannel_h
#define G4TauLeptonicDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4TauLeptonicDecayChannel :public G4VDecayChannel
{
  // Class Decription
  //  This class describes tau leptonic decay kinemtics.
  //  This version assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
  //              without tau polarization 

  public:  // With Description
    //Constructors 
      G4TauLeptonicDecayChannel(const G4String& theParentName,
				G4double        theBR,
				const G4String& theLeptonName);
    //  Destructor
      virtual ~G4TauLeptonicDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);   
  
  private:
     static G4double   spectrum(G4double momentum,
				G4double energy,
				G4double mtau,
				G4double ml);
};  


#endif

