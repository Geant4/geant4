//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
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

  protected:
    // Copy constructor and assignment operator
      G4TauLeptonicDecayChannel(const G4TauLeptonicDecayChannel &);
      G4TauLeptonicDecayChannel & operator=(const G4TauLeptonicDecayChannel &);

  protected:
      G4TauLeptonicDecayChannel();


  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);   
  
  private:
     static G4double   spectrum(G4double momentum,
				G4double energy,
				G4double mtau,
				G4double ml);
};  


#endif

