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
// $Id: G4NeutronBetaDecayChannel.hh,v 1.1 2001-09-18 02:15:10 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      18 Sep. 2001 H.Kurashige
// ------------------------------------------------------------
#ifndef G4NeutronBetaDecayChannel_h
#define G4NeutronBetaDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4NeutronBetaDecayChannel :public G4VDecayChannel
{
  // Class Decription
  //  This class describes free neutron beta decay  kinemtics.
  //  This version neglects neutron/electron polarization  
  //  without Coulomb effect

  public:  // With Description
    //Constructors 
      G4NeutronBetaDecayChannel(const G4String& theParentName,
				G4double        theBR);
    //  Destructor
      virtual ~G4NeutronBetaDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);     
  
  protected: 
  // e-neutrino angular correlation parameter 
     const G4double aENuCorr;
};  


#endif
