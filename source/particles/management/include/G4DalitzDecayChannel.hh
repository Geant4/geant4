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
// $Id: G4DalitzDecayChannel.hh,v 1.4 2001-07-11 10:01:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4DalitzDecayChannel_h
#define G4DalitzDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4DalitzDecayChannel :public G4VDecayChannel
{
 // Class Description
 //   This class describes kinematics in Dalitz decay
 //       parent -> lepton + anti_lepton   
 //

 public: // With Description
    //Constructors 
      G4DalitzDecayChannel(const G4String& theParentName,
			   G4double        theBR,
			   const G4String& theLeptonName,
			   const G4String& theAntiLeptonName);
    //  Destructor
      virtual ~G4DalitzDecayChannel();

  public: // With Description
     virtual G4DecayProducts *DecayIt(G4double);     

  private:
     enum{idGamma=0, idLepton=1, idAntiLepton=2}; 

};  


#endif
