// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DalitzDecayChannel.hh,v 1.1 1999-01-07 16:10:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
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
  public:
    //Constructors 
      G4DalitzDecayChannel(const G4String& theParentName,
			   G4double        theBR,
			   const G4String& theLeptonName,
			   const G4String& theAntiLeptonName);
    //  Destructor
      virtual ~G4DalitzDecayChannel();

  public:
     virtual G4DecayProducts *DecayIt(G4double);     

  private:
     enum{idGamma=0, idLepton=1, idAntiLepton=2}; 

};  


#endif
