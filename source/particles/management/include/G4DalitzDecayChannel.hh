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
// G4DalitzDecayChannel
//
// Class description:
//
// Class describing the kinematics in Dalitz decay:
// parent -> lepton + anti_lepton

// Author: H.Kurashige, 30 May 1997 
// --------------------------------------------------------------------
#ifndef G4DalitzDecayChannel_hh
#define G4DalitzDecayChannel_hh 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4DalitzDecayChannel : public G4VDecayChannel
{

  public:

    G4DalitzDecayChannel(const G4String& theParentName,
                               G4double  theBR,
                         const G4String& theLeptonName,
                         const G4String& theAntiLeptonName);

    virtual ~G4DalitzDecayChannel();

    virtual G4DecayProducts* DecayIt(G4double);     

  protected:

    G4DalitzDecayChannel(const G4DalitzDecayChannel&);
    G4DalitzDecayChannel& operator=(const G4DalitzDecayChannel&);
      // Copy constructor and assignment operator

  private:

    G4DalitzDecayChannel();

    enum{idGamma=0, idLepton=1, idAntiLepton=2};
};

#endif
