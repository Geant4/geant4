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
// Author: Ivana Hrivnacova, 11/09/2018  (ivana@ipno.in2p3.fr)

#ifndef G4VScoreNtupleWriter_h
#define G4VScoreNtupleWriter_h 1

#include "G4Threading.hh"
#include "globals.hh"

class G4HCofThisEvent;

// class description:
//
// This class implements the interface for storing hits collections of G4THitsMap<G4double>
// type vith Geant4 analysis tools.

class G4VScoreNtupleWriter
{
  public: 
    virtual ~G4VScoreNtupleWriter();

    // static methods
    static G4VScoreNtupleWriter* Instance();

    // methods
    virtual G4bool Book(G4HCofThisEvent* hce) = 0;
    virtual void   OpenFile() = 0;
    virtual void   Fill(G4HCofThisEvent* hce, G4int eventNumber) = 0;
    virtual void   Write() = 0;

  protected:
    G4VScoreNtupleWriter();
    virtual G4VScoreNtupleWriter* CreateInstance() const = 0;

    // static data members
    static G4VScoreNtupleWriter* fgMasterInstance;    
    static G4ThreadLocal G4VScoreNtupleWriter* fgInstance;    
};

#endif
