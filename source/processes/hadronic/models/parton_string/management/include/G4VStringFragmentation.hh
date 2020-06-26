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
#ifndef G4VStringFragmentation_h
#define G4VStringFragmentation_h 1

#include "G4HadronicInteraction.hh"
#include "G4ExcitedStringVector.hh"

class G4KineticTrackVector;

class G4VStringFragmentation : public G4HadronicInteraction
{
  public:
    G4VStringFragmentation(const G4String& name = "StringFragmentation");
    ~G4VStringFragmentation() override;

    G4VStringFragmentation(const G4VStringFragmentation &right) = delete;
    const G4VStringFragmentation & operator=
    (const G4VStringFragmentation &right) = delete;
    G4bool operator==(const G4VStringFragmentation &right) const = delete;
    G4bool operator!=(const G4VStringFragmentation &right) const = delete;

    virtual G4KineticTrackVector * FragmentStrings(const G4ExcitedStringVector * theStrings)=0;

};

#endif

