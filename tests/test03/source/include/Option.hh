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
// $Id$
//
/// \file Options.hh
/// \brief Definition of the Options class

#ifndef Options_h
#define Options_h 1

#include "G4AnyType.hh"

template<class T>
class Option<T>
{
  public:
    template<class T>
    Option<t> (const G4String& name,
               const G4String& shortName,
               const G4String& description,
               T defaultValue)
     : fName(name),
       fShortName(shortName),
       fDescription(description),
       fValue(defaultValue) {}          
    
    // methods
    G4String GetName() const { return fName; }
    G4String GetShortName() const { return fShortName; }
    G4String GetDescription() const { return fDescription; }

    template<class T>
    T GetValue() const { return fValue; }

  private:
    template<class T>
    Option() {}

    G4String fName;
    G4String fShortName;
    G4String fDescription;
    T        fValue;
};

#endif

    
