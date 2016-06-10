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
// $Id: G4VAttValueFilter.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Abstract base class for G4AttValue filters
//
// Jane Tinslay, September 2006
//
#ifndef G4VATTVALUEFILTER_HH
#define G4VATTVALUEFILTER_HH

#include "globals.hh"
#include "G4String.hh"
#include "G4VFilter.hh"

class G4AttValue;

class G4VAttValueFilter : public G4VFilter<G4AttValue> {

public:

  // Constructor
  G4VAttValueFilter(const G4String& name = "G4AttValueFilter")
    :G4VFilter<G4AttValue>(name){}

  // Destructor
  virtual ~G4VAttValueFilter() {}
  
  // Filter methods
  virtual G4bool Accept(const G4AttValue&) const = 0;
  virtual G4bool GetValidElement(const G4AttValue&, G4String&) const = 0;

  // Print configuration
  virtual void PrintAll(std::ostream& ostr) const = 0;
  
  // Reset 
  virtual void Reset() = 0;

  // Load filter data
  virtual void LoadIntervalElement(const G4String&) = 0;
  virtual void LoadSingleValueElement(const G4String&) = 0;

};

#endif
