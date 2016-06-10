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
// $Id: G4AttDefT.hh 78955 2014-02-05 09:45:46Z gcosmo $
//
// Templated G4AttDef. Generates type key for given template parameter.
//
// Jane Tinslay, September 2006
//
#ifndef G4ATTDEFT_HH
#define G4ATTDEFT_HH

#include "globals.hh"
#include "G4AttDef.hh"
#include "G4TypeKeyT.hh"

template <typename T>  
class G4AttDefT : public G4AttDef {

public:
  
  typedef T Type;

  // Constructor
  G4AttDefT(const G4String& name,
            const G4String& desc,
            const G4String& category,
            const G4String& extra="")    
    :G4AttDef(name, desc, category, extra, G4TypeKeyT<T>())
  {}
  
  // Destructor
  virtual ~G4AttDefT() {};
  
private:

};

#endif
