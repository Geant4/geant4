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

#ifndef G4ATTVALUE_HH
#define G4ATTVALUE_HH

#include "globals.hh"

// Class Description:
//
// @class G4AttValue
//
// @brief This class represents a HepRep-style Attribute Value.
// G4AttValues can be attached to a Trajectory, Trajectory Point or Sensitive
// Detector Hit.  These attributes are then made available to the end user
// in interactive visualization systems (such as WIRED).
// The G4AttValue is further defined in a G4AttDef object.
// The association between the G4AttValue and the G4AttDef object is
// made through the data member "name".
// For details, see the HepRep home page at http://heprep.freehep.org
//  
// @author M.Frailis 
// @author R.Giannitrapani
// @author J.Perl
// Class Description - End:


  class G4AttValue {
    
  public:
    G4AttValue(const G4String& name,
               const G4String& value,
               const G4String& showLabel): 
      m_name(name),m_value(value),
      m_showLabel(showLabel){};
    G4AttValue()= default;
    
    const G4String& GetName()const{return m_name;};
    const G4String& GetValue()const{return m_value;};
    const G4String& GetShowLabel()const{return m_showLabel;};

    void SetName(const G4String& name){m_name = name;};
    void SetValue(const G4String& val){m_value = val;};
    void SetShowLabel(const G4String& lab){m_showLabel = lab;};

  private:
    /// The name of the attribute 
    G4String m_name;
    /// The value of the attribute
    G4String m_value;
    /// The bitmap for the label display
    G4String m_showLabel;
  };

#endif //G4ATTVALUE_H
