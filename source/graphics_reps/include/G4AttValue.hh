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
// $Id: G4AttValue.hh,v 1.3 2002-10-24 14:35:20 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
    G4AttValue(){};
    
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
