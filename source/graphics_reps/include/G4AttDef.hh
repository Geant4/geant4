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
// $Id: G4AttDef.hh,v 1.3 2002-10-24 14:35:14 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4ATTDEF_HH
#define G4ATTDEF_HH

#include "globals.hh"

// Class Description:
//
// @class G4AttDef
//
// @brief This class represents a HepRep-style Attribute Definition.
// The G4AttDef is used to define new kinds of attributes that can
// then have values set for a Trajectory, Trajectory Point or Sensitive
// Detector Hit.  These attributes are then made available to the end user
// in an interactive visualization system (such as WIRED).
// Values are set by creating G4AttValue objects and attaching them to the
// relevant Trajectory, Trajectory Point or Sensitive Detector Hit.
// The association between the G4AttValue and the G4AttDef object is
// made through the data member "name".
// For details, see the HepRep home page at http://heprep.freehep.org
//  
// @author M.Frailis 
// @author R.Giannitrapani
// @author J.Perl
// Class Description - End:

  
  class G4AttDef{

  public:
    G4AttDef(const G4String& name,
	     const G4String& desc,
	     const G4String& category,
	     const G4String& extra,
	     const G4String& valueType):
      m_name(name),m_desc(desc),
      m_category(category),
      m_extra(extra),m_valueType(valueType){};

    G4AttDef(){};
    
    const G4String& GetName()const{return m_name;};
    const G4String& GetDesc()const{return m_desc;};
    const G4String& GetCategory()const{return m_category;};
    const G4String& GetExtra()const{return m_extra;};
    const G4String& GetValueType()const{return m_valueType;};

    void SetName(const G4String& name){m_name = name;};
    void SetDesc(const G4String& desc){m_desc = desc;};
    void SetCategory(const G4String& cat){m_category = cat;};
    void SetExtra(const G4String& extra){m_extra = extra;};
    void SetValueType(const G4String& type){m_valueType = type;};

  private:
    /// The name of the attribute
    G4String m_name;
    /// A short description of the attribute
    G4String m_desc;
    /// The category (Draw, Physics, PickAction, Association) 
    G4String m_category;
    /// Some extra property of the attribute (units etc)
    G4String m_extra;
    /// The type of the value of the attribute (int, double, string)
    G4String m_valueType;
  };
#endif //G4ATTDEF_H
