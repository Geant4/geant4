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
// $Id: G4AttDef.hh,v 1.1 2002-10-23 12:08:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4ATTDEF_HH
#define G4ATTDEF_HH

#include <string>

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
    G4AttDef(G4std::string name, G4std::string desc, 
	   G4std::string category, G4std::string extra, 
	   G4std::string valueType):
      m_name(name),m_desc(desc),
      m_category(category),
      m_extra(extra),m_valueType(valueType){};

    G4AttDef():
      m_name(""),m_desc(""),
      m_category(""),
      m_extra(""),m_valueType(""){};
    
    G4std::string getName(){return m_name;};
    G4std::string getDesc(){return m_desc;};
    G4std::string getCategory(){return m_category;};
    G4std::string getExtra(){return m_extra;};
    G4std::string getValueType(){return m_valueType;};

    void setName(G4std::string name){m_name = name;};
    void setDesc(G4std::string desc){m_desc = desc;};
    void setCategory(G4std::string cat){m_category = cat;};
    void setExtra(G4std::string extra){m_extra = extra;};
    void setValueType(G4std::string type){m_valueType = type;};

  private:
    /// The name of the attribute
    G4std::string m_name;
    /// A short description of the attribute
    G4std::string m_desc;
    /// The category (Draw, Physics, PickAction, Association) 
    G4std::string m_category;
    /// Some extra property of the attribute (units etc)
    G4std::string m_extra;
    /// The type of the value of the attribute (int, double, string)
    G4std::string m_valueType;
  };
#endif //G4ATTDEF_H
