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

#ifndef G4ATTDEF_HH
#define G4ATTDEF_HH

#include "globals.hh"
#include "G4TypeKey.hh"
#include <map>

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

    // G4Typekey based constructor
    G4AttDef(const G4String& name,
             const G4String& desc,
             const G4String& category,
             const G4String& extra,
             const G4TypeKey& typeKey):
      m_name(name),m_desc(desc),
      m_category(category),
      m_extra(extra),m_valueType("Null"), 
      m_typeKey(typeKey)
    {};

    G4AttDef()= default;
    virtual ~G4AttDef()= default;
    
    const G4String& GetName()const{return m_name;};
    const G4String& GetDesc()const{return m_desc;};
    const G4String& GetCategory()const{return m_category;};
    const G4String& GetExtra()const{return m_extra;};
    const G4String& GetValueType()const{return m_valueType;};
    const G4TypeKey& GetTypeKey()const{return m_typeKey;};

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
    /// The category (Draw, Physics, PickAction, Association, etc.) 
    G4String m_category;
    /// Some extra property of the attribute (units, etc.)
    G4String m_extra;
    /// The type of the value of the attribute (int, double, vector, etc.)
    G4String m_valueType;
    // Type key
    G4TypeKey m_typeKey;

  };

// Deprecated.  It is not a good idea to output a pointer since failure to
// include this prototype will not cause a compilation error - it will merely
// cause your code to use a default function that outputs the pointer as an
// address.
std::ostream& operator<<
(std::ostream& os, const std::map<G4String,G4AttDef>* definitions);

// Use this instead.
std::ostream& operator<<
(std::ostream& os, const std::map<G4String,G4AttDef>& definitions);

#endif //G4ATTDEF_H
