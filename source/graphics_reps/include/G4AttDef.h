#ifndef G4ATTDEF_H
#define G4ATTDEF_H
#include <string>

/** 
 * @class G4AttDef
 *
 * @brief This class represents a HepRep-style Attribute Definition.
 * The G4AttDef is used to define new kinds of attributes that can
 * then have values set for a Trajectory, Trajectory Point or Sensitive
 * Detector Hit.  These attributes are then made available to the end user
 * in an interactive visualization system (such as WIRED).
 * Values are set by creating G4AttValue objects and attaching them to the
 * relevant Trajectory, Trajectory Point or Sensitive Detector Hit.
 * The association between the G4AttValue and the G4AttDef object is
 * made through the data member "name".
 * For details, see the HepRep home page at http://heprep.freehep.org
 *  
 * @author M.Frailis 
 * @author R.Giannitrapani
 * @author J.Perl
 *    
 * \$Header\$ 
 */
  
  class G4AttDef{

  public:
    G4AttDef(std::string name, std::string desc, 
	   std::string category, std::string extra, 
	   std::string valueType):
      m_name(name),m_desc(desc),
      m_category(category),
      m_extra(extra),m_valueType(valueType){};

    G4AttDef():
      m_name(""),m_desc(""),
      m_category(""),
      m_extra(""),m_valueType(""){};
    
    std::string getName(){return m_name;};
    std::string getDesc(){return m_desc;};
    std::string getCategory(){return m_category;};
    std::string getExtra(){return m_extra;};
    std::string getValueType(){return m_valueType;};

    void setName(std::string name){m_name = name;};
    void setDesc(std::string desc){m_desc = desc;};
    void setCategory(std::string cat){m_category = cat;};
    void setExtra(std::string extra){m_extra = extra;};
    void setValueType(std::string type){m_valueType = type;};

  private:
    /// The name of the attribute
    std::string m_name;
    /// A short description of the attribute
    std::string m_desc;
    /// The category (Draw, Physics, PickAction, Association) 
    std::string m_category;
    /// Some extra property of the attribute (units etc)
    std::string m_extra;
    /// The type of the value of the attribute (int, double, string)
    std::string m_valueType;
  };
#endif //G4ATTDEF_H
