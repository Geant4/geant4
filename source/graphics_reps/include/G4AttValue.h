#ifndef G4ATTVALUE_H
#define G4ATTVALUE_H
#include <string>

/** 
 * @class G4AttValue
 *
 * @brief This class represents a HepRep-style Attribute Value.
 * G4AttValues can be attached to a Trajectory, Trajectory Point or Sensitive
 * Detector Hit.  These attributes are then made available to the end user
 * in interactive visualization systems (such as WIRED).
 * The G4AttValue is further defined in a G4AttDef object.
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

  class G4AttValue {
    
  public:
    G4AttValue(std::string name, std::string value, 
             std::string showLabel): 
      m_name(name),m_value(value),
      m_showLabel(showLabel){};
    G4AttValue():
      m_name(""),m_value(""),
      m_showLabel(""){};
    
    virtual std::string getName(){return m_name;};
    virtual std::string getValue(){return m_value;};
    virtual std::string getShowLabel(){return m_showLabel;};

    virtual void setName(std::string name){m_name = name;};
    virtual void setValue(std::string val){m_value = val;};
    virtual void setShowLabel(std::string lab){m_showLabel = lab;};

  private:
    /// The name of the attribute 
    std::string m_name;
    /// The value of the attribute
    std::string m_value;
    /// The bitmap for the label display
    std::string m_showLabel;
  };
#endif //G4ATTVALUE_H
