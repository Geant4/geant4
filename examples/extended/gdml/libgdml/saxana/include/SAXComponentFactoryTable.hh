#ifndef SAX_COMPONENT_FACTORY_TABLE_H
#define SAX_COMPONENT_FACTORY_TABLE_H 1

#include "SAXComponentObject.hh"

#include <vector>

class SAXComponentFactoryBase;

class SAXComponentFactoryTable {
public:
  typedef std::vector<SAXComponentFactoryBase*> Components;
  
  const SAXComponentFactoryTable::Components* GetComponents( SAXComponentObject::EType type );
    
  void Register( SAXComponentFactoryBase* c );
  
  static SAXComponentFactoryTable* GetInstance();

  ~SAXComponentFactoryTable();
  
protected:
  SAXComponentFactoryTable();
  SAXComponentFactoryTable( const SAXComponentFactoryTable& );

private:
  SAXComponentFactoryTable::Components* fProcesses;
  SAXComponentFactoryTable::Components* fActions;
  SAXComponentFactoryTable::Components* fSubscribers;
};

#endif // SAX_COMPONENT_FACTORY_TABLE_H

