#include "SAXComponentFactoryBase.hh"
#include "SAXComponentFactoryTable.hh"

static SAXComponentFactoryTable* sComponentsFactoryTable = 0;

// Returns a vector of components according to the given component type
const SAXComponentFactoryTable::Components*
SAXComponentFactoryTable::GetComponents( SAXComponentObject::EType type )
{
  SAXComponentFactoryTable::Components* ret = 0;
  
  switch( type )
  {
    case SAXComponentObject::eProcess   : ret = fProcesses; break;
    case SAXComponentObject::eAction    : ret = fActions; break;
    case SAXComponentObject::eSubscriber: ret = fSubscribers; break;
  };
  
  return ret;
}
  
void SAXComponentFactoryTable::Register( SAXComponentFactoryBase* c )
{
  SAXComponentObject::EType type = c->Type();
  
  switch( type )
  {
    case SAXComponentObject::eProcess   : fProcesses->push_back( c ); break;
    case SAXComponentObject::eAction    : fActions->push_back( c ); break;
    case SAXComponentObject::eSubscriber: fSubscribers->push_back( c ); break;
  };
}

SAXComponentFactoryTable* SAXComponentFactoryTable::GetInstance()
{
  if( 0 == sComponentsFactoryTable )
  {
    sComponentsFactoryTable = new SAXComponentFactoryTable();
  }
  
  return sComponentsFactoryTable;
}

SAXComponentFactoryTable::~SAXComponentFactoryTable()
{
  if( fProcesses != 0 )
  {
    fProcesses->clear();
    delete fProcesses;
    fProcesses = 0;
  }
  if( fActions != 0 )
  {
    fActions->clear();
    delete fActions;
    fActions = 0;
  }
  if( fSubscribers != 0 )
  {
    fSubscribers->clear();
    delete fSubscribers;
    fSubscribers = 0;
  }
}

#include <iostream>

SAXComponentFactoryTable::SAXComponentFactoryTable()
: fProcesses( 0 ), fActions( 0 ), fSubscribers( 0 )
{
  fProcesses   = new SAXComponentFactoryTable::Components();
  fActions     = new SAXComponentFactoryTable::Components();
  fSubscribers = new SAXComponentFactoryTable::Components();
  std::cout << "SAXComponentFactoryTable created" << std::endl;
}


