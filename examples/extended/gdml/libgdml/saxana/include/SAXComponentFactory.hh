#ifndef SAX_COMPONENT_FACTORY_H
#define SAX_COMPONENT_FACTORY_H 1

#include "SAXComponentFactoryBase.hh"
#include "SAXComponentFactoryTable.hh"

template <class TYPE>
class SAXComponentFactory : virtual public SAXComponentFactoryBase
{
public:
  SAXComponentFactory( SAXComponentObject::EType type )
  : SAXComponentFactoryBase(), fType( type )
  {
    SAXComponentFactoryTable::GetInstance()->Register( this );
  }
  
  ~SAXComponentFactory()
  {
  }
  virtual SAXComponentObject* Build() const;
  virtual SAXComponentObject::EType Type() const;

private:
  SAXComponentObject::EType fType;
};

template <class TYPE>
SAXComponentObject* SAXComponentFactory<TYPE>::Build() const
{
  SAXComponentObject* obj = new TYPE();
  return obj;
}

template <class TYPE>
SAXComponentObject::EType SAXComponentFactory<TYPE>::Type() const
{
  return fType;
}

#define DECLARE_PROCESS_FACTORY(x) \
static const SAXComponentFactory<x> s_##x##Factory( SAXComponentObject::eProcess ); \
const SAXComponentFactoryBase& x##FactoryRef = s_##x##Factory;

#define DECLARE_SUBSCRIBER_FACTORY(x) \
static const SAXComponentFactory<x> s_##x##Factory( SAXComponentObject::eSubscriber ); \
const SAXComponentFactoryBase& x##FactoryRef = s_##x##Factory;

#define LOAD_COMPONENT(x)         extern const SAXComponentFactoryBase& x##FactoryRef; x##FactoryRef.Load();

#endif // SAX_COMPONENT_FACTORY_H

