#ifndef QUANTITYTYPEPROCESS_H
#define QUANTITYTYPEPROCESS_H 1

#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXObjectHandle.hh"
#include "SAXComponentFactory.hh"

#include "QuantityType.hh"

#include <cstdlib>
#include <iostream>

class QuantityTypeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  QuantityTypeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ) {
  }
  
  virtual ~QuantityTypeProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {    
    QuantityType* qobj = dynamic_cast<QuantityType*>( m_obj );
    
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    std::string qunit  = attrs.getValue( "unit" );
    std::string qtype  = attrs.getValue( "type" );
    std::string qvalue = attrs.getValue( "value" );
    
    qobj->set_unit( qunit );
    qobj->set_type( qtype );
    qobj->set_value( qvalue );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try {
      QuantityType* saxobj = dynamic_cast<QuantityType*>( m_obj );
      
      if( saxobj != 0 )
      {
      //  std::cout << "PROCESS END OF TAG:: " << name
      //            << saxobj->get_type() << ": "
      //            << saxobj->get_value() << "["
      //            << saxobj->get_unit() << "]"                  
      //            << " looks OK" << std::endl;
      }
      else
      {
        std::cerr << "PROCESS END OF TAG:: GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG " << name << " ERROR: "
                << " Cannot cast properly the data object!" << std::endl;
    }
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name )
  {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name )
  {
  }
  
protected:
  SAXObject* m_obj;
};

#endif // QUANTITYTYPEPROCESS_H
