#ifndef REFERENCETYPEPROCESS_H
#define REFERENCETYPEPROCESS_H 1

#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXObjectHandle.hh"
#include "SAXComponentFactory.hh"

#include "ReferenceType.hh"

#include <cstdlib>
#include <iostream>

class ReferenceTypeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  ReferenceTypeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ), m_obj( 0 ) {
  }
  
  virtual ~ReferenceTypeProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {    
    //std::cout << "ReferenceTypeProcess PROCESS::START OF TAG  : " << name << std::endl;
    
    std::string ref  = attrs.getValue( "ref" );
    ReferenceType* refobj = dynamic_cast<ReferenceType*>(m_obj);
    refobj->set_ref( ref );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    //std::cout << "ReferenceTypeProcess PROCESS::END OF TAG  : " << name << std::endl;
    try {
      ReferenceType* saxobj = dynamic_cast<ReferenceType*>( m_obj );
      
      if( saxobj != 0 ) {
        //std::cout << "PROCESS END OF TAG:: "
        //          << name << " ref: " << saxobj->get_ref()
        //          << " looks OK" << std::endl;
      } else  {
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

#endif // REFERENCETYPEPROCESS_H
