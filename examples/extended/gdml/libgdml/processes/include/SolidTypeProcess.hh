#ifndef SOLIDTYPEPROCESS_H
#define SOLIDTYPEPROCESS_H 1

#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXObjectHandle.hh"
#include "SAXComponentFactory.hh"

#include "SolidType.hh"

#include <cstdlib>
#include <iostream>

class SolidTypeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  SolidTypeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ), m_obj( 0 ) {
  }
  
  virtual ~SolidTypeProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    std::string lunit  = attrs.getValue( "lunit" );
    std::string aunit  = attrs.getValue( "aunit" );
    std::string sname   = attrs.getValue( "name" );
    
    SolidType* sobj  = dynamic_cast<SolidType*>( m_obj );
    
    sobj->set_lunit( lunit );
    sobj->set_aunit( aunit );
    sobj->set_name( sname );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
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

#endif // SOLIDTYPEPROCESS_H
