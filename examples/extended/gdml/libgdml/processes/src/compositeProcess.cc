#include "ReferenceTypeProcess.hh"

#include "composite.hh"

#include <cstdlib>
#include <iostream>

class compositeProcess : public ReferenceTypeProcess
{
public:
  compositeProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~compositeProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    composite* co = new composite;
    
    std::string cn = attrs.getValue( "n" );
    co->set_n( cn );
    
    *obj = co;
    m_obj = co;
    
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    ReferenceTypeProcess::EndElement( name );
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "composite";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(compositeProcess)

