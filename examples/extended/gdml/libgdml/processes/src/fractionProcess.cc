#include "ReferenceTypeProcess.hh"

#include "fraction.hh"

#include <cstdlib>
#include <iostream>

class fractionProcess : public ReferenceTypeProcess
{
public:
  fractionProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~fractionProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    fraction* fr = new fraction;
    
    std::string fn = attrs.getValue( "n" );
    fr->set_n( fn );
    
    *obj = fr;
    m_obj = fr;
    
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
    static std::string tag = "fraction";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(fractionProcess)

