#include "ReferenceTypeProcess.hh"

#include "BooleanSolidType.hh"

#include <cstdlib>
#include <iostream>

class secondProcess : public ReferenceTypeProcess
{
public:
  secondProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~secondProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    BooleanSolidType::second* co = new BooleanSolidType::second;
    
    *obj = co;
    m_obj = co;
    
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    ReferenceTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "second";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(secondProcess)

