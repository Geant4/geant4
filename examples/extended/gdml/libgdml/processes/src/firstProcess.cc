#include "ReferenceTypeProcess.hh"

#include "BooleanSolidType.hh"

#include <cstdlib>
#include <iostream>

class firstProcess : public ReferenceTypeProcess
{
public:
  firstProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~firstProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    BooleanSolidType::first* co = new BooleanSolidType::first;
    
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
    static std::string tag = "first";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(firstProcess)

