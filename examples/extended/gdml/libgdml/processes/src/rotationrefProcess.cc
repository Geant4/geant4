#include "ReferenceTypeProcess.hh"

#include "BooleanSolidType.hh"

#include <cstdlib>
#include <iostream>

class rotationrefProcess : public ReferenceTypeProcess
{
public:
  rotationrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~rotationrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    BooleanSolidType::rotationref* co = new BooleanSolidType::rotationref;
    
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
    static std::string tag = "rotationref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(rotationrefProcess)

