#include "ReferenceTypeProcess.hh"

#include "DorDref.hh"

#include <cstdlib>
#include <iostream>

class DrefProcess : public ReferenceTypeProcess
{
public:
  DrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~DrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    *obj = new Dref;
    m_obj = *obj;
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "Dref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(DrefProcess)

