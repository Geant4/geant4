#include "ReferenceTypeProcess.hh"

#include "SinglePlacementType.hh"

#include <cstdlib>
#include <iostream>

class volumerefProcess : public ReferenceTypeProcess
{
public:
  volumerefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~volumerefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    SinglePlacementType::volumeref* vo = new SinglePlacementType::volumeref;
    
    *obj  = vo;
    m_obj = vo;
    
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    ReferenceTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "volumeref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(volumerefProcess)

