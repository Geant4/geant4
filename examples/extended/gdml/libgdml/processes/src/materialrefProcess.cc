#include "ReferenceTypeProcess.hh"

#include "VolumeType.hh"

#include <cstdlib>
#include <iostream>

class materialrefProcess : public ReferenceTypeProcess
{
public:
  materialrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~materialrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    VolumeType::materialref* mo = new VolumeType::materialref;
    
    *obj  = mo;
    m_obj = mo;
    
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    ReferenceTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "materialref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(materialrefProcess)

