#include "ReferenceTypeProcess.hh"

#include "VolumeType.hh"

#include <cstdlib>
#include <iostream>

class solidrefProcess : public ReferenceTypeProcess
{
public:
  solidrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~solidrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    
    VolumeType::solidref* so = new VolumeType::solidref;
    
    *obj  = so;
    m_obj = so;
    
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    ReferenceTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "solidref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(solidrefProcess)

