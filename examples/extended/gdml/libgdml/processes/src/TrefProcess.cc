#include "ReferenceTypeProcess.hh"

#include "TorTref.hh"

class TrefProcess : public ReferenceTypeProcess
{
public:
  TrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~TrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    *obj = new Tref;
    m_obj = *obj;
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "Tref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(TrefProcess)

