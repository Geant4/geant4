#include "ReferenceTypeProcess.hh"

#include "PorPref.hh"

class PrefProcess : public ReferenceTypeProcess
{
public:
  PrefProcess( const ProcessingContext* context = 0 )
  : ReferenceTypeProcess( context ) {
  }
  
  virtual ~PrefProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    SAXObject** obj = Context()->GetTopObject();
    *obj = new Pref;
    m_obj = *obj;
    ReferenceTypeProcess::StartElement( name, attrs );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string tag = "Pref";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(PrefProcess)

