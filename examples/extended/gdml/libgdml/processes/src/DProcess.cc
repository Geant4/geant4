#include "QuantityTypeProcess.hh"

#include "DorDref.hh"

class DProcess : public QuantityTypeProcess
{
public:
  DProcess( const ProcessingContext* context = 0 )
  : QuantityTypeProcess( context ) {
  }
  
  virtual ~DProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {    

    SAXObject** obj = Context()->GetTopObject();
    *obj = new D;
    m_obj = *obj;

    QuantityTypeProcess::StartElement( name, attrs );    
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    QuantityTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string m_tag = "D";
    return m_tag;
  }
};

DECLARE_PROCESS_FACTORY(DProcess)

