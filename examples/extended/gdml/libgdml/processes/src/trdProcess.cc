#include "SolidTypeProcess.hh"

#include "trd.hh"

class trdProcess : public SolidTypeProcess
{
public:
  trdProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }

  virtual ~trdProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();

    trd* trd_element = new trd;
    
    trd_element->set_x1( attrs.getValue( "x1" ) );
    trd_element->set_y1( attrs.getValue( "y1" ) );
    trd_element->set_x2( attrs.getValue( "x2" ) );
    trd_element->set_y2( attrs.getValue( "y2" ) );
    trd_element->set_z( attrs.getValue( "z" ) );

    m_obj = trd_element;
    *obj  = trd_element;
    
    SolidTypeProcess::StartElement( name, attrs );
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    SolidTypeProcess::EndElement( name );
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    //std::cout << "PROCESS::" << name << " NOTIFIED AFTER THE TAG: " << name << std::endl;
    SolidTypeProcess::StackPopNotify( name );
  }

  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "trd";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(trdProcess)
