#include "SolidTypeProcess.hh"

#include "trap.hh"

class trapProcess : public SolidTypeProcess
{
public:
  trapProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }

  virtual ~trapProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();

    trap* trap_element = new trap;
    
    trap_element->set_x1( attrs.getValue( "x1" ) );
    trap_element->set_x2( attrs.getValue( "x2" ) );
    trap_element->set_x3( attrs.getValue( "x3" ) );
    trap_element->set_x4( attrs.getValue( "x4" ) );
    trap_element->set_y1( attrs.getValue( "y1" ) );
    trap_element->set_y2( attrs.getValue( "y2" ) );
    trap_element->set_z( attrs.getValue( "z" ) );
    trap_element->set_theta( attrs.getValue( "theta" ) );
    trap_element->set_phi( attrs.getValue( "phi" ) );
    trap_element->set_alpha1( attrs.getValue( "alpha1" ) );
    trap_element->set_alpha2( attrs.getValue( "alpha2" ) );

    m_obj = trap_element;
    *obj  = trap_element;
    
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
    static std::string tag = "trap";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(trapProcess)
