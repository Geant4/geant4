#include "SolidTypeProcess.hh"

#include "para.hh"

class paraProcess : public SolidTypeProcess
{
public:
  paraProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }

  virtual ~paraProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();

    para* para_element = new para;
    
    para_element->set_x( attrs.getValue( "x" ) );
    para_element->set_y( attrs.getValue( "y" ) );
    para_element->set_z( attrs.getValue( "z" ) );
    para_element->set_alpha( attrs.getValue( "alpha" ) );
    para_element->set_theta( attrs.getValue( "theta" ) );
    para_element->set_phi( attrs.getValue( "phi" ) );

    m_obj = para_element;
    *obj  = para_element;
    
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
    static std::string tag = "para";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(paraProcess)
