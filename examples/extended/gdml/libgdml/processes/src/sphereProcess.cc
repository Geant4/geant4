#include "SolidTypeProcess.hh"

#include "sphere.hh"

class sphereProcess : public SolidTypeProcess
{
public:
  sphereProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }

  virtual ~sphereProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();

    sphere* sphere_element = new sphere;
    
    sphere_element->set_rmin( attrs.getValue( "rmin" ) );
    sphere_element->set_rmax( attrs.getValue( "rmax" ) );
    sphere_element->set_startphi( attrs.getValue( "startphi" ) );
    sphere_element->set_deltaphi( attrs.getValue( "deltaphi" ) );
    sphere_element->set_starttheta( attrs.getValue( "starttheta" ) );
    sphere_element->set_deltatheta( attrs.getValue( "deltatheta" ) );

    m_obj = sphere_element;
    *obj  = sphere_element;
    
    SolidTypeProcess::StartElement( name, attrs );
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    SolidTypeProcess::EndElement( name );
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    std::cout << "PROCESS::" << name << " NOTIFIED AFTER THE TAG: " << name << std::endl;
    SolidTypeProcess::StackPopNotify( name );
  }

  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "sphere";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(sphereProcess)
