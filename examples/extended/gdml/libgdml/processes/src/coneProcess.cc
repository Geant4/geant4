#include "SolidTypeProcess.hh"

#include "cone.hh"

class coneProcess : public SolidTypeProcess
{
public:
  coneProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }

  virtual ~coneProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();

    cone* cone_element = new cone;
    
    cone_element->set_rmin1( attrs.getValue( "rmin1" ) );
    cone_element->set_rmax1( attrs.getValue( "rmax1" ) );
    cone_element->set_rmin2( attrs.getValue( "rmin2" ) );
    cone_element->set_rmax2( attrs.getValue( "rmax2" ) );
    cone_element->set_z( attrs.getValue( "z" ) );
    cone_element->set_startphi( attrs.getValue( "startphi" ) );
    cone_element->set_deltaphi( attrs.getValue( "deltaphi" ) );

    m_obj = cone_element;
    *obj  = cone_element;
    
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
    static std::string tag = "cone";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(coneProcess)
