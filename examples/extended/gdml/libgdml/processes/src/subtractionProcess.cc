#include "BooleanSolidTypeProcess.hh"

#include "boolean_subtraction.hh"

class subtractionProcess : public BooleanSolidTypeProcess
{
public:
  subtractionProcess( const ProcessingContext* context = 0 )
  : BooleanSolidTypeProcess( context ) {
  }

  virtual ~subtractionProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    std::cout << "PROCESS::START OF TAG  : " << name << std::endl;

    SAXObject** obj = Context()->GetTopObject();

    boolean_subtraction* subtraction_lement = new boolean_subtraction;

    m_obj = subtraction_lement;
    *obj  = subtraction_lement;
    
    BooleanSolidTypeProcess::StartElement( name, attrs );
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    BooleanSolidTypeProcess::EndElement( name );
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    std::cout << "PROCESS::subtraction NOTIFIED AFTER THE TAG: " << name << std::endl;
    BooleanSolidTypeProcess::StackPopNotify( name );
  }

  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "subtraction";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(subtractionProcess)
