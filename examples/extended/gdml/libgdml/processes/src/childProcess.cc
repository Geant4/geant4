#include "SinglePlacementTypeProcess.hh"

#include "child.hh"

class childProcess : public SinglePlacementTypeProcess
{
public:
  childProcess( const ProcessingContext* context = 0 )
  : SinglePlacementTypeProcess( context ) {
  }

  virtual ~childProcess() {
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    std::cout << "PROCESS::START OF TAG  : " << name << std::endl;

    SAXObject** obj = Context()->GetTopObject();

    child* child_lement = new child;

    m_obj = child_lement;
    *obj  = child_lement;
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    std::cout << "PROCESS::child NOTIFIED AFTER THE TAG: " << name << std::endl;
    SinglePlacementTypeProcess::StackPopNotify( name );
  }

  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "child";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(childProcess)
