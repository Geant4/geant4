#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXObjectHandle.hh"
#include "SAXComponentFactory.hh"

#include "define.hh"

#include <cstdlib>
#include <iostream>

class positionProcess : virtual public SAXStateProcess, virtual public SAXComponentObject
{
public:
  positionProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ) {
  }
  
  virtual ~positionProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    std::string pname = attrs.getValue( "name" );
    std::string ptype = attrs.getValue( "type" );
    std::string px    = attrs.getValue( "x" );
    std::string py    = attrs.getValue( "y" );
    std::string pz    = attrs.getValue( "z" );
    std::string punit = attrs.getValue( "unit" );

    SAXObject** obj = Context()->GetTopObject();
    *obj = new define::position;
    define::position* position_element = dynamic_cast<define::position*>(*obj);
    position_element->set_name( pname );
    position_element->set_type( ptype );
    position_element->set_unit( punit );
    position_element->set_x( px );
    position_element->set_y( py );
    position_element->set_z( pz );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try
    {
      SAXObject** obj = Context()->GetTopObject();
      define::position* saxobj = dynamic_cast<define::position*>( *obj );
      
      if( saxobj != 0 )
      {
        //std::cout << "PROCESS END OF TAG:: position " << saxobj->get_name() << " looks OK" << std::endl;
        ;
      }
      else
      {
        std::cerr << "PROCESS END OF TAG:: GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG " << name << " ERROR: "
                << " Cannot cast properly the data object!" << std::endl;
    }
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name )
  {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name )
  {
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string m_tag = "position";
    return m_tag;
  }
};

DECLARE_PROCESS_FACTORY(positionProcess)

