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

class quantityProcess : virtual public SAXStateProcess, virtual public SAXComponentObject
{
public:
  quantityProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ) {
  }
  
  virtual ~quantityProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )  {    
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    std::string qname  = attrs.getValue( "name" );
    std::string qvalue = attrs.getValue( "value" );
    std::string qunit  = attrs.getValue( "unit" );
    std::string qtype  = attrs.getValue( "type" );

    SAXObject** obj = Context()->GetTopObject();
    *obj = new define::quantity;
    define::quantity* quantity_element = dynamic_cast<define::quantity*>(*obj);
    quantity_element->set_name( qname );
    quantity_element->set_value( qvalue );
    quantity_element->set_unit( qunit );
    quantity_element->set_type( qtype );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try
    {
      SAXObject** obj = Context()->GetTopObject();
      define::quantity* saxobj = dynamic_cast<define::quantity*>( *obj );
      
      if( saxobj != 0 )
      {
        //std::cout << "PROCESS END OF TAG:: quantity " << saxobj->get_name() << " looks OK" << std::endl;
      }
      else
      {
        std::cerr << "PROCESS END OF TAG:: quantity:: GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG " << name << " ERROR::"
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
    static std::string tag = "quantity";
    return tag;
  }
};

DECLARE_PROCESS_FACTORY(quantityProcess)

