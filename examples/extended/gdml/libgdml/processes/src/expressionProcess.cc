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

class expressionProcess : virtual public SAXStateProcess, virtual public SAXComponentObject
{
public:
  expressionProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ) {
  }
  
  virtual ~expressionProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    m_expression_value = "(";
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    std::string ename  = attrs.getValue( "name" );

    SAXObject** obj = Context()->GetTopObject();
    *obj = new define::expression;
    define::expression* expression_element = dynamic_cast<define::expression*>(*obj);
    expression_element->set_name( ename );
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    m_expression_value += ")";
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try
    {
      SAXObject** obj = Context()->GetTopObject();
      define::expression* saxobj = dynamic_cast<define::expression*>( *obj );
      saxobj->set_text( m_expression_value );
      
      if( saxobj != 0 )
      {
        //std::cout << "PROCESS END OF TAG:: expression " << saxobj->get_name()
        //          << " has value:: "                    << saxobj->get_text() << std::endl;
      }
      else
      {
        std::cerr << "PROCESS END OF TAG:: expression:: GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG " << name << " ERROR: "
                << " Cannot cast properly the data object!" << std::endl;
    }
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& text )
  {
    m_expression_value += text;
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name )
  {
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string m_tag = "expression";
    return m_tag;
  }
  
private:
  std::string m_expression_value;
};

DECLARE_PROCESS_FACTORY(expressionProcess)

