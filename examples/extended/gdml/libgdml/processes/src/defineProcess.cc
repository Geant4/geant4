#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXComponentFactory.hh"

#include "define.hh"

#include <cstdlib>
#include <iostream>

class defineProcess : virtual public SAXStateProcess, virtual public SAXComponentObject
{
public:
  defineProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context )
  {
    m_tag = "define";
  }

  virtual ~defineProcess()
  {
  }

  virtual const SAXComponentObject* Build() const
  {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {
    std::cout << "PROCESS::START OF TAG  : " << name << std::endl;

    SAXObject** obj = Context()->GetTopObject();
    m_define = new define;
    *obj = m_define;
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try
    {
      SAXObject** obj = Context()->GetTopObject();
      defineType* saxobj = dynamic_cast<defineType*>( *obj );

      if( saxobj != 0 )
      {
        std::cout << "PROCESS END OF TAG:: define" << std::endl;
      }
      else
      {
        std::cout << "PROCESS END OF TAG:: define:: GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG:: define:: BAD CAST exception caught...\a" << std::endl;
    }
  }

  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name )
  {
    // Irrelevant in this process
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name )
  {
    // Get the defineType interface
    defineType* define_element = dynamic_cast<defineType*>( m_define );

    // We can expect the 5 inner elements: constant, quantity, expression, position, rotation
    // Whenever one of these appears we have to add it as a choice content group into define element
    SAXObject** so = Context()->GetTopObject();
    if( name == "constant" || name == "quantity" || name == "expression" || name == "position" || name == "rotation" ) {
      // Add the element into choice content group of define element
      define_element->add_to_content( name, *so );
    } else {
      // This should not happen, as it's not fatal problem just warn about it
      std::cerr << "Got unexpected element \"" << name << "\"\a" << std::endl;
    }
  }

  // The name of the state this object will process
  virtual const std::string& State() const
  {
    return m_tag;
  }

private:
  std::string m_tag;
  define*     m_define;
};

DECLARE_PROCESS_FACTORY(defineProcess)

