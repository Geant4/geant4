#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXComponentFactory.hh"

#include "setup.hh"

#include <cstdlib>
#include <iostream>

class setupProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  setupProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ), m_ignoreTheRest( false ) {
  }

  virtual ~setupProcess() {
  }

  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;

    SAXObject** obj = Context()->GetTopObject();
    
    m_setup = new setup;
    m_setup->set_name( attrs.getValue( "name" ) );
    m_setup->set_version( attrs.getValue( "version" ) );
    
    if( Context()->GetConfig()->SetupName().empty() ) {
      std::cerr << "PROCESS:: " << name
                << " WARNING:: No setup name given, taking the first found..." << std::endl;
      m_setup->set_match( true );
    } else {
      if( m_setup->get_name() == Context()->GetConfig()->SetupName() ) {
        if( Context()->GetConfig()->SetupVersion().empty() ) {
          std::cerr << "PROCESS:: " << name
                    << " WARNING:: No setup version given, taking the first found..." << std::endl;
          m_setup->set_match( true );
        } else {
          if( m_setup->get_version() == Context()->GetConfig()->SetupVersion() ) {
            m_setup->set_match( true );
          }
        }
      }
    }
    
    *obj = m_setup;
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    //std::cout << "PROCESS::END OF TAG  : " << name << std::endl;
    try {
      SAXObject** obj = Context()->GetTopObject();
      setup* saxobj = dynamic_cast<setup*>( *obj );

      if( saxobj != 0 ) {
//         std::cout << "PROCESS END OF TAG:: " << name;
//         std::cout << " name: "               << m_setup->get_name()
//                   << " version: "            << m_setup->get_version()
//                   << " found: "              << (int)m_setup->get_match() << std::endl;
      } else {
        std::cerr << "PROCESS END OF TAG:: " << name << " GOT ZERO DATA POINTER! " << std::endl;
      }
    } catch( ... ) {
      std::cerr << "PROCESS END OF TAG:: " << name << " BAD CAST exception caught...\a" << std::endl;
    }
  }

  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
    // Irrelevant in this process
  }

  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    // Get the defineType interface

    // We can expect the 5 inner elements: constant, quantity, expression, position, rotation
    // Whenever one of these appears we have to add it as a choice content group into define element
    SAXObject** so = Context()->GetTopObject();
    if( !m_ignoreTheRest ) {
      setup::world* world_element = dynamic_cast<setup::world*>( *so );
      if( world_element != 0 ) {
        if( m_setup->get_match() ) {
          m_setup->set_world( *world_element );
          m_ignoreTheRest = true;
        }
      } else {
        std::cerr << "NOTIFY AFTER TAG:: " << name << " ERROR:: BAD OBJECT TYPE...\a" << std::endl;      
      }
    } else {
      std::cout << "NOTIFY AFTER TAG:: " << name << " Ignoring world volume..." << std::endl;
    }
  }

  // The name of the state this object will process
  virtual const std::string& State() const {
    static std::string m_tag = "setup";
    return m_tag;
  }

private:
  bool       m_ignoreTheRest;
  setup*     m_setup;
};

DECLARE_PROCESS_FACTORY(setupProcess)

