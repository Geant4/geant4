#include "MaterialTypeProcess.hh"

#include "material.hh"

#include <cstdlib>
#include <iostream>



class materialProcess : public MaterialTypeProcess
{
public:
  materialProcess( const ProcessingContext* context = 0 )
  : MaterialTypeProcess( context ), m_sequence( 0 ) {
  }
  
  virtual ~materialProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {    
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    
    // Reset the sequence storage
    m_sequence = 0;
    
    std::string mname = attrs.getValue( "name" );
    std::string mf    = attrs.getValue( "formula" );
    std::string mz    = attrs.getValue( "Z" );
    std::string ms    = attrs.getValue( "state" );
    
    SAXObject** obj = Context()->GetTopObject();

    material* mat = new material;
    
    mat->set_name( mname );
    mat->set_formula( mf );
    mat->set_state( ms );
    
    if( !mz.empty() )
      mat->set_Z( mz );
    
    m_obj = mat;
    *obj  = mat;
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    //std::cout << "PROCESS::END OF TAG  : " << name << " ";
    try
    {
      SAXObject** obj = Context()->GetTopObject();
      material* saxobj = dynamic_cast<material*>( *obj );
      
      std::cout << saxobj->get_name() << std::endl;
      
      if( saxobj != 0 ) {
        std::string tag;
        size_t count = 0;

        ContentChoice* cc = saxobj->get_choice();
        if( cc->content().object->type() == SAXObject::element ) {
          // Must be an atom
          tag = cc->content().tag;
          //std::cout << "Simple definition by atom" << std::endl;
        } else {
          // Must be a sequence of composites or fractions
          //std::cout << "Complex definition by ";
          const ContentSequence* cseq = dynamic_cast<const ContentSequence*>( cc->content().object );
          count = cseq->size();
          std::string tag = cseq->content( 0 ).tag;
          if( tag == "composite" ) {
            // Composition by atoms of elements
            //std::cout << count << " composites" << std::endl;
          } else {
            // Composition by fraction of mass
            //std::cout << count << " fractions" << std::endl;
          }
        }
      
        //std::cout << " looks OK" << std::endl;
      } else {
        std::cerr << "PROCESS END OF TAG:: material GOT ZERO DATA POINTER! " << std::endl;
      }
    }
    catch( ... )
    {
      std::cerr << "PROCESS END OF TAG " << name << " ERROR: "
                << " Cannot cast properly the data object!" << std::endl;
    }
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name )
  {
    //std::cout << "PROCESS::material NOTIFIED AFTER THE TAG: " << name << std::endl;
       
    SAXObject** obj = Context()->GetTopObject();
    material*   mat = dynamic_cast<material*>( m_obj );
    
		if( name == "atom" ) {
      mat->set_choice( name, *obj );
    } else if( name == "fraction" || name == "composite" ) {
      if( m_sequence == 0 ) {
        m_sequence = new ContentSequence;
        mat->set_choice( "sequence", m_sequence );
      }
      ContentGroup::ContentItem ci = { name, *obj };
      m_sequence->add_content( ci );
		}	else if( name == "D" || name == "Dref" ) {
      mat->set_DorDref( name, *obj );
    } else {
      MaterialTypeProcess::StackPopNotify( name );
    }
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "material";
    return tag;
  }
  
private:
  // Sequence of either composites or fractions
  ContentSequence* m_sequence;  
};

DECLARE_PROCESS_FACTORY(materialProcess)

