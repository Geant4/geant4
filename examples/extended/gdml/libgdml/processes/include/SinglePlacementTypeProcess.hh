#ifndef SINGLEPLACEMENTTYPEPROCESS_H
#define SINGLEPLACEMENTTYPEPROCESS_H 1

#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXComponentFactory.hh"

#include "SinglePlacementType.hh"

class SinglePlacementTypeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  SinglePlacementTypeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ) {
  }
  
  virtual ~SinglePlacementTypeProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
  }

  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    SAXObject** so = Context()->GetTopObject();
    SinglePlacementType* sptobj = dynamic_cast<SinglePlacementType*>( m_obj );
    sptobj->add_content( name, *so );
  }

protected:
  SAXObject* m_obj;
};

#endif // SINGLEPLACEMENTTYPEPROCESS_H
