#ifndef SAX_STATE_PROCESS_H
#define SAX_STATE_PROCESS_H 1

#include <string>

#include "ASCIIAttributeList.hh"
#include "SAXComponentObject.hh"

class ProcessingContext;

class SAXStateProcess : virtual public SAXComponentObject
{
public:
  virtual const SAXComponentObject* Build() const
  {
    return 0;
  }
  virtual SAXComponentObject::EType Type() const
  {
    return SAXComponentObject::eProcess;
  }

public:
  SAXStateProcess( const ProcessingContext* context )
  : fContext( context )
  {
  }
  
  virtual ~SAXStateProcess()
  {
    fContext = 0;
  }
  
  const ProcessingContext* Context() const
  {
    return fContext;
  }
  
  void SetContext( const ProcessingContext* context )
  {
    fContext = context;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) = 0;
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) = 0;
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) = 0;
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) = 0;
  // The name of the state this object will process
  virtual const std::string& State() const = 0;

private:
  const ProcessingContext* fContext;
};

#endif // SAX_STATE_PROCESS_H

