#ifndef GDML_LINK_RESOLVER_H
#define GDML_LINK_RESOLVER_H 1

#include "SAXEventGun.hh"

class GDMLLinkResolver : virtual public SAXEventGun
{
public:
  StatusCode Run();

public:
  virtual void characters( const   XMLCh* const    chars, const unsigned int    length );
  // Receive notification of the #PCDATA characters in element content.

  virtual void endDocument();
  // Receive notification of the end of the document.

  virtual void endElement(const XMLCh* const name);
  // Receive notification of the end of an element.

  virtual void ignorableWhitespace( const   XMLCh* const  chars, const unsigned int  length );
  // Receive notification of ignorable whitespace in element content.

  virtual void startDocument();
  // Receive notification of the beginning of the document.

  virtual void startElement( const   XMLCh* const    name,       AttributeList&  attributes  );
  // Receive notification of the start of an element.

  GDMLLinkResolver( SAXProcessor* target );
  ~GDMLLinkResolver();
  
  bool Done()
  {
    return fIsDone;
  }
  
private:
  bool fFound;
  bool fIsDone;
};

#endif // GDML_LINK_RESOLVER_H

