//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: SAXEventGun.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_EVENT_GUN_H
#define SAX_EVENT_GUN_H 1

#include <xercesc/sax/HandlerBase.hpp>

#include "StatusCode.hh"

class SAXParser;

class SAXProcessor;
class ProcessingConfigurator;

class SAXEventGun : virtual public HandlerBase
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

  virtual void processingInstruction( const   XMLCh* const  target, const XMLCh* const  data );
  // Receive notification of a processing instruction.

  virtual void resetDocument();
  // Reset the Docuemnt object on its reuse

  virtual void setDocumentLocator( const Locator* const locator );
  // Receive a Locator object for document events.

  virtual void startDocument();
  // Receive notification of the beginning of the document.

  virtual void startElement( const   XMLCh* const    name,       AttributeList&  attributes  );
  // Receive notification of the start of an element.

  virtual InputSource* resolveEntity( const XMLCh* const  publicId, const XMLCh* const systemId );
  // Resolve an external entity.

  virtual void error(const SAXParseException& exception);
  // Receive notification of a recoverable parser error.

  virtual void fatalError(const SAXParseException& exception);
  // Report a fatal XML parsing error.

  virtual void warning(const SAXParseException& exception);
  // Receive notification of a parser warning.

  virtual void resetErrors();
  // Reset the Error handler object on its reuse

  virtual void notationDecl( const XMLCh* const name
                           , const XMLCh* const publicId
                           , const XMLCh* const    systemId );
  // Receive notification of a notation declaration.

  virtual void resetDocType();
  // Reset the DTD object on its reuse

  virtual void unparsedEntityDecl( const XMLCh* const name
                                 , const XMLCh* const publicId
                                 , const XMLCh* const systemId
                                 , const XMLCh* const notationName );
  // Receive notification of an unparsed entity declaration.

  SAXEventGun( SAXProcessor* target );
  ~SAXEventGun();
  
  const SAXProcessor*    Target() const;
  void SetTarget( SAXProcessor* target );
  
  void Configure( ProcessingConfigurator* config );

private:
  SAXProcessor*           fTarget;

protected:
  SAXParser*              fParser;
  ProcessingConfigurator* fConfig;
};

#endif // SAX_EVENT_GUN_H

