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
// $Id: GDMLLinkResolver.cc,v 1.2 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/XMLPScanToken.hpp>
#include <xercesc/parsers/SAXParser.hpp>

#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <xercesc/sax/SAXParseException.hpp>

#include "GDMLLinkResolver.hh"

#include "SAXProcessor.hh"
#include "SAXEvents.hh"
#include "ASCIIAttributeList.hh"
#include "ProcessingConfigurator.hh"


//#include <strstream>

StatusCode GDMLLinkResolver::Run()
{
  StatusCode sc;

  fParser = new SAXParser();
  
  fParser->setDocumentHandler( this );
  fParser->setErrorHandler( this );
  fParser->setValidationScheme( SAXParser::Val_Auto );

  const std::string& uri = fConfig->URI();
  const char* xmlfile = uri.c_str();

  try
  {
    // Create a progressive scan token
    XMLPScanToken token;

    if (!fParser->parseFirst(xmlfile, token))
    {
      std::cerr << "scanFirst() failed\n" << std::endl;
      sc = StatusCode::eError;
    }
    else
    {
      //
      //  We started ok, so lets call scanNext() until we find what we want
      //  or hit the end.
      //
      bool gotMore = true;
      while (gotMore && !Done())
          gotMore = fParser->parseNext(token);

      //
      //  Reset the parser. In this simple progrma, since we just exit
      //  now, its not technically required. But, in programs which
      //  would remain open, you should reset after a progressive parse
      //  in case you broke out before the end of the file. This insures
      //  that all opened files, sockets, etc... are closed.
      //
      fParser->parseReset(token);
    }
  }
  catch (const XMLException& toCatch)
  {
    const char* msg = XMLString::transcode( toCatch.getMessage() );
    std::cerr << "\nAn error occured: '" << xmlfile << "'\n"
              << "Exception message is: \n"
              << msg
              << "\n" << endl;
    if( msg != 0 )
    {
      delete [] msg;
    }
    sc = StatusCode::eError;
  }
    
  return sc;
}

void GDMLLinkResolver::characters( const   XMLCh* const    chars
                               , const unsigned int    length )
// Receive notification of #PCDATA characters in element content.
{
  if( fFound )
  {
    SAXEventGun::characters( chars, length );
  }
}

void GDMLLinkResolver::endDocument()
// Receive notification of the end of the document.
{
  if( !fFound && !fIsDone )
  {
    std::cerr << "Invalid GDML link!" << std::endl;
    std::cerr << "Link was: " << fConfig->URI() << " " << fConfig->Type() << " " << fConfig->SetupName() << std::endl;
  }
}

void GDMLLinkResolver::endElement( const XMLCh* const name )
// Receive notification of the end of an element.
{
  const char* aName = XMLString::transcode( name );
  
  if( fFound && fConfig->SetupName() == aName )
  {
    SAXEventGun::endElement( name );
    fIsDone = true;
  }
}

void GDMLLinkResolver::ignorableWhitespace( const   XMLCh* const    chars
                                        , const unsigned int    length )
// Receive notification of ignorable whitespace in element content.
{
  if( fFound )
  {
    SAXEventGun::ignorableWhitespace( chars, length );
  }
}

void GDMLLinkResolver::startDocument()
// Receive notification of the beginning of the document.
{
  fFound  = false;
  fIsDone = false;
}

void GDMLLinkResolver::startElement( const   XMLCh* const    name
                                 , AttributeList&  attributes )
// Receive notification of the start of an element.
{
  const char* cType = XMLString::transcode( name );

  if( fConfig->Type() == cType )
  {
    XMLCh* xName  = XMLString::transcode( fConfig->SetupName().c_str() );
    const XMLCh* xValue = attributes.getValue( xName );
    const char* aName = XMLString::transcode( xValue );
    
    if( fConfig->SetupName() == aName )
    {
      fFound = true;
      SAXEventGun::startElement( name, attributes );
    }
    
    if( xName != 0 )
    {
      delete [] xName;
    }

    if( aName != 0 )
    {
      delete [] aName;
    }
  }
  
  if( cType != 0 )
  {
    delete [] cType;
  }
}

GDMLLinkResolver::GDMLLinkResolver( SAXProcessor* target )
: SAXEventGun( target ), fFound( false ), fIsDone( false )
{
}

GDMLLinkResolver::~GDMLLinkResolver()
{
}



