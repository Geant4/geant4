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
// $Id: SAXStateProcess.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
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

