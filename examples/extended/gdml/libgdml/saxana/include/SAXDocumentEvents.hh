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
// $Id: SAXDocumentEvents.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_DOCUMENT_EVENTS_H
#define SAX_DOCUMENT_EVENTS_H 1

#include "SAXEvent.hh"

class SAXEventStartDocument : virtual public SAXEvent
{
public:
  SAXEventStartDocument() {}
  ~SAXEventStartDocument() {}
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eStartDocument;
  }
};

class SAXEventEndDocument : virtual public SAXEvent
{
public:
  SAXEventEndDocument() {}
  ~SAXEventEndDocument() {}
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eEndDocument;
  }
};

#endif // SAX_DOCUMENT_EVENTS_H

