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
// $Id: ProcessingContext.hh,v 1.2 2002-06-03 12:09:32 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef PROCESSING_CONTEXT_H
#define PROCESSING_CONTEXT_H 1

class StateStack;
class SAXEvent;
class SAXEventGun;
class ProcessingConfigurator;
class SAXObject;

class ProcessingContext
{
public:
  virtual const StateStack*             GetStack() const = 0;
  virtual const SAXEvent*               GetLastEvent() const = 0;
  virtual const SAXEventGun*            GetSAXEventGun() const = 0; 
  virtual void                          SetSAXEventGun( const SAXEventGun* gun ) = 0;
  virtual const ProcessingConfigurator* GetConfig() const = 0;
  virtual SAXObject**                   GetTopObject() const = 0;
};

#endif // PROCESSING_CONTEXT_H

