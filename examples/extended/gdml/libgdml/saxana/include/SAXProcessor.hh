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
// $Id: SAXProcessor.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_PROCESSOR_H
#define SAX_PROCESSOR_H 1

#include "ProcessingContext.hh"
#include "StatusCode.hh"

class SAXEvent;
class SAXEventGun;
class StateStack;
class StateProcessMap;
class SAXSubscriberPool;
class ProcessingConfigurator;
class SAXObject;

class SAXProcessor : virtual public ProcessingContext
{
public:
  SAXProcessor();
  ~SAXProcessor();
  
  // Mandatory interface from ProcessingContext
  virtual const StateStack*             GetStack() const;
  virtual const SAXEvent*               GetLastEvent() const;
  virtual const SAXEventGun*            GetSAXEventGun() const; 
  virtual void                          SetSAXEventGun( const SAXEventGun* gun );
  virtual const ProcessingConfigurator* GetConfig() const;
  virtual SAXObject**                   GetTopObject() const;
  
  StatusCode Initialize();
  StatusCode Configure( ProcessingConfigurator* config );
  StatusCode Run();
  StatusCode Finalize();
  
  void ProcessEvent( const SAXEvent* const event );
  
private:
  StateProcessMap*        fMap;
  SAXSubscriberPool*      fPool;
  StateStack*             fStack;
  StateStack*             fNotifyStack;
  ProcessingConfigurator* fConfig;
  const SAXEvent*         fCurrentEvent;
  const SAXEventGun*      fCurrentGun;
  bool                    fIgnoring;
};

#endif // SAX_PROCESSOR_H

