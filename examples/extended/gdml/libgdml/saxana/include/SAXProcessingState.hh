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
// $Id: SAXProcessingState.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_PROCESING_STATE_H
#define SAX_PROCESING_STATE_H 1

#include "ProcessingState.hh"
#include "SAXObject.hh"

class SAXStateProcess;

class SAXProcessingState : virtual public ProcessingState
{
public:
  //typedef RCObjectHandle<SAXStateProcess> Process;
  
  //SAXProcessingState( const SAXObjectHandle& obj, const SAXProcessingState::Process& process )
  SAXProcessingState( SAXObject** obj, SAXStateProcess* process )
  : fObject( obj ), fProcess( process )
  {
  }
  
  ~SAXProcessingState()
  {
  }
  
  virtual ProcessingState::EType Type() const
  {
    return ProcessingState::eState;
  }
  
  //const SAXObjectHandle& GetObject() const
  SAXObject** GetObjectRef() const
  {
    return fObject;
  }
  
  //const SAXObjectHandle& GetObject() const
  //void SetObject( const SAXObject::Ref objPtrRef )
  //{
  //  fObject = objPtrRef;
  //}
  
  //const SAXProcessingState::Process& GetProcess() const
  SAXStateProcess* GetProcess() const
  {
    return fProcess;
  }

private:
  //SAXObjectHandle              fObject;
  //SAXProcessingState::Process  fProcess;
  SAXObject**       fObject;
  SAXStateProcess*  fProcess;
};

#endif // SAX_PROCESING_STATE_H

