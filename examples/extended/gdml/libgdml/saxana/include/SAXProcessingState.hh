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

