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
// $Id: SAXProcessor.cc,v 1.3 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

#include "ProcessingConfigurator.hh"
#include "SAXEvents.hh"
#include "StateStack.hh"
#include "SAXEventGun.hh"
#include "StateProcessMap.hh"
#include "SAXStateProcess.hh"
#include "SAXSubscriber.hh"
#include "SAXSubscriberPool.hh"
#include "SAXObjectHandle.hh"
#include "SAXProcessingState.hh"
#include "ActionProcessingState.hh"

#include "SAXProcessor.hh"

#include <iostream>

SAXProcessor::SAXProcessor()
: fMap( 0 ), fPool( 0 ), fStack( 0 ), fNotifyStack( 0 ),
  fConfig( 0 ), fCurrentEvent( 0 ), fIgnoring( false )
{
  fMap         = new StateProcessMap();
  fPool        = new SAXSubscriberPool();
  fStack       = new StateStack();
  fNotifyStack = new StateStack();
}

SAXProcessor::~SAXProcessor()
{
  if( fMap != 0 )
  {
    delete fMap;
    fMap = 0;
  }  
  if( fPool != 0 )
  {
    delete fPool;
    fPool = 0;
  }  
  if( fStack != 0 )
  {
    delete fStack;
    fStack = 0;
  }
  if( fNotifyStack != 0 )
  {
    delete fNotifyStack;
    fStack = 0;
  }
}

// Mandatory interface from ProcessingContext
const StateStack* SAXProcessor::GetStack() const
{
  return fStack;
}

const SAXEvent* SAXProcessor::GetLastEvent() const
{
  return fCurrentEvent;
}

const SAXEventGun* SAXProcessor::GetSAXEventGun() const
{
  return fCurrentGun;
}
 
void SAXProcessor::SetSAXEventGun( const SAXEventGun* gun )
{
  fCurrentGun = gun;
}

const ProcessingConfigurator* SAXProcessor::GetConfig() const
{
  return fConfig;
}

SAXObject** SAXProcessor::GetTopObject() const
{
  SAXProcessingState* state = dynamic_cast<SAXProcessingState*>( fStack->Top() );
  SAXObject** obj = state->GetObjectRef();
  return obj;
}


StatusCode SAXProcessor::Initialize()
{
  StatusCode sc;
  
  fMap->Initialize();
  fPool->Initialize();
  
  // Initialize the Xerces-C system
  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch( const XMLException& toCatch )
  {
    char* msg = XMLString::transcode( toCatch.getMessage() );
    std::cerr << "Error during initialization! :\n" << msg << std::endl;
    sc = StatusCode::eError;
    if( msg != 0 )
    {
      delete [] msg;
    }
  }
  
  return sc;
}

StatusCode SAXProcessor::Configure( ProcessingConfigurator* config )
{
  fConfig = config;
  return StatusCode::eOk;
}

StatusCode SAXProcessor::Run()
{
  StatusCode sc;
  SAXEventGun saxgun( this );
  
  saxgun.Configure( fConfig );
  sc = saxgun.Run();
  
  return sc;
}

StatusCode SAXProcessor::Finalize()
{
  StatusCode sc;
  
  // And call the termination method
  XMLPlatformUtils::Terminate();
  
  return sc;
}

void SAXProcessor::ProcessEvent( const SAXEvent* const event )
{
  fCurrentEvent = event;
  SAXProcessingState* top = 0;
  SAXProcessingState* notifytop = 0;
  SAXObject** objectRef = 0;
    
  if( !fStack->Empty() )
  {
    //top = dynamic_cast<SAXProcessingState*>( (fStack->Top())() );
    top = dynamic_cast<SAXProcessingState*>( fStack->Top() );
  }
  
  if( !fNotifyStack->Empty() )
  {
    //top = dynamic_cast<SAXProcessingState*>( (fStack->Top())() );
    notifytop = dynamic_cast<SAXProcessingState*>( fNotifyStack->Top() );
  }
  
  if( SAXEvent::eStartElement == event->Type() )
  {
    const SAXEventStartElement* ev = dynamic_cast<const SAXEventStartElement*>(event);

    std::string tagname = ev->Name();
    //std::cout << "SXP::PE:: Got start element event for tag: " << tagname << std::endl;
        
    if( fMap->Check( ev->Name() ) )
    {
      fIgnoring = false;
      
      // We've got the state process in place, let's activate it
      objectRef = new SAXObject*;
      StateProcessMap::Process process = fMap->GetProcess( tagname );
      process->SetContext( this );
      const StateStack::State state   = new SAXProcessingState( objectRef, process() );
      
      // Setup properly the stacks
      if( fStack->Empty() && fNotifyStack->Empty() )
      {
        fStack->Push( state );
	    }
	    else
	    {
        StateStack::State notifystate = top;
	      fNotifyStack->Push( notifystate );
	      fStack->Push( state );
	    }

      process->StartElement( tagname, ev->Attributes() );
    }
    else
    {
      fIgnoring = true;
    }
  }
  else if( SAXEvent::eEndElement == event->Type() )
  {
    const SAXEventEndElement* ev = dynamic_cast<const SAXEventEndElement*>(event);

    std::string name = ev->Name();
    //std::cout << "SXP::PE:: Got end element event for tag: " << name << std::endl;

    if( !fStack->Empty() ) {
        if( ev->Name() == top->GetProcess()->State() ) {
        // Activate state process
        SAXStateProcess* process = top->GetProcess();
        process->EndElement( ev->Name() );

        // Grab a new object from the stack and play with it a bit
        objectRef = top->GetObjectRef();

        // Locate the subcriber(s) and send them the created object
        const SAXSubscriberPool::Subscribers* actors = fPool->GetSubscribers( ev->Name() );

        if( actors != 0 ) {
          // There are some guys waiting out there for a gift
          SAXSubscriberPool::Subscribers::const_iterator subscriberRef;

          for( subscriberRef = actors->begin(); subscriberRef != actors->end(); subscriberRef++ ) {
            // Now we call only subscribers processing only this single element
            if( (*subscriberRef)->GetSubscriptions()->size() == 1 ) {
              //std::cout << "SXP::PE:: Executing subscriber(s) for element: " << ev->Name() << std::endl;
              (*subscriberRef)->Activate( *objectRef );
            }
          }

          delete actors;
        }

        if( !fNotifyStack->Empty() ) {
          // Let the parent process know about its child element
          notifytop->GetProcess()->StackPopNotify( ev->Name() );

          // Now we call only the subscribers who need to know about their children elements
          // This workaround is here due to the single-type tree grammar. This type of tree grammar
          // forbids conflicting non-terminal symbols in its productions and this is here to remove
          // the ambiguity caused by the conflicting non-terminals
          // In other words it means whenever we have to process two XML elements which have the same
          // child element in their content models we notify only the one whose child is now created
          SAXSubscriber::Subscription ParentTag = notifytop->GetProcess()->State();
          const SAXSubscriberPool::Subscribers* actors = fPool->GetSubscribers( ParentTag );
          if( actors != 0 ) {
            SAXSubscriberPool::Subscribers::const_iterator subscriberRef;

            for( subscriberRef = actors->begin(); subscriberRef != actors->end(); subscriberRef++ ) {
              if( (*subscriberRef)->IsSubscribedTo( ev->Name() ) ) {
                //std::cout << "SXP::PE:: Executing PARENT subscriber(s) for: "
                //          << ParentTag << " for element: " << ev->Name() << std::endl;
                (*subscriberRef)->Activate( *objectRef );
              }
            }
            
            delete actors;
          }

          // Pop this guy from notify stack
          fNotifyStack->Pop();
        }

        if( objectRef != 0 ) {
          delete objectRef;
        }

        fStack->Pop();
        //std::cout << "--------------------------------------------------------------" << std::endl;
      }
    }
  }
  else
  {
    switch( event->Type() )
    {
      case SAXEvent::eCharacters :
      {
        if( !fIgnoring ) {
          const SAXEventCharacters* ev = dynamic_cast<const SAXEventCharacters*>(event);
          //std::cout << "SXP::PE:: Got characters event >>" << ev->Data() << "<<" << std::endl;
          if( !fStack->Empty() ) {
            top->GetProcess()->Characters( ev->Data() );
          }
        }
      }
                                    break;
      case SAXEvent::ePI         :
      {
      }
                                    break;
      case SAXEvent::eWarning    :
      case SAXEvent::eError      :
      case SAXEvent::eFatalError :
      {
        const SAXErrorEventBase* ev = dynamic_cast<const SAXErrorEventBase*>(event);
        std::cerr << ev->Message() << std::endl;
      }
                                    break;
      default                    :
                                    break;
    };
  }
}



