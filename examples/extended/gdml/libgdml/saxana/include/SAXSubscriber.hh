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
// $Id: SAXSubscriber.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SAX_SUBSCRIBER_H
#define SAX_SUBSCRIBER_H 1

#include "SAXComponentObject.hh"
#include "SAXObject.hh"

#include <string>
#include <vector>

class SAXSubscriber : virtual public SAXComponentObject
{
public:
  typedef std::string               Subscription;
  typedef std::vector<Subscription> SubscriptionSet;
  
public:
  virtual const SAXComponentObject* Build() const
  {
    return 0;
  }
  virtual SAXComponentObject::EType Type() const
  {
    return SAXComponentObject::eSubscriber;
  }

public:
  SAXSubscriber()
  : fObject( 0 ), fSet( 0 )
  {
    fSet = new SubscriptionSet();
  }

  virtual ~SAXSubscriber()
  {
    if( fSet != 0 )
    {
      delete fSet;
      fSet = 0;
    }
  }

  SAXObject* GetObject()
  {
    return fObject;
  }

  SAXObject* GetObject() const
  {
    return fObject;
  }

  void SetObject( SAXObject* object )
  {
    fObject = object;
  }

  void Subscribe( const Subscription s )
  {
    fSet->push_back( s );
  }
  void Subscribe( const char* s )
  {
    Subscription ss = s;
    Subscribe( ss );
  }

  const SubscriptionSet* GetSubscriptions() const
  {
    return fSet;
  }

  bool IsSubscribedTo( const Subscription& s ) const
  {
    size_t howmany = fSet->size();
    if( howmany > 1 ) {
      for( size_t idx = 0; idx < howmany; idx++ ) {
        if( (*fSet)[idx] == s ) {
          return true;
        }
      }
    } else if( howmany == 1 ) {
      return( (*fSet)[0] == s );
    } else {
      return false;
    }
    return false;
  }

  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) = 0;

protected:
  SAXObject*       fObject;

private:
  SubscriptionSet* fSet;
};

#endif // SAX_SUBSCRIBER_H


