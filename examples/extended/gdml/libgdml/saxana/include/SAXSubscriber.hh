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


