#ifndef SAX_BUSCRIBER_POOL_H
#define SAX_BUSCRIBER_POOL_H 1

#include "RCObjectHandle.hh"
#include "SAXSubscriber.hh"

#include <vector>
#include <map>

class SAXSubscriberPool
{
public:
  typedef RCObjectHandle<SAXSubscriber>                     Subscriber;
  typedef std::vector<Subscriber>                           Subscribers;
  typedef std::map<SAXSubscriber::Subscription,Subscribers> Pool;
  
  SAXSubscriberPool();
  ~SAXSubscriberPool();

  void AddSubscriber( const SAXSubscriber::Subscription& type, Subscriber client );
  void RemoveSubscriber( const SAXSubscriber::Subscription& type, Subscriber client );
  const Subscribers* GetSubscribers( const SAXSubscriber::Subscription& subscription );
  void Initialize();
  void Reset();
  
private:
  SAXSubscriberPool::Pool fPool;
};

#endif // SAX_BUSCRIBER_POOL_H


