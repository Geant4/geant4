#ifndef __COLLISIONS__
#define __COLLISIONS__

#include "globals.hh"
#include <list>
#include <vector>
#include "reactionChannels.hh"
#include "ParticleBase.hh"
#include "ParticleType.hh"
#include "Tree.hh"
#include "EventHandling.hh"

class CollisionTab;

class CollisionType : public Knot<CollisionType>
{
  friend class Knot<CollisionType>;
  friend class CollisionTab;
  friend class ParticleType;
  enum { MAX_TRIES = 20 };
  std::vector<ParticleType*> incoming;
  decayMode* minimalMass;
  CollisionType(const std::vector<ParticleBase*>& P);
  CollisionType(const ParticleBase& P);
  CollisionType(const ParticleType& P);
  CollisionType();
public:
  CollisionType(std::istream& in);
  double Crossection(double s) const;
  double checkProducts(const CollisionType& x,bool&) const;
  virtual double isEqualTo(const CollisionType& x) const;
  void perform(const std::vector<ParticleBase*>&,selection = ALL,bool force = false) const;
  static CollisionType* checkCollision(const std::vector<ParticleBase*>&);
  std::vector<decayMode*> channels;
  void print(std::ostream& o) const;
private:
  decayMode& chooseMode(double Emax,selection = ALL,bool = false) const;
  static double FindDecomposition(int,const CollisionType&,const std::vector<ParticleBase*>&);
  //  static bool checkProducts;
};

class CollisionTab 
{
public:
  friend bool operator==(const CollisionTab& x,const CollisionTab& y);
  double time;
  std::vector<ParticleBase*> incoming;
  static std::vector<CollisionTab*> Root;
  const CollisionType& coll;
  selection select;
  CollisionTab(double t,const CollisionType&,const std::vector<ParticleBase*>&,
	       selection);
  CollisionTab(double t,const CollisionType&,const ParticleBase&,
	       selection); 
public:
  typedef std::vector<CollisionTab*>::iterator Entry;
  //  CollisionTab() : time(0),coll((const CollisionType&)(*Knot<CollisionType>::Root)) {}
  ~CollisionTab() {}
  operator double() { return time; }
  void perform();


  static Entry find(double t);
  static Entry exists(const ParticleBase*);
  static Entry next() { return Root.begin(); }
  static void addEntry(double t,const std::vector<ParticleBase*>& in,selection = ALL) ;
  static void addEntry(double t,const ParticleBase& in,selection = ALL) ;
  static void addEntry(double t,const CollisionType& T) ;
  static Entry  remove(Entry);
  static void perform(double time,bool force = false);
  static void erase();
};

#endif
