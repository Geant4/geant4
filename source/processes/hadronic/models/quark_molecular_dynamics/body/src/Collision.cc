#include <algo.h>
#include "globals.hh"
#include "Collision.hh"
#include "Permutations.hh"
#include "ParticleProperties.hh"

//bool CollisionType::checkProducts = false;

Knot<CollisionType>* Knot<CollisionType>::Root = new CollisionType;

bool operator==(const CollisionTab& x,const CollisionTab& y) 
{ 
  return x.time == y.time; 
}

vector<CollisionTab*> CollisionTab::Root;

CollisionTab::Entry CollisionTab::find(double t) {
  Entry X = Root.begin();
  while ( X<Root.end() && (*X)->time < t ) 
    ++X;
  return X;
}

void CollisionTab::perform() {
  coll.perform(incoming,select);
}

void CollisionTab::erase() {
  while ( !Root.empty() ) {
    CollisionTab* X = Root.front();
    delete X;
    Root.erase(Root.begin());
  }
}

void CollisionTab::perform(double t,bool force) {
  vector<CollisionTab*> toPerform;
  vector<CollisionTab*>::iterator W;
  for ( W = Root.begin(); W!=Root.end() && (*W)->time <= t; W++) {
    toPerform.insert(toPerform.end(),*W);
  }
  Root.erase(Root.begin(),W);
  for (vector<CollisionTab*>::iterator X = toPerform.begin(); X!=toPerform.end(); X++) {
    CollisionTab* Y = *X;
      if ( Y->incoming.size() ) {
	try {
	  Y->coll.perform(Y->incoming,Y->select,force);
	  for (unsigned int i=0; i<Y->incoming.size(); i++) {
	    delete Y->incoming[i];
	  }
	}
	catch ( char* s) {
	  G4cerr << "Could not perform decay: " << s << G4endl;
	}
	delete Y;
      }
      else
	throw "COLLTAB-ERROR!!";
  }
  toPerform.erase(toPerform.begin(),toPerform.end());
}

CollisionTab::Entry CollisionTab::exists(const ParticleBase* in)
{
  Entry X = Root.begin();
  while ( X != Root.end() ) {
    bool found = true;
    for (unsigned int i=0; i<(*X)->incoming.size(); i++)
      found = found && ( in == (*X)->incoming[i] );
    if ( found ) return X;
    ++X;
  }
  return X;
}

void CollisionTab::addEntry(double t,const ParticleBase& in,selection which) {
  //  G4cerr << "CollTab: " << in.Name() << G4endl;
  Entry y = exists(&in);
  Entry x = find(t);
  if ( x<y || y == Root.end() ) {
    if ( y < Root.end() )
      Root.erase(y);
    CollisionType& coll = Knot<CollisionType>::FindKnot(CollisionType(in));
    CollisionType* T = &coll;
    if ( which == DECOMPOSITION ) {
      bool decFound = false;
      while ( !T->isRoot() && !decFound ) {
	for (int i=0; i<T->channels.size(); i++) 
	  if ( T->channels[i]->isDecomposition() ) {
	    decFound = true;
	    break;
	  }
	if ( !decFound ) 
	  T = &(T->goUp());
      }
      if ( !decFound ) 
	throw "No collision found...";
    }
    Root.insert(x,new CollisionTab(t,*T,in,which));
  }
}

void CollisionTab::addEntry(double t,const vector<ParticleBase*>& in,selection which) {
  Entry x = find(t);
  vector<Entry> existing;
  unsigned int i;
  for ( i=0; i<in.size(); i++) {
    Entry y = exists(in[i]);
    if ( y<x ) 
      break;
    existing.insert(existing.end(),y);
  }
  if ( i<in.size() ) {
    for (int j=existing.size(); j; j--) 
      Root.erase(existing[j-1]);
    CollisionType& coll = Knot<CollisionType>::FindKnot(CollisionType(in));
    Root.insert(x,new CollisionTab(t,coll,in,which));
  }
}


CollisionTab::Entry CollisionTab::remove(CollisionTab::Entry x) {
  Entry y = x++;
  delete (*y);
  Root.erase(y);
  return x;
}

CollisionTab::CollisionTab(double t,const CollisionType& type,const vector<ParticleBase*>& in,selection which) 
    : time(t),incoming(in),coll(type),select(which) {}

CollisionTab::CollisionTab(double t,const CollisionType& type,const ParticleBase& in,selection which) 
    : time(t),incoming(vector<ParticleBase*>(1)),coll(type),select(which) 
{ 
  incoming[0] = &(ParticleBase&)in; 
}

CollisionType::CollisionType() : minimalMass(0) {}

CollisionType::CollisionType(const vector<ParticleBase*>& P) : minimalMass(0)
{
  for (unsigned int i=0; i<P.size(); i++)
    incoming.insert(incoming.end(),&((ParticleType&)(P[i]->getType())));
}

CollisionType::CollisionType(const ParticleBase& P) : minimalMass(0)
{
  incoming.insert(incoming.end(),&(ParticleType&)(P.getType()));
}

CollisionType::CollisionType(const ParticleType& P) : minimalMass(0)
{
  incoming.insert(incoming.end(),&(ParticleType&)P);
}

CollisionType::CollisionType(G4std::istream& in) : minimalMass(0)
{
  String Name,CollName;
  try {
    while ( in >> Name && Name != "->" ) {
      ParticleType& X = Knot<ParticleType>::FindKnot(Name);
      incoming.insert(incoming.end(),&X);
      if ( length(CollName) == 0 )
	CollName = Name;
      else
	CollName += ","+Name;
    }
    while ( in && Name == "->" ) {
      decayMode* y = NEW decayMode(in,incoming);
      channels.insert(channels.end(),y);
      if ((!minimalMass || y->SumMass<minimalMass->SumMass) && y->partialCrossection(1)>0.0 )
	minimalMass = y;
      in >> Name;
    }
    Name.putToStream(in);
    Insert(*this,CollName);
  }
  catch ( ... ) {
    throw String("CollisionTable: unknown particle type \""+Name+"\"");
  }
}

void CollisionType::print(G4std::ostream& o) const {
  for (int i=0; i<incoming.size(); i++) { o << i << ". " << *incoming[i] << G4endl; }
  o << "channels: ";
  for (int j=0; j<channels.size(); j++) { o << j << ". " << *channels[j] << G4endl; }
}

double CollisionType::Crossection(double s) const 
{
  double S = 0;
  for (int i=0; i<channels.size(); i++)
    S += channels[i]->partialCrossection(s);
  return S;
}

CollisionType* CollisionType::checkCollision(const vector<ParticleBase*>& L)
{
  double rp = length(L[0]->Momentum()-L[1]->Momentum());
  if ( rp<1.0 ) return 0;
  double d = length(L[0]->Coordinates()-L[1]->Coordinates());
  if ( d>1.0 ) return 0;
  try {
    CollisionType& T = Knot<CollisionType>::FindKnot(CollisionType(L));
    double s = sqrt(sqr(L[0]->E()+L[1]->E())-square(L[0]->Momentum()+L[1]->Momentum()));
    if ( T.Crossection(s) > M_PI*d*d )
      return &T;
    else
      return 0;
  }
  catch ( ... ) { return 0; }
}

double CollisionType::FindDecomposition(int C,const CollisionType& X,const vector<ParticleBase*>& x)
{
  CollisionType y;
  y.channels.insert(y.channels.end(),new decayMode(C,x));
  bool decFound = false;
  CollisionType* T = (CollisionType*)&X;
  double l;
  while ( !T->isRoot() && (l=T->checkProducts(y,decFound))<0 && !decFound )
    T = &(T->goUp());
  if ( T->isRoot() || (decFound && l<0) )
    return -1;
  else
    return fabs(T->channels[int(l*T->channels.size())]->partialCrossection(1.0));
}

double CollisionType::checkProducts(const CollisionType& x,bool& decFound) const
{
  for (int l=0; l<channels.size(); l++) {
    decFound = decFound || channels[l]->isDecomposition();
    vector< vector<isotop*> > Perm = Permutations(channels[l]->products);
    for (unsigned int i=0; i<Perm.size(); i++) {
      bool y = false,s=true;
      for (unsigned int j=0; j<channels[l]->products.size(); j++) {
	y = (*Perm[i][j] == *x.channels[0]->products[j]);
	if ( !y ) {
	  s = false;
	  break;
	}
      }
      if ( s ) 
	return double(l)/channels.size();
    }
  }
  return -1.0;
}

double CollisionType::isEqualTo(const CollisionType& x) const
{
    if ( incoming.size() != x.incoming.size() )
      return -1;
    vector< vector<ParticleType*> > Perm = Permutations(incoming);
    int best = -1;
    double val = -1;
    for (unsigned int i=0; i<Perm.size(); i++) {
      double s = 0,y = 0;
      for (unsigned int j=0; j<incoming.size(); j++) {
	y = Perm[i][j]->isEqualTo(*x.incoming[j]);
	if ( y == 1.0 && *Perm[i][j] != *x.incoming[j] ) 
	  y = -1; 
	if ( y<0 ) 
	  break;
	else 
	  s += y;
      }
      // find possible "best match"
      if ( y>=0 && s > val ) {
	best = i;
	val = s;
      }
    }
    return val;
}

decayMode& CollisionType::chooseMode(double Emax,selection which,bool force) const
{
  if ( which == LOWEST )
    if ( minimalMass->SumMass > Emax && !force )
      throw "No valid decay mode...";
    else
      return *minimalMass;
  double* y = new double[channels.size()+1];
  decayMode** M = new decayMode*[channels.size()+1];
  y[0] = 0;
  int k=0;
  decayMode* act;
  for (unsigned int i=0; i<channels.size(); i++) {
   act = channels[i];
   if ( ( act->SumMass<Emax || force ) && ( (!which && !(act->select & DECOMPOSITION)) || which & act->select) ) {
     //	 (act->Fraction > 0 && which>=0 
     //	  || act->Fraction < 0 && which & act->select ) ) {
      ++k;
      y[k] = y[k-1]+fabs(act->partialCrossection(Emax));
      M[k] = channels[i];
    }
  }
  if ( k>0 ) {
    double r = rand_gen::Random()*y[k];
    int kmin=1;
    while ( kmin<=k && r>y[kmin] ) 
      kmin++;
    act = M[kmin];
    delete [] y;
    delete [] M;
  }
  else {
    delete [] y;
    delete [] M;
    throw "DECAY ERROR: Mass too low!";
  }
  return *act;
}

void CollisionType::perform(const vector<ParticleBase*>& p,selection which,bool force) const
{
  if ( p.empty() ) 
    return;
  int valid = MAX_TRIES;
  G4cerr << *p[0] << G4endl;
  double sqrt_s = p[0]->Mass();
  decayMode* y;
  do {
    try {
      if ( valid )
	y = &chooseMode(sqrt_s,which,force);
      else 
	y = &chooseMode(sqrt_s,selection(LOWEST+which),force);
      Vektor3 beta=(-1/p[0]->E())*p[0]->Momentum();
      G4cerr << Name() << " --> (" << which << ") ";
      y->performDecay(p,sqrt_s,beta,p[0]->Coordinates(),force);
      valid = -1;
    }
    catch ( char* s ) {
      if ( valid) {
	G4cerr << valid << ". DECAY ERROR: " 
	     << Name() << ": " << s << G4endl;
	--valid;
      }
      else {
	G4cerr << "FATAL DECAY ERROR! :";
	G4cerr << Name() << ": " << s << G4endl;
	throw;
      }
    }
  }
  while ( valid>=0 );
}
