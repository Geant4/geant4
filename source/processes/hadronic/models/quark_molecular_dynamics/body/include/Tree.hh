#ifndef __TREE__
#define __TREE__

#ifdef IS_GCC
#pragma interface
#endif

#include "globals.hh"
#include <vector>
#include "String.hh"


template<class t> class Knot;

template<class t>
class Knot 
{
public:
  static Knot<t>* Root;
private:
  t* ancestor;
  std::vector<t*> successors;
  String name;

  void printTree(std::ostream&,int) const;
  void find(const t&,std::vector<t*>& ) const;
  void find(const String&,std::vector<t*>& ) const;
  void findKnot(const t&,t*&) const;
  void findKnot(const String&,t*&) const;
  void flatten(std::vector<t*>&) const;
  void count(int& n) const;
protected:
  Knot() : ancestor(0),name("root") {}
  virtual ~Knot();
  void insert(t& x);
  t& Successor(int i) const { return *successors[i]; }
  int NSuccessors() const { return successors.size(); }
  virtual double isEqualTo(const t& x) const { return 0; }
  virtual void ClassInfo(std::ostream& o) const { }
public:
  std::vector<t*> getList() const;
  void printTree(std::ostream& o) ;
  int countEntries();
  String Name() const { return name; }
  bool hasSuccessors() const { return successors.size() > 0; }
  bool isRoot() const { return !ancestor; }
  t& goUp() const { return *ancestor; }
  static std::vector<t*> Find(const t&) ;
  static std::vector<t*> Find(const String&) ;
  static t& FindKnot(const t&) ;
  static t& FindKnot(const String&) ;
  static void Insert(t& x,const String& name_);
};

// -----------------------------------------------------
// implementation from Tree.tcc:

template<class t>
void Knot<t>::printTree(std::ostream& o,int n) const
{
  String s;
  for (int i=0; i<n; i++)
    s += "  ";
  o << Name();
  ClassInfo(o);
  o << G4endl;
  for (int j=0; j<successors.size(); j++) {
    o << s << j+1 << ".";
    successors[j]->printTree(o,n+1);
  }
}

template<class t>
void Knot<t>::find(const t& x,std::vector<t*>& list) const
{
  if ( isEqualTo(x)>=0 || !ancestor ) {
    int l = successors.size();
    if ( l>0 ) {
      for (int i=0; i<successors.size(); i++)
        successors[i]->find(x,list);
      if ( successors.size() == l ) {
        std::vector<t*> new_list;
        flatten(new_list);
        list.insert(list.end(),new_list.begin(),new_list.end());
      }
    }
    else
      list.insert(list.end(),(t*)this);
  }
}

template<class t>
void Knot<t>::find(const String& s,std::vector<t*>& list) const
{
  if ( Name() == s ) 
    flatten(list);
  else
    for (int i=0; i<successors.size(); i++)
      successors[i]->find(s,list);
}

template<class t>
void Knot<t>::findKnot(const t& x,t*& found) const
{
  double y = isEqualTo(x);
  if ( y>=0 || !ancestor ) {
    for (unsigned int i=0; i<successors.size() && !found && y<1; i++)
      successors[i]->findKnot(x,found);
    if ( !found )
      found = (t*)this;
  }
}

template<class t>
void Knot<t>::findKnot(const String& x,t*& found) const
{
  if ( Name() == x ) 
    found = (t*)this;
  else
    for (unsigned int i=0; i<successors.size() && !found; i++)
      successors[i]->findKnot(x,found);
}

template<class t>
void Knot<t>::flatten(std::vector<t*>& list) const
{
  if ( !successors.size() ) 
    list.insert(list.end(),(t*)this);
  else
    for (int i=0; i<successors.size(); i++)
      successors[i]->flatten(list);
}

template<class t>
void Knot<t>::count(int& n) const
{
  ++n;
  for (int i=0; i<successors.size(); i++)
    successors[i]->count(n);
}

template<class t>
Knot<t>::~Knot()
{
  for (unsigned int i=0; i<successors.size(); i++) {
    delete [] successors[i];
  }
}

template<class t>
void Knot<t>::insert(t& x)
{
  unsigned int num = successors.size();
  unsigned int i;
  for (i=0; i<successors.size(); i++) {
    double u = successors[i]->isEqualTo(x);
    if ( u>=0 ) {
      if ( u<1 ) {
        successors[i]->insert(x);
        break;
      }
      else 
        --num;
    }
    else {
      //      break;
    }
  }
  if ( i == successors.size() && !x.ancestor ) {
    if ( num== successors.size() || num == 0 ) {
      successors.insert(successors.end(),&x);
      x.ancestor = (t*)this;
    }
    else {
      for (i=0; i<successors.size(); i++) 
        if ( successors[i]->isEqualTo(x) == 1 ) {
          successors[i]->insert(x);
          break;
        }
    }
  }
}

template<class t>
std::vector<t*> Knot<t>::getList() const
{
  std::vector<t*> list;
  flatten(list);
  return list;
}

template<class t>
void Knot<t>::printTree(std::ostream& o) 
{
  printTree(o,1);
}

template<class t>
int Knot<t>::countEntries()
{
  int n = -1;
  count(n);
  return n;
}

template<class t>
std::vector<t*> Knot<t>::Find(const t& x) 
{
  std::vector<t*> list;
  Root->find(x,list);
  return list;
}

template<class t>
std::vector<t*> Knot<t>::Find(const String& x) 
{
  std::vector<t*> list;
  Root->find(x,list);
  return list;
}

template<class t>
t& Knot<t>::FindKnot(const t& x) 
{
  t* found = 0;
  Root->findKnot(x,found);
  if ( !found ) 
    throw "Knot not found...";
  return *found;
}

template<class t>
t& Knot<t>::FindKnot(const String& x) 
{
  t* found = 0;
  Root->findKnot(x,found);
  if ( !found ) 
    throw String("'"+x+"' not found...");
  return *found;
}

template<class t>
void Knot<t>::Insert(t& x,const String& name_)
{
  x.name = name_;
  Root->insert(x);
}

// ------------------------------------------------

#endif
