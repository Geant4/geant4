#ifndef __TREE__
#define __TREE__

#ifdef IS_GCC
#pragma interface
#endif

#include "globals.hh"
#include <vector.h>
#include "String.hh"

template<class t>
class Knot 
{
public:
  static Knot<t>* Root;
private:
  t* ancestor;
  vector<t*> successors;
  String name;

  void printTree(ostream&,int) const;
  void find(const t&,vector<t*>& ) const;
  void find(const String&,vector<t*>& ) const;
  void findKnot(const t&,t*&) const;
  void findKnot(const String&,t*&) const;
  void flatten(vector<t*>&) const;
  void count(int& n) const;
protected:
  Knot() : ancestor(0),name("root") {}
  virtual ~Knot();
  void insert(t& x);
  t& Successor(int i) const { return *successors[i]; }
  int NSuccessors() const { return successors.size(); }
  virtual double isEqualTo(const t& x) const { return 0; }
  virtual void ClassInfo(ostream& o) const { }
public:
  vector<t*> getList() const;
  void printTree(ostream& o) ;
  int countEntries();
  String Name() const { return name; }
  bool hasSuccessors() const { return successors.size() > 0; }
  bool isRoot() const { return !ancestor; }
  t& goUp() const { return *ancestor; }
  static vector<t*> Find(const t&) ;
  static vector<t*> Find(const String&) ;
  static t& FindKnot(const t&) ;
  static t& FindKnot(const String&) ;
  static void Insert(t& x,const String& name_);
};

#ifndef IS_GCC
#include "Tree.tcc"
#endif

#endif
