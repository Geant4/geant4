#ifndef __PATTERN__
#define __PATTERN__
#ifdef IS_GCC
#pragma interface
#endif

#include "globals.hh" 

template<class t>
class Pattern
{
  friend bool operator==(const Pattern<t>& x,const Pattern<t>& y) { return x.matches(y) >= 0; }
  enum { DUMMY=32767 };
  int N;
  t* field;
protected:
  Pattern(int n);
  Pattern(const Pattern<t>&);
  virtual ~Pattern();
  t getEntry(int i) const { return field[i]; }
  void setEntry(int i,t x) { field[i] = x; }
public:
  double matches(const Pattern<t>& x) const;
  int joker() const;
  Pattern<t>& operator=(const Pattern<t>& x);
};

#ifndef IS_GCC
#include "Pattern.tcc"
#endif

#endif
