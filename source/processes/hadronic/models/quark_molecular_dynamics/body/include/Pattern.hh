#ifndef __PATTERN__
#define __PATTERN__
#ifdef IS_GCC
#pragma interface
#endif

#include "globals.hh" 


template<class t> class Pattern;

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


// -----------------------------------------------------
// implementation from pattern.tcc:

template<class t>
Pattern<t>::Pattern(int n) : N(n),field(new t[n]) 
{
  for (int i=0; i<N; i++)
    field[i] = DUMMY;
}

template<class t>
Pattern<t>::Pattern(const Pattern<t>& x) : N(x.N),field(new t[N])
{
  for (int i=0; i<N; i++)
    field[i] = x.field[i];
}

template<class t>
Pattern<t>::~Pattern()
{
  delete [] field;
}

template<class t>
double Pattern<t>::matches(const Pattern<t>& x) const
{
  int count=0;
  for (int i=0; i<N; i++ ) {
    if ( field[i] == x.field[i] && field[i] != DUMMY )
      ++count;
    else if ( field[i] != DUMMY )
      return -1;
  }
  return double(count)/N;
}

template<class t>
int Pattern<t>::joker() const
{
  int count = 0;
  for (int i=0; i<N; i++)
    if ( field[i] == DUMMY )
      ++count;
  return count;
}

template<class t>
Pattern<t>& Pattern<t>::operator=(const Pattern<t>& x)
{
  for (int i=0; i<N; i++)
    field[i] = x.field[i];
  return *this;
}

// -----------------------------------------------------


#endif
