#ifndef __UNDEF__
#define __UNDEF__
#ifdef IS_GCC
#pragma interface
#endif

#include "g4std/iostream"

enum joker { DUMMY = 32767 };

template<class t>
class undef
{
  friend G4std::istream& operator>>(G4std::istream& in,undef<t>&);
  t val,def;
  bool valid;
public:
  undef() : def(DUMMY),valid(false) {}
  undef(const t& x) : val(x),def(DUMMY),valid(true) {}
  void setDefault(const t& x) { def = x; }
  bool isValid() const { return valid; }
  operator t() const { return valid ? val : def; }
};

#ifndef IS_GCC
#include "undef.tcc"
#endif

#endif
