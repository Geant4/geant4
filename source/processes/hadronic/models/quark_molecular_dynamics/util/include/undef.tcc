#include "String.hh"

template<class t>
istream& operator>>(istream& in,undef<t>& x)
{
  char c;
  in >> c;
  if ( c == '?' )
    x.valid = false;
  else {
    in.putback(c);
    if ( !(in >> x.val || in.eof()) )
      throw "Error reading in undef<t>!";
    x.valid = true;
  }
  return in;
}

