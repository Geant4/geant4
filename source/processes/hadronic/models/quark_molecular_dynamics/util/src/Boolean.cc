/*
Example programs from the book Scientific and Engineering Programming
in C++: An Introduction with Advanced Techniques and Examples,
Addison-Wesley, 1994.
 
                 (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
                 CORPORATION 1994.  ALL RIGHTS RESERVED.

See README file for further details.
*/
#include "Boolean.h"
#include <string.h>
#include "g4std/iostream"

G4std::ostream& operator<<(G4std::ostream& s, Boolean b) {
  return s << (b ? "true" : "false");
}

G4std::istream& operator>>(G4std::istream& s, Boolean& b) {
  const size_t bufsize = 1 + sizeof("false");
  char inbuf[bufsize];
  s.width(int(bufsize));     // Set maximum # of characters to be read
  s >> inbuf;		      // Read character string (also resets width)
  if (strcmp(inbuf, "true") == 0)       b = Boolean::True;
  else if (strcmp(inbuf, "false") == 0) b = Boolean::False;
  else s.clear(G4std::ios::failbit);    // Indicate extraction failed
  return s;
}
