#ifndef __STREAMBUFFER__
#define __STREAMBUFFER__

#include "String.hh"

class G4std::istream;
class G4std::ostream;

class StreamBuffer
{
  friend G4std::ostream& operator<<(G4std::ostream&,const StreamBuffer&);
  G4std::istream* in;
  String buffer;
public:
  StreamBuffer();
  StreamBuffer(G4std::istream&);
  StreamBuffer& operator<<(const String& s);
  operator G4std::istream&();
  void reset();
  void createStream(int n,char** args);
};


#endif
