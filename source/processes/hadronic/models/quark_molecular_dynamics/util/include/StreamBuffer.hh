#ifndef __STREAMBUFFER__
#define __STREAMBUFFER__

#include "String.hh"

class istream;
class ostream;

class StreamBuffer
{
  friend ostream& operator<<(ostream&,const StreamBuffer&);
  istream* in;
  String buffer;
public:
  StreamBuffer();
  StreamBuffer(istream&);
  StreamBuffer& operator<<(const String& s);
  operator istream&();
  void reset();
  void createStream(int n,char** args);
};


#endif
