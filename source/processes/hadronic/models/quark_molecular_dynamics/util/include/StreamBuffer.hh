#ifndef __STREAMBUFFER__
#define __STREAMBUFFER__

#include "String.hh"
#include <iostream>

class StreamBuffer
{
  friend std::ostream& operator<<(std::ostream&,const StreamBuffer&);
  std::istream* in;
  String buffer;
public:
  StreamBuffer();
  StreamBuffer(std::istream&);
  StreamBuffer& operator<<(const String& s);
  operator std::istream&();
  void reset();
  void createStream(int n,char** args);
};


#endif
