#include <iostream.h>
#include <strstream.h>
#include "StreamBuffer.hh"
#include "array.hh"

ostream& operator<<(ostream& o,const StreamBuffer& x) { return o << x.buffer; }

StreamBuffer::StreamBuffer() : in(0) {}

StreamBuffer::StreamBuffer(istream& i) : in(0) 
{
  String s;
  while ( i >> s ) 
    *this << s;
}

StreamBuffer& StreamBuffer::operator<<(const String& s) 
{ 
  buffer += " "+s; 
  return *this; 
}
  
StreamBuffer::operator istream&() 
{ 
  return ( in ) ? *in : *(in = new istrstream((char*)buffer,length(buffer))); 
}

void StreamBuffer::reset() 
{ 
  if (in) 
    delete in; 
  in = new istrstream((char*)buffer,length(buffer)); 
}

void StreamBuffer::createStream(int n,char** args) 
{
  for (int i=0; i<n; i++) {
    String s(args[i]);
    if ( s.findFirst("--")>=0 ) {
      s.erase("--");
      if ( s.empty() )
	s = args[++i];
      Array<String> X = s.divide("=");
      for (int j=0; j<X.size(); j++)
	(*this) << X[j];
      (*this) << "\n";
    }
  }
}

