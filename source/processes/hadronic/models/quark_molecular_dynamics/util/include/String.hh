#ifndef __STRING__
#define __STRING__

#include "globals.hh"
#include "array.hh"
#include "g4std/strstream"

class G4std::ostream;
class G4std::istream;

class String 
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const String& s);
  friend G4std::istream& operator>>(G4std::istream& i,String& s);
  friend bool operator==(const String& s1,const String& s2);
  friend bool operator==(const String& s1,const char* s2);
  friend bool operator==(const char* s1,const String& s2);
  friend bool operator!=(const String& s1,const char* s2);
  friend bool operator!=(const char* s1,const String& s2);
  friend String operator+(const String& s1,const String& s2);
  friend String operator+(const char* s1,const String& s2);
  friend String operator+(const String& s1,const char* s2);
  friend int length(const String& s) { return s.len; }
  friend String eatwhite(const String& s);
  char* data;
  G4std::istream* str;
  int len;
  String(int l,int);
public:
  String(); 
  String(const char& c);
  String(char* s);
  String(const String& s);
  ~String();
  operator char*() const { return data; }
  int getLength() const { return len; }
  bool empty() const { return bool (len == 0); }
  String& operator=(const String& s);
  String& operator=(const char*);
  String& operator+=(const String& s);
  String& operator+=(const char& c);
  String& erase(const String&);
  String& erase(int start,int end);

  int findFirst(const char& c) const; 
  int findFirst(const String& s) const;
  String subString(int start) const;
  String subString(int start,int end) const;
  
  Array<String> divide(String = ",");
  void putToStream(G4std::istream&) const;
};

#endif
