#include "String.hh"
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <ctype.h>

typedef Array<String> AString;

ostream& operator<<(ostream& o,const String& s)
{
  return (o << s.data);
}

istream& operator>>(istream& in,String& s)
{
  char p[8000];
  in >> p;
  s = p;
  s.len = strlen(p);
  return in;
}

bool operator==(const String& s1,const String& s2)
{
  return !strcmp(s1.data,s2.data);
}

bool operator==(const String& s1,const char* s2)
{
  return !strcmp(s1.data,s2);
}

bool operator==(const char* s2,const String& s1)
{
  return !strcmp(s1.data,s2);
}

bool operator!=(const String& s1,const char* s2)
{
  return strcmp(s1.data,s2);
}

bool operator!=(const char* s1,const String& s2)
{
  return strcmp(s1,s2.data);
}

String operator+(const String& s1,const String& s2)
{
  String s(length(s1)+length(s2),0);
  strcpy(s.data,s1.data);
  strcat(s.data,s2.data);
  return s;
}

String operator+(const char* s1,const String& s2)
{
  String s(strlen(s1)+length(s2),0);
  strcpy(s.data,s1);
  strcat(s.data,s2.data);
  return s;
}

String operator+(const String& s1,const char* s2)
{
  String s(length(s1)+strlen(s2),0);
  strcpy(s.data,s1.data);
  strcat(s.data,s2);
  return s;
}

String eatwhite(const String& s)
{
  String t;
  for (int i=0; i<s.len; i++) 
    if ( !isspace(s[i]) )
      t += s[i];
  return t;
}

String::String() : data((char*)malloc(sizeof(char))),len(0) { data[0] = '\0';}

String::String(int l,int) : data((char*)malloc((l+1)*sizeof(char))),len(l) { data[0] = 0;}

String::String(const char& c) : data((char*)malloc(2*sizeof(char))),len(1) { data[0] = c; data[1] = 0; }

String::String(char* s) : data(strdup(s)),len(strlen(s)) {}

String::String(const String& s) : data(strdup(s.data)),len(s.len) {}

String::~String() { free(data); }

String& String::operator=(const String& s) 
{
  if ( *this == s ) 
    return *this;
  free(data);
  data = strdup(s.data);
  len = s.len; 
  return *this;
}

String& String::operator=(const char* s) 
{
  free(data);
  data = strdup(s);
  len = strlen(s); 
  return *this;
}

String& String::operator+=(const String& s) 
{
  *this = *this + s;
  return *this;
}

String& String::operator+=(const char& c) 
{
  *this = *this + String(c);
  return *this;
}

int String::findFirst(const char& c) const
{
  char* p = strchr(data,c);
  if ( p ) 
    return int(p-data);
  else
    return -1;
}

int String::findFirst(const String& s) const
{
  char* p = strstr(data,s);
  if ( p ) 
    return int(p-data);
  else
    return -1;
}

String String::subString(int start) const
{
  if ( start<len )
    return String(data+start);
  else
    return String();
}

String String::subString(int start,int end) const
{
  if ( start<len && end<len && start>=0 && end >= 0 ) {
    String s(end-start+1,0);
    strncpy(s.data,data+start,s.len);
    s[s.getLength()] = '\0';
    return s;
  }
  else
    return String();
}

String& String::erase(int start,int end)
{
  if ( start<len && end<len && start>=0 && end >= 0 ) {
    char* s = (char*)malloc((len-(end-start))*sizeof(char));
    strncpy(s,data,start);
    strncpy(s+start,data+end+1,len-end-1);
    s[len = len-(end-start)-1] = '\0';
    free(data);
    data = s;
  }
  return *this;
}

String& String::erase(const String& s)
{
  int p = findFirst(s);
  if (p>=0) {
    erase(p,p+s.getLength()-1);
  }
  return *this;
}

void String::putToStream(istream& in) const
{
  in.putback(' ');
  for (int i=len; i; i--) 
    in.putback(data[i-1]);
}

Array<String> String::divide(String div) 
{
  String s = *this;
  Array<String> List;
  int pos,l = div.getLength();
  do {
    pos = s.findFirst(div);
    if ( pos >= 0 ) {
      List.insert(eatwhite(s.subString(0,pos-1)));
      s.erase(0,pos+l-1);
    }
  }
  while ( pos >= 0 );
  List.insert(eatwhite(s));
  return List;
}

