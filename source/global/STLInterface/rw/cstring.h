#ifndef __cstring
#define __cstring

#include <string>

#if defined(G4USE_NAMESPACE)
  using std::string;
#endif

#ifndef G4USE_OLDSTL
#include <cstring>
#else
#include <string.h>
#endif
#include "rw/defs.h"

#if defined(__KCC)
#include <istream>
#define std_string  std::string
#define std_istream std::istream
#define std_ws      std::ws
#else
#define std_string  string
#define std_istream istream
#define std_ws      ws
#endif

#if defined(WIN32)
#define strcasecmp _stricmp
#endif

class RWCString;

class RWCSubString
{
public:

  inline RWCSubString(const RWCSubString&);

  inline RWCSubString& operator=(const char*);         

  inline RWCSubString& operator=(const RWCString&s);
  inline RWCSubString& operator=(const RWCSubString&s);
 
  inline char& operator()(size_t i);
  inline char  operator()(size_t i) const;
  inline char& operator[](size_t i);
  inline char  operator[](size_t i) const;
  size_t        length() const{return extent;}
  size_t        start() const{return mystart;}
  //void          toLower();
  //void          toUpper();

  // For detecting null substrings:
  RWBoolean     isNull() const
    {return extent==0 ? true : false;}

  int operator!() const
    {return extent==0 ? 1 : 0;}

  inline RWBoolean operator==(RWCString) const;
  inline RWBoolean operator==(const char*) const;

  inline RWBoolean operator!=(RWCString) const;
  inline RWBoolean operator!=(const char*) const;

private:

  RWCSubString(RWCString&str,size_t s,size_t e):
    mystring(&str),mystart(s),extent(e){} 

  RWCString*    mystring;     
  size_t        mystart;
  size_t        extent;

  friend class RWCString;

};
 
//#define std 
 
class RWCString : public std_string {
  
public: 

  enum caseCompare { exact, ignoreCase };
  enum stripType { leading, trailing, both };

  inline RWCString ( const char * );

  RWCString ( char c )
    {
      char str[2];
      str[0]=c;
      str[1]='\0';
      std_string::operator=(str);
    }

  RWCString (){}
  RWCString ( const RWCString& s):std_string(s){}
  RWCString( const RWCSubString&s ):std_string(*(s.mystring),s.mystart,s.extent){}
  RWCString ( const std_string& s ):std_string(s){}

  RWCString& operator=(const RWCString&s)
    {
      std_string::operator=(s);
      return *this;
    }
  RWCString& operator=(const std_string&s)
    {
      std_string::operator=(s);
      return *this;
    }
  RWCString& operator=(const char* s)
    {
      std_string::operator=(s);
      return *this;
    }
  virtual ~RWCString (){}
  inline char operator () (size_t) const; 
  char& operator () (size_t i){return std_string::operator[](i);}
  int compareTo(const char*s,caseCompare mode=exact)
    {
      //GB:OSF1-cxx if(mode=exact)
      if(mode==exact)
        return strcmp(c_str(),s);
      else
        return strcasecmp(c_str(),s);
    }
  int compareTo(const RWCString&s,caseCompare mode=exact)
    {
      //GB:OSF1-cxx if(mode=exact)
      if(mode==exact)
        return strcmp(data(),s.data());
      else
        return strcasecmp(data(),s.data());     
    }
  RWCString& operator+=(const RWCSubString&s)
    {
      RWCString tmp(s);
      std_string::operator+=(tmp);
      return *this;
    }
  RWCString& operator+=(const char*s)
    {
      std_string::operator+=(s);
      return *this;
    }
  RWCString& operator+=(const std_string &s)
    {
      std_string::operator+=(s);
      return *this;
    }
  RWCString& operator+=(const char &c)
    {
      std_string::operator+=(c);
      return *this;
    }
  RWBoolean operator==(const RWCString &s) const
    {
      const string *a=this;
      const string *b=&s;
      return *a==*b;
    }
  RWBoolean operator==(const char* s) const
    {
      const string *a=this;
      const string b=s;
      return *a==b;
    }
  RWBoolean operator!=(const RWCString &s) const
    {
      const string *a=this;
      const string *b=&s;
      return *a!=*b;
    }
  RWBoolean operator!=(const char* s) const
    {
      const string *a=this;
      const string b=s;
      return *a!=b;
    }

  //inline RWCString operator () (unsigned int, unsigned int);
  inline operator const char*() const {return c_str();};
  RWCString& prepend (const char*str)
    {
      insert(0,str);
      return *this;
    }

  RWCString& append(const RWCString&s)
    {
      std_string::operator+=(s);
      return *this;
    }
  inline std_istream& readLine ( std_istream&, RWBoolean skipWhite=true);
  inline RWCString& replace (unsigned int, unsigned int, 
                             const char*, unsigned int );
  RWCString& replace(size_t pos,size_t n,const char*s)
    {
      std_string::replace(pos,n,s);
      return *this;
    }
  RWCString& remove(size_t n)
    {
      if(n<size())
	erase(n,size()-n);
      return *this;
    }
  RWCString& remove(size_t pos,size_t N)
    {
      erase(pos,N+pos);
      return *this;
    }
  int first(char c) const{return find(c);}
  int last(char c) const{return rfind(c);}

  RWBoolean contains(std_string s) const
    {
      int i=find(s);
      if(i!=std_string::npos) return true;
      return false;
    }

  RWBoolean contains(char c) const
    {
      int i=find(c);
      if(i!=std_string::npos) return true;
      return false;
    }
      
  //
  // stripType = 0 beginning
  // stripType = 1 end
  // stripType = 2 both
  //
  inline RWCString strip ( int stripType, char );
  inline void toLower ( void );
  inline void toUpper ( void );
  inline RWBoolean isNull() const {return empty ();}
  inline size_t index (const char*,int pos=0) const; 
  inline size_t index (char,int pos=0) const; 
  inline size_t index( const RWCString&, size_t, size_t, caseCompare ) const;
  RWCSubString operator()(size_t start, size_t extent)
    {
      return RWCSubString(*this,start,extent);
    }

  const char* data() const{return c_str();}

  unsigned int hash( caseCompare cmp = exact ) const
    {
      const char*s=c_str();
      unsigned long h = 0;
      for ( ; *s; ++s)
	h = 5*h + *s;
      
      return size_t(h);
    }

  unsigned int stlhash() const
    {
      const char*s=c_str();
      unsigned long h = 0;
      for ( ; *s; ++s)
	h = 5*h + *s;
      
      return size_t(h);
    }

  // useful for supplying hash functions to template hash collection ctors:
  static unsigned hash(const RWCString&s)
    {
      return s.hash();
    }

};


RWCSubString::RWCSubString(const RWCSubString&s)
{
  mystring=s.mystring;
  mystart=s.mystart;
  extent=s.extent;
}

RWCSubString& RWCSubString::operator=(const RWCString&s)
{
  RWCString str(s);
  return operator=(str);
}
RWCSubString& RWCSubString::operator=(const RWCSubString&s)
{
  mystring->replace(mystart,extent,s.mystring->data(),s.length());
  extent=s.length();
  return *this;
}

RWCSubString& RWCSubString::operator=(const char*s)         
{
  mystring->replace(mystart,extent,s,strlen(s));
  extent=strlen(s);
  return *this;
}

char& RWCSubString::operator()(size_t i)
{return mystring->operator[](mystart+i);}

char  RWCSubString::operator()(size_t i) const
{return mystring->operator[](mystart+i);}

char&  RWCSubString::operator[](size_t i)
{return mystring->operator[](mystart+i);}

char  RWCSubString::operator[](size_t i) const
{return mystring->operator[](mystart+i);}

RWBoolean RWCSubString::operator==(RWCString s) const
{
  return mystring->substr(mystart,extent)==s ? true:false;
}

RWBoolean RWCSubString::operator==(const char* s) const
{
  return mystring->substr(mystart,extent)==string(s) ? true:false;
}

RWBoolean RWCSubString::operator!=(RWCString s) const
{
  return mystring->substr(mystart,extent)!=string(s) ? true:false;
}

RWBoolean RWCSubString::operator!=(const char* s) const
{
  return mystring->substr(mystart,extent)!=string(s) ? true:false;
}

//#include <rw/cstring.icc>

//GB:OSF1-cxx : #include <cstring>

RWCString::RWCString ( const char * astring ) : std_string ( astring ) {}

void RWCString::toLower ( void ) {
  for (int i=0; i<=size();i++) {
    //GB:HP-UX-aCC,Linux-KCC 
    std_string::operator[](i) = tolower(std_string::operator[](i));
    //at(i) = tolower(at(i)); 
  } 
}

void RWCString::toUpper ( void ) {
  for (int i=0; i<=size();i++) {
    //GB:HP-UX-aCC,Linux-KCC 
    std_string::operator[](i) = tolower(std_string::operator[](i));
    //at(i) = toupper(at(i)); 
  }
}

std_istream& RWCString::readLine ( std_istream& s, RWBoolean skipWhite/*GB:OSF1-cxx =true*/ ) {
  char tmp[256];
  if ( skipWhite ) {
    s >> std_ws;
    s.getline(tmp,256);
    *this=tmp;
  }
  else {
    s.getline(tmp,256);    
    *this=tmp;
  } 
  return s;
}
RWCString& RWCString::replace (unsigned int start,  unsigned int nbytes, 
                               const char* buff, unsigned int n2 ) {
  std_string::replace ( start, nbytes, buff, n2 ); 
  return *this;                                                              
}                                                                          

RWCString RWCString::strip ( int stripType=trailing, char c=' ' ) {
  RWCString retVal = *this;
  if(length()==0) return retVal;
  int i;
  switch ( stripType ) {
  case leading: 
    {
      for(i=0;i<length();i++)
	if (std_string::operator[](i) != c) break;
      retVal = substr(i,length()-i);
    }
    break;
  case trailing:
    { 
      for(i=length()-1;i>=0;i--)
	if (std_string::operator[](i) != c) break; 
      retVal = substr(0,i+1);
    }
    break;
  case both:
    { 
      for(i=0;i<length();i++)
	if (std_string::operator[](i) != c) break;
      RWCString tmp(substr(i,length()-i));
      for(i=tmp.length()-1;i>=0;i--)
	if (tmp.std_string::operator[](i) != c) break; 
      retVal = tmp.substr(0,i+1);
    }
  }
  return retVal;
}

// "cs" optional parameter is NOT implemented ! 
size_t RWCString::index( const RWCString& s, size_t ln, 
                         size_t st, RWCString::caseCompare cs ) const {
  return std_string::find( s.c_str(), st, ln );
}

// "cmp" optional parameter is NOT implemented ! 
// N.B.: The hash value returned is generally DIFFERENT from the
//       one returned by the original RW function.
//       Users should not rely on the specific return value.

char RWCString::operator () (size_t i) const
{
  return operator[](i);
}

size_t RWCString::index (const char* str,int pos) const
{
  return std_string::find(str,pos);
}

size_t RWCString::index (char c,int pos) const
{
  return std_string::find(c,pos);
}


#endif









