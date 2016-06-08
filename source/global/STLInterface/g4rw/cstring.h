// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: cstring.h,v 1.6 1999/11/29 10:17:39 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4String
//
//  Class description:
//
//  Definition of a Geant4 string.
//  It's currently implemented as an STL wrapper class for Strings.
//  It implements Rogue Wave RWCString signature but intrinsically
//  using STL string.

//---------------------------------------------------------------

#ifndef __cstring
#define __cstring

#include <string>
#include "G4Types.hh"
#include "g4rw/defs.h"

#ifdef __KCC
  #include <istream>
#endif

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

class G4String;

class G4SubString
{
public:

  inline G4SubString(const G4SubString&);

  inline G4SubString& operator=(const char*);         

  inline G4SubString& operator=(const G4String&);
  inline G4SubString& operator=(const G4SubString&);
 
  inline char& operator()(size_t);
  inline char  operator()(size_t) const;
  inline char& operator[](size_t);
  inline char  operator[](size_t) const;

  inline int operator!() const;

  inline G4bool operator==(G4String) const;
  inline G4bool operator==(const char*) const;
  inline G4bool operator!=(G4String) const;
  inline G4bool operator!=(const char*) const;

  inline size_t length() const;
  inline size_t start() const;

  // For detecting null substrings
  //
  inline G4bool isNull() const;

private:

  inline G4SubString(G4String&, size_t, size_t);

  G4String*    mystring;     
  size_t        mystart;
  size_t        extent;

  friend class G4String;

};
 

class G4String : public G4std::string
{

  typedef G4std::string std_string;

public: 

  enum caseCompare { exact, ignoreCase };
  enum stripType { leading, trailing, both };

  inline G4String ();
  inline G4String ( char );
  inline G4String ( const char * );
  inline G4String ( const G4String& );
  inline G4String ( const G4SubString& );
  inline G4String ( const std_string& );
  virtual ~G4String () {}

  inline G4String& operator=(const G4String&);
  inline G4String& operator=(const std_string&);
  inline G4String& operator=(const char*);

  inline char operator () (size_t) const; 
  inline char& operator () (size_t);

  inline G4String& operator+=(const G4SubString&);
  inline G4String& operator+=(const char*);
  inline G4String& operator+=(const std_string&);
  inline G4String& operator+=(const char&);
  inline G4bool operator==(const G4String&) const;
  inline G4bool operator==(const char*) const;
  inline G4bool operator!=(const G4String&) const;
  inline G4bool operator!=(const char*) const;

  //inline G4String operator () (unsigned int, unsigned int);
  inline operator const char*() const;
  inline G4SubString operator()(size_t, size_t);

  inline int compareTo(const char*, caseCompare mode=exact);
  inline int compareTo(const G4String&, caseCompare mode=exact);

  inline G4String& prepend (const char*);
  inline G4String& append (const G4String&);

  inline G4std::istream& readLine (G4std::istream&, G4bool skipWhite=true);
  
  inline G4String& replace (unsigned int, unsigned int, 
                             const char*, unsigned int );
  inline G4String& replace(size_t, size_t, const char*);

  inline G4String& remove(size_t);
  inline G4String& remove(size_t, size_t);

  inline int first(char) const;
  inline int last(char) const;

  inline G4bool contains(std_string) const;
  inline G4bool contains(char) const;

  // stripType = 0 beginning
  // stripType = 1 end
  // stripType = 2 both
  //
  inline G4String strip (int stripType=trailing, char c=' ');

  inline void toLower ( void );
  inline void toUpper ( void );

  inline G4bool isNull() const;

  inline size_t index (const char*, int pos=0) const; 
  inline size_t index (char, int pos=0) const; 
  inline size_t index (const G4String&, size_t, size_t, caseCompare) const;

  inline const char* data() const;

  inline unsigned int hash( caseCompare cmp = exact ) const;
  inline unsigned int stlhash() const;

  // useful for supplying hash functions to template hash collection ctors
  //
  static inline unsigned hash(const G4String&);

};

#include "g4rw/cstring.icc"

#endif
