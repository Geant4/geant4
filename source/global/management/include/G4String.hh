//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4String.hh,v 1.3 2002-03-25 15:32:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//  Derived from the Rogue Wave implementation of RWCString;
//  it uses intrinsically STL string.

//---------------------------------------------------------------

#ifndef __G4String
#define __G4String

#include <stdio.h>
#include <string>
#include "G4Types.hh"
#include "g4std/iostream"

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#ifdef G4USE_STD_NAMESPACE
  typedef G4std::string::size_type str_size;
#else
  typedef size_t str_size;
#endif

class G4String;

class G4SubString
{
public:

  inline G4SubString(const G4SubString&);

  inline G4SubString& operator=(const char*);         

  inline G4SubString& operator=(const G4String&);
  inline G4SubString& operator=(const G4SubString&);
 
  inline char& operator()(str_size);
  inline char  operator()(str_size) const;
  inline char& operator[](str_size);
  inline char  operator[](str_size) const;

  inline G4int operator!() const;

  inline G4bool operator==(G4String) const;
  inline G4bool operator==(const char*) const;
  inline G4bool operator!=(G4String) const;
  inline G4bool operator!=(const char*) const;

  inline str_size length() const;
  inline str_size start() const;

  // For detecting null substrings
  //
  inline G4bool isNull() const;

private:

  inline G4SubString(G4String&, str_size, str_size);

  G4String*    mystring;     
  str_size     mystart;
  str_size     extent;

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
  inline G4String ( const G4std::string & );
  virtual ~G4String () {}

  inline G4String& operator=(const G4String&);
  inline G4String& operator=(const G4std::string &);
  inline G4String& operator=(const char*);

  inline char operator () (str_size) const; 
  inline char& operator () (str_size);

  inline G4String& operator+=(const G4SubString&);
  inline G4String& operator+=(const char*);
  inline G4String& operator+=(const G4std::string &);
  inline G4String& operator+=(const char&);
  inline G4bool operator==(const G4String&) const;
  inline G4bool operator==(const char*) const;
  inline G4bool operator!=(const G4String&) const;
  inline G4bool operator!=(const char*) const;

  //inline G4String operator () (unsigned int, unsigned int);
  inline operator const char*() const;
  inline G4SubString operator()(str_size, str_size);

  inline G4int compareTo(const char*, caseCompare mode=exact);
  inline G4int compareTo(const G4String&, caseCompare mode=exact);

  inline G4String& prepend (const char*);
  inline G4String& append (const G4String&);

  inline G4std::istream& readLine (G4std::istream&, G4bool skipWhite=true);
  
  inline G4String& replace (unsigned int, unsigned int, 
                             const char*, unsigned int );
  inline G4String& replace(str_size, str_size, const char*);

  inline G4String& remove(str_size);
  inline G4String& remove(str_size, str_size);

  inline G4int first(char) const;
  inline G4int last(char) const;

  inline G4bool contains(G4std::string) const;
  inline G4bool contains(char) const;

  // stripType = 0 beginning
  // stripType = 1 end
  // stripType = 2 both
  //
  inline G4String strip (G4int stripType=trailing, char c=' ');

  inline void toLower ( void );
  inline void toUpper ( void );

  inline G4bool isNull() const;

  inline str_size index (const char*, G4int pos=0) const; 
  inline str_size index (char, G4int pos=0) const; 
  inline str_size index (const G4String&, str_size, str_size, caseCompare) const;

  inline const char* data() const;

  inline G4int strcasecompare(const char*, const char*) const;

  inline unsigned int hash( caseCompare cmp = exact ) const;
  inline unsigned int stlhash() const;

  // useful for supplying hash functions to template hash collection ctors
  //
  static inline unsigned hash(const G4String&);

};

#include "G4String.icc"

#endif
