

//



//
// $Id: scl_string.cc,v 1.2 1999-05-21 20:21:13 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST Utils Class Library
* clutils/scl_string.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <scl_string.h>
#include <stdio.h>
/*******************************************************************/

int SCLstring::newBufSize (int len) const
{
  int rval =  (((len + 1)/32) +1) * 32;

  // if there\'s a max and rval > max return max +1
  return (MaxLength () && (rval > MaxLength ())) ?  MaxLength ()+1 : rval;
}

//constructor
SCLstring::SCLstring (const char * str, int lim)
:  _max_len (lim)
{

  if (!str) { _strBufSize = 0; _strBuf = 0;  return;  }

  // check the length
  int len;
  if (lim && (strlen (str) > lim)) len = lim;
  else len = strlen (str);

  _strBufSize = newBufSize (len);
  _strBuf = new char[_strBufSize];
  strncpy (_strBuf, str, len);
  _strBuf [len] = '\0';

}	

SCLstring::SCLstring (const SCLstring& s)  
:  _max_len (s.MaxLength ())
{
  if (s.is_null ()) { _strBufSize = 0; _strBuf = 0;  return;  }
//  int len = newBufSize (s.Length ());
//  _strBuf = new char [len];
//  strncpy (_strBuf, s, len);
//  _strBuf [len] = '\0';
  _strBufSize = newBufSize (s.Length ());
  _strBuf = new char [_strBufSize];
  strncpy (_strBuf, s, _strBufSize);
  _strBuf [_strBufSize] = '\0';
}  

//destructor
SCLstring::~SCLstring()
{
    if(_strBufSize)
    {
	delete [] _strBuf;
//	_strbuf = 0;
    }
}

//operators

SCLstring& SCLstring::operator= (const char* str) 
{
  if (!str) { 
    _strBufSize = 0; delete [] _strBuf; _strBuf =0;  
  }
  else {
    int slen = strlen (str);
    if(!_strBuf || (StrBufSize () < newBufSize(slen+1) ) )  {
      // make more room
      // first clean up the space
      if (_strBuf)  delete [] _strBuf;
      _strBufSize = newBufSize (slen+1);
      _strBuf = new char[_strBufSize];
      }
    strncpy (_strBuf, str, _strBufSize -1);
    _strBuf [_strBufSize -1] = '\0';
  }
    return *this;
}

// Josh L, 4/11/95
SCLstring& SCLstring::operator= (const SCLstring& s)
{
    // changed from sending s.chars() so undefined will be set correctly - DAS
  operator=(s.rep());
  _strBufSize = s.StrBufSize();
  _max_len = s.MaxLength();
  return *this;
}

SCLstring::operator const char * ()  const
{
    if (_strBuf) return _strBuf;
    else return "";
}

SCLstring& 
SCLstring::Append (const char * s)  
{
  if (!s) return *this;
  int olen = Length (), slen = strlen (s);
  int len = olen + slen +1;
  if (_strBuf && (len < _strBufSize))  strncat (_strBuf, s, slen+1);
  else {  // make more space
    int sz = newBufSize (len);
    char * tmp = new char [sz];  tmp [0] ='\0';
    _strBufSize = sz;
    if (_strBuf) strcpy (tmp, _strBuf );
    len = ((slen > sz-olen) ? sz-olen : slen);
    strncat (tmp, s, len);
    *(tmp+olen+len) = '\0';  
    delete [] _strBuf;
    _strBuf = tmp;
    }
  return *this;
}

SCLstring& 
SCLstring::Append (const char c)  
{
  char tmp [2];
  tmp [0] = c;
  tmp [1] = '\0';
  return Append (tmp);
}
   
SCLstring& 
SCLstring::Append (const long int i)
{
  char tmp [BUFSIZ];
  sprintf (tmp, "%ld", i);
  return Append (tmp);
}

SCLstring& 
SCLstring::Append (const int i)
{
  char tmp [BUFSIZ];
  sprintf (tmp, "%d", i);
  return Append (tmp);
}

SCLstring& 
SCLstring::Append (const double i, const int prec)
{
  char tmp [BUFSIZ];
  sprintf (tmp, "%.*g", prec, i);
  return Append (tmp);
} 

SCLstring& 
SCLstring::Prepend (const char * s)  
{
  if (!_strBuf)  {  // no string 
    int sz = newBufSize (strlen (s));
    _strBuf = new char [sz];  
    sz = sz-1;
    strncpy (_strBuf, s, sz);
    _strBuf [sz] = '\0';
    return *this;
  }
  //  otherwise make some space
  int slen = strlen (s);
  int sz = newBufSize (slen + Length () );
  char * nstr = new char [sz];  
  strncpy (nstr, s, sz-1);
  if (slen < sz)  // copy the rest
    strncat (nstr, _strBuf, sz-slen);
  nstr [sz -1] = '\0';  // make sure it\'s null terminated
  delete [] _strBuf;
  _strBuf = nstr;
  return *this;
}
   

int 
SCLstring::operator== (const char* s) const
{
  return !strcmp (_strBuf, s);  
}

// is_null returns true if the string doesn't exist or if it is empty.
// If is_null() returns true you may also want to call is_undefined() to 
// see if it is undefined.
int 
SCLstring::is_null() const
{
     if (_strBuf) return _strBuf[0] == '\0';
     else return 1;
}

// If the buffer doesn't exist SCLstring returns true for is_null() and 
// and is_undefined().
int
SCLstring::is_undefined() const
{
  return !_strBuf;
}

// see also set_null()
void
SCLstring::set_undefined()
{
    delete [] _strBuf;
    _strBuf = 0;
    _strBufSize = 0;
}

// see also set_undefined()
int 
SCLstring::set_null() 
{
  for (int i =0; i < _strBufSize; ++i)
    *(_strBuf +i) = '\0';
  return 1;
/*
  delete [] _strBuf;
  _strBuf =0;
  _strBufSize = 0;
  return 1;
  */
}

int
SCLstring::StrBufSize() const
{ return _strBufSize; }

int
SCLstring::Length() const
{ 
  if (_strBuf) return strlen (_strBuf); 
  else return 0;
}
