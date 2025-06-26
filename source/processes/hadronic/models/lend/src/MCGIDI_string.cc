/*
 # <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
 # <<END-copyright>>
 * 
 * IMPLEMENTATION NOTES
 *
 * My goal was to create a small, and conceptually simple
 * string class.  I allocate unique memory for every object
 * created, instead of doing anything fancy.
 *
 * I've taken care about the order in which I free memory,
 * in case we are copying overlapping regions.
 *
 */

#include <stdio.h>
#include <stdlib.h> // malloc, realloc
#include <string.h>
#include <assert.h>

#include "MCGIDI_string.hpp"

namespace MCGIDI {

#define MCGIDI_MIN(x,y) ( ( (x) < (y) ) ? (x) : (y) )
#define MCGIDI_SWAP(a,b,type) {type ttttttttt=a;a=b;b=ttttttttt;}
LUPI_HOST_DEVICE int MCGIDI_strcmp (const char *p1, const char *p2);
LUPI_HOST_DEVICE size_t MCGIDI_strlen (const char *str);
LUPI_HOST_DEVICE void MCGIDI_memmove(char *dest, const char *src, size_t n);
LUPI_HOST_DEVICE int MCGIDI_strncmp( const char * s1, const char * s2, size_t n );

  const String::size_type String::npos = static_cast<size_t>(-1);

  /*
   * Like the 'new' operator, we want to guarantee that we NEVER
   * return nullptr. Loop until there is free memory.
   *
   */
  LUPI_HOST_DEVICE static char* malloc_never_null(const size_t b)
  {
      char *p;

      do {
          p = static_cast<char*>(malloc(b));
      } while ( p == nullptr );

      return p;
  }

  /**
   * Allocates memory for the copy the string at the same time sets "this->allocated".
   * @param s
   * @return
   */
  LUPI_HOST_DEVICE char* String::strdup_never_null(const char* s)
  {
      const size_t len = MCGIDI_strlen(s)+1;
      char *p2 = malloc_never_null(len);
      memcpy(p2, s, len);
      allocated_=len;
      size_=len-1;
      return p2;
  }

  LUPI_HOST_DEVICE String::String() : p( nullptr ), allocated_(0), size_(0) { }

  LUPI_HOST_DEVICE String::~String()
  {
      free(p);
  }

  LUPI_HOST_DEVICE String::String(const String& s)
    : p(nullptr)
  {
    p = malloc_never_null( s.size_ + 1 );  // copy only used part
    allocated_ = s.size_ + 1;
    size_      = s.size_;
    memcpy(p, s.p, size_ + 1);
  }

  LUPI_HOST_DEVICE String::String(const char* s)
  : p(strdup_never_null(s))
  {
  }

  LUPI_HOST_DEVICE String& String::operator=(const char* s)
  {
      if ( p != s ) {
        // s could point into our own string, so we have to allocate a new string
        const size_t len = MCGIDI_strlen(s);
        char* copy = (char*) malloc( len + 1);
        MCGIDI_memmove(copy, s, len+1); // trailing 0
        free( p );
        p = copy;
        size_ = len;
        allocated_ = len+1;
      }

      return *this;
  }

  LUPI_HOST_DEVICE String& String::operator=(const String& s)
  {
      return operator=(s.p);
  }

  LUPI_HOST_DEVICE String& String::operator+=(const String& s)
  {
    if (s.size_ > 0){
      this->reserve(size_ + s.size_);
      MCGIDI_memmove(p+size_, s.p, s.size_+1); // trailing 0
      size_ += s.size_;
    }
    return *this;
  }

  // since p and s may overlap, we have to copy our own string first
  LUPI_HOST_DEVICE String& String::operator+=(const char* s)
  {
      const size_type lens = MCGIDI_strlen(s);
      if (lens > 0){
        if (size_ + lens + 1 <= allocated_) {
          MCGIDI_memmove(p+size_, s, lens+1); // trailing 0
          size_ += lens;
        }
        else {
          String s2( *this );  // copy own data
          s2.reserve(size_ + lens);
          MCGIDI_memmove(s2.p+size_, s, lens+1); // trailing 0
          s2.size_ = size_ + lens;
          this->swap( s2 );
        }
      }
      return *this;
  }

  LUPI_HOST_DEVICE String& String::operator+=(const char c)
  {
      push_back(c);
      return *this;
  }


  LUPI_HOST_DEVICE void String::push_back(const char c) {

    if (size_ == allocated_ - 1) {
      size_t more =  (allocated_* 3) / 2; // factor 1.5
      if ( more < 4 ) more = 4;
      reserve( size_ + more );
    }

    p[size_] = c;
    size_++;
    p[size_] = 0;
  }

  LUPI_HOST_DEVICE bool String::operator==(const char* s) const
  {
      return !MCGIDI_strcmp(p, s);
  }

  LUPI_HOST_DEVICE bool String::operator==(const String& s) const
  {
      return !MCGIDI_strcmp(p, s.p);
  }

  LUPI_HOST_DEVICE void String::clearMemory()
  {
    String s;
    this->swap( s );
  }

  LUPI_HOST_DEVICE void String::clear()
  {
      size_ = 0;
      p[0]  = 0;
  }

  LUPI_HOST_DEVICE String operator+(const String& lhs, const String& rhs)
  {
      return String(lhs) += rhs;
  }

  LUPI_HOST_DEVICE String String::substr(const size_type pos, size_type in_length) const
  {
      String s;
      const size_type len = size_;

      if ( pos > len )
        LUPI_THROW("MCGIDI::String::substr: pos index out of range");

      size_type remain = len - pos;

      if ( in_length > remain )
          in_length = remain;

      s.reserve( in_length );

      memcpy(s.p, p + pos, in_length);
      s.p[in_length] = '\0';
      s.size_ = in_length;

      return s;
  }


  // checked access, accessing the NUL at end is allowed
  LUPI_HOST_DEVICE char& String::at(const size_type i)
  {
      if ( i > MCGIDI_strlen(p) )
        LUPI_THROW( "MCGIDI::String::at(): index out_of_range");

      return p[i];
  }
  LUPI_HOST_DEVICE char String::at(const size_type i) const
  {
      if ( i > MCGIDI_strlen(p) )
        LUPI_THROW("MCGIDI::String::at(): index out_of_range");

      return p[i];
  }

  LUPI_HOST_DEVICE String& String::erase(size_type pos, size_type len)
  {
    if (len > 0) {

      if ( pos >= size_ ) // user must not remove trailing 0
        LUPI_THROW("MCGIDI::String::erase: pos index out_of_range");

      long s2 = size_;
      long remain = s2 - (long) pos - len;

      if (remain > 0) {
        // erase by overwriting
        MCGIDI_memmove(p + pos, p + pos + len, remain);
      }

      if ( remain < 0 ) remain = 0;

      // remove unused space
      this->resize( pos+remain );

    }
    return *this;
  }

  LUPI_HOST_DEVICE String& String::append( const char* str, size_type n) {
    if (str && n > 0) {
      size_t lens = MCGIDI_strlen(str);
      if (n > lens)
        n = lens;
      size_t newlen = size_ + n;
      this->reserve( newlen );
      MCGIDI_memmove(p+size_, str, n); // p and s.p MAY overlap
      p[newlen] = 0; // add NUL termination
      size_ = newlen;
    }
    return *this;
  }

  LUPI_HOST_DEVICE int String::compare( size_type pos, size_type len, const String& str ) const {
    if (pos > size_)
      LUPI_THROW("MCGIDI::String::compare: pos index out of range");

    if ( len > size_ - pos)
      len = size_ - pos; // limit len to available length

    const size_type osize = str.size();
    const size_type len2   = MCGIDI_MIN(len, osize);
    int r = MCGIDI_strncmp( p + pos, str.p, len2);
    if (r==0) // equal so far, now compare sizes
    r = len < osize ? -1 : ( len == osize ? 0 : +1 );
    return r;
  }

  LUPI_HOST_DEVICE int String::compare( size_type pos, size_type len, const char* str ) const {
    if (pos > size_)
      LUPI_THROW("MCGIDI::String::compare: pos index out of range");

    if ( len > size_ - pos)
      len = size_ - pos; // limit len to available length

    const size_type osize = MCGIDI_strlen(str);
    const size_type len2   = MCGIDI_MIN(len, osize);
    int r = MCGIDI_strncmp( p + pos, str, len2);
    if (r==0) // equal so far, now compare sizes
        r = len < osize ? -1 : ( len == osize ? 0 : +1 );
    return r;
  }


  LUPI_HOST_DEVICE void String::my_realloc( size_type n, char ** address) {
    if (address != nullptr && *address != nullptr) {
        p = *address;
        long delta = sizeof(char) * n;
        long sub = delta % 8;
        if (sub != 0) delta += (8-sub);
        *address += delta;
        return;
    }
    if (n > 0 ) {
      char* pnew = static_cast<char*>(malloc(n)); // could return nullptr
      if (pnew) {
        free(p);
        p = pnew;
      }
      else
        LUPI_THROW("MCGIDI::String::my_realloc out of memory");
    }
   }


  LUPI_HOST_DEVICE void String::reserve( const size_type n, char ** address) {
    if (n >= allocated_ ) {
      this->my_realloc(n + 1, address);
      allocated_ = n + 1;
    }
  }

  LUPI_HOST_DEVICE void String::resize( const size_type n, char ** address) {
    this->resize( n, 0, address );
  }

  LUPI_HOST_DEVICE void String::resize( const size_type n, const char c, char ** address) {
    if (n < allocated_ ) {
      p[n] = 0;
      size_ = n;
    }
    else if (n >= allocated_ ) {
      this->reserve( n, address );
      for (size_type i=size_; i < n; ++i )
        p[i] = c;
      p[n] = 0;
      size_ = n;
    }
  }

  LUPI_HOST_DEVICE void String::swap( String& s ) {
    MCGIDI_SWAP( allocated_, s.allocated_, size_t );
    MCGIDI_SWAP( size_,      s.size_,      size_t );
    MCGIDI_SWAP( p,          s.p,          char * );
  }


  // Comparison
  LUPI_HOST_DEVICE bool operator<( const String& s1, const String& s2 ) {
    return MCGIDI_strcmp( s1.c_str(), s2.c_str() ) < 0;
  }

  /* Compare S1 and S2, returning less than, equal to or
     greater than zero if S1 is lexicographically less than,
     equal to or greater than S2.  */
  LUPI_HOST_DEVICE int MCGIDI_strcmp (const char *p1, const char *p2)
  {
    const unsigned char *s1 = (const unsigned char *) p1;
    const unsigned char *s2 = (const unsigned char *) p2;
    unsigned char c1, c2;

    do
      {
        c1 = (unsigned char) *s1++;
        c2 = (unsigned char) *s2++;
        if (c1 == '\0')
          return c1 - c2;
      }
    while (c1 == c2);

    return c1 - c2;
  }

  LUPI_HOST_DEVICE size_t MCGIDI_strlen (const char *str) {
    size_t len = 0;
    while (*str != '\0') {
        str++;
        len++;
    }
    return len;
  }

  // A function to copy block of 'n' bytes from source
  // address 'src' to destination address 'dest'.
  LUPI_HOST_DEVICE void MCGIDI_memmove(char *dest, const char *src, size_t n)
  {
     // Typecast src and dest addresses to (char *)
     char *csrc = (char *)src;
     char *cdest = (char *)dest;
 
     // Create a temporary array to hold data of src
     char* temp = (char*) malloc( n);
 
     // Copy data from csrc[] to temp[]
     for (size_t i=0; i<n; i++)
         temp[i] = csrc[i];
 
     // Copy data from temp[] to cdest[]
     for (size_t i=0; i<n; i++)
         cdest[i] = temp[i];
 
     free(temp);
  }

  LUPI_HOST_DEVICE int MCGIDI_strncmp( const char * s1, const char * s2, size_t n )
  {
      while ( n && *s1 && ( *s1 == *s2 ) )
      {
          ++s1;
          ++s2;
          --n;
      }
      if ( n == 0 )
      {
          return 0;
      }
      else
      {
          return ( *(unsigned char *)s1 - *(unsigned char *)s2 );
      }
  }

}
