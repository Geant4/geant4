/*
 # <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
 # <<END-copyright>>
 */

#ifndef MCGIDI_STRING_HPP
#define MCGIDI_STRING_HPP

/* Modified from Karsten Burger's version 2017
 * Made changes to make it more compatible with GPUs.
 * Allow for data to be initialized to nullptr.
 *
 * Modified from public domain software:
 * Karsten Burger 2014
 *
 * Sourceforge project "Simple C++ String Class"
 * http://sourceforge.net/projects/simplecstringclass/

 * This a simple C++ string class based on class my_string by
 * Christian Stigen Larsen, 2007, http://csl.name/programming/my_string/
 *
 * It only uses the C-string functions and is thus independent of the
 * standard C++ library.
 *
 * It is public domain, in the hope, that you find it useful.
 * Please note that there is no guarantee of any kind: it is supplied
 * without any warranty; without even the implied warranty of
 * merchantability or fitness for a particular purpose.
 *
 * You can probably replace std::string with this one in many
 * cases, but a lot of stuff is missing, and I would recommend
 * you stick to std::string anyway.
 *
 * I want to point out that there is nothing fancy about this class.
 * It keeps every string in its own buffer, and copies as often as
 * needed.
 * The data always contains a trailing NUL char.
 *
 * However, I believe that is a good approach.  For instance, it
 * uses malloc rather than new, which makes it possible to use
 * realloc.  On many systems, realloc will try to use up "invisible"
 * space that was used by malloc to pad a string for memory alignment.
 * That makes it potentially fast for small concatenations.
 *
 * I don't propose to use this class for anything practical, since
 * we already have std::string, but it may be an interesting read
 * for C++ novices at the very least. Also, additional functions can
 * easily be expanded.
 *
 * Also I met a case, where I had to avoid std::string because of
 * link problems with an application using mixed libraries, especially
 * one compiled with an old Intel compiler icc 7.
 *
 * Bugs/suggestions to info [at) dr-burger ]dot[ com
 * or via the Sourceforge project page.
 */

#include <sys/types.h> // size_t
#include <stdexcept>
#include <LUPI_declareMacro.hpp>

  /** @brief Simple C++ string class, useful as replacement for
   std::string if this cannot be used, or just for fun.

   */
namespace MCGIDI {

class String
{

    char*     p;           ///< The data
    size_t    allocated_;  ///< The allocated memory size (including trailing NUL)
    size_t    size_;       ///< The currently used memory size (excluding trailing NUL)

  public:
    typedef size_t size_type;
    static const   size_type npos;

    LUPI_HOST_DEVICE String();
    LUPI_HOST_DEVICE ~String();
    LUPI_HOST_DEVICE String(const String&);
    LUPI_HOST_DEVICE String(const char*);

    LUPI_HOST_DEVICE String&  operator=(const char*);
    LUPI_HOST_DEVICE String&  operator=(const String&);

    LUPI_HOST_DEVICE String&  operator+=(const String&);
    LUPI_HOST_DEVICE String&  operator+=(const char*);
    LUPI_HOST_DEVICE String&  operator+=(char);
    LUPI_HOST_DEVICE void     push_back(char);

    friend String
    LUPI_HOST_DEVICE operator+(const String& lhs, const String& rhs);

    LUPI_HOST_DEVICE bool     operator==(const char*) const;
    LUPI_HOST_DEVICE bool     operator==(const String&) const;

    LUPI_HOST_DEVICE void     clear();       // set string to empty string (memory remains reserved)
    LUPI_HOST_DEVICE void     clearMemory(); // set string to empty string (memory is free'd)

    LUPI_HOST_DEVICE size_type size()   const   { return size_; }   ///< size without terminating NUL
    LUPI_HOST_DEVICE size_type length() const   { return size_; }   ///< as size()

    // size if fully used
    LUPI_HOST_DEVICE size_type capacity() const { return allocated_-1; }

    // 8 byte alligned size
    LUPI_HOST_DEVICE long internalSize() const {
        long delta = allocated_;
        long sub = delta % 8;
        if (sub != 0) delta += (8-sub);
        return delta * sizeof(char);
    }

    LUPI_HOST_DEVICE bool      empty() const    { return size_ == 0; }

    LUPI_HOST_DEVICE const char*  c_str() const { return p; } ///< raw data

    /** Reserve internal string memory so that n characters can be put into the
        string (plus 1 for the NUL char). If there is already enough memory,
        nothing happens, if not, the memory will be realloated to exactly this
        amount.
        */
    LUPI_HOST_DEVICE void reserve( size_type n, char ** address = nullptr);

    /** Resize string. If n is less than the current size, the string will be truncated.
        If n is larger, then the memory will be reallocated to exactly this amount, and
        the additional characters will be NUL characters.
        */
    LUPI_HOST_DEVICE void resize( size_type n, char ** address = nullptr);

    /** Resize string. If n is less than the current size, the string will be truncated.
        If n is larger, then the memory will be reallocated to exactly this amount, and
        the additional characters will be c characters.
        */
    LUPI_HOST_DEVICE void resize( size_type n, char c, char ** address = nullptr);

    /// swap contents
    LUPI_HOST_DEVICE void swap( String& );

    LUPI_HOST_DEVICE String   substr(const size_type pos, size_type length) const;

    // unchecked access:
    LUPI_HOST_DEVICE char&    operator[](const size_type i)       { return p[i]; }
    LUPI_HOST_DEVICE char     operator[](const size_type i) const { return p[i]; }
    // checked access:
    LUPI_HOST_DEVICE char&    at(const size_type i);
    LUPI_HOST_DEVICE char     at(const size_type i) const;

    /// erase len characters at position pos
    LUPI_HOST_DEVICE String& erase(size_type pos, size_type len);
    /// Append n characters of a string
    LUPI_HOST_DEVICE String& append(const char* str, size_type n);

    LUPI_HOST_DEVICE int     compare( size_type pos, size_type len, const String& str ) const;
    LUPI_HOST_DEVICE int     compare( size_type pos, size_type len, const char*   str ) const;

  private:
    // reallocate the internal memory
    LUPI_HOST_DEVICE void  my_realloc( size_type n, char ** address = nullptr);
    LUPI_HOST_DEVICE char* strdup_never_null(const char* other);

};
// class

LUPI_HOST_DEVICE bool operator<(const String&, const String&);

}

#endif
