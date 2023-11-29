//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4String
//
// Class description:
//
// Common string type for Geant4
//
// Provides a `std::string` compliant implementation of basic
// strings. It currently inherits from `std::string`, but this should
// not be assumed other than that it will implement the same interface
// as `std::string` as of the minimum C++ standard required by Geant4
// (currently C++17).
//
// It can be passed to any function accepting `std::string` or `std::string_view` 
// arguments. Whilst it currently implements a conversion operator to `const char*`,
// this should be considered deprecated for use outside of Geant4. Passing to
// non-Geant4 functions taking `const char*` arguments should use the
// `std::string::c_str` member function to explicitly convert.
// 
// See `std::string` for primary interfaces, `G4StrUtil` for additional query/manipulation functions
//
// Author: G.Cosmo, 11 November 1999
//---------------------------------------------------------------------

#ifndef G4String_hh
#define G4String_hh 1

#include "G4Types.hh"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <string_view>

class G4String : public std::string
{
 public:
  /// @brief
  /// @deprecated Will be removed in future release
  enum caseCompare
  {
    exact,
    ignoreCase
  };

  /// @brief
  /// @deprecated Will be removed in future release
  enum stripType
  {
    leading,
    trailing,
    both
  };

  using std::string::string;
  using std::string::operator=;

  inline G4String() = default;
  inline G4String(const std::string&);
  inline G4String(const G4String&);
  inline G4String(std::string&&);
  inline G4String(G4String&&);
  inline G4String& operator=(const G4String&);
  inline G4String& operator=(G4String&&);

  /// @brief Implicitly convert to `const char*`
  /// @deprecated Will be removed in future releases for `std::string` compliance
  ///   If passing `G4String` to functions requiring `const char*`, use `std::string::c_str`
  ///   to explicitly convert. `G4String` also implicitly converts to `std::string_view`
  ///   to match the `std::string` interface.
  inline operator const char*() const;

  /// @brief Override of subscript operator for `int` to suppress C2666 errors with MSVC
  /// @deprecated Will be removed at the same time as `operator const char*` that requires it
  ///
  /// This override is required because of G4String's provision of an implicit conversion
  /// operator to `const char*`. Together with the subscript operator and C++'s built-in
  /// `operator[](const char*, int) operator, use of G4String::operator[] will trigger
  /// [MSVC error C2666](https://docs.microsoft.com/en-us/cpp/error-messages/compiler-errors-2/compiler-error-c2666?view=msvc-170)
  /// This is a known issue with mixing implicit conversion to `const char*` and subscript
  /// operators. Provision of the override with `int` argument is thus a workaround
  /// until the conversion operator is removed.
  inline reference operator[](int);

  /// @overload  
  inline const_reference operator[](int) const;

  /// @brief Deprecated function
  /// @deprecated Use `std::string::compare` or `G4StrUtil::icompare` instead
  [[deprecated("Use std::string::compare, or G4StrUtil::icompare for case-insensitive comparison")]]
  inline G4int compareTo(std::string_view, caseCompare mode = exact) const;

  /// @brief Deprecated function
  /// @deprecated Use `std::getline` plus `G4StrUtil::lstrip` instead
  [[deprecated("Use std::getline instead, plus G4StrUtil::lstrip if required")]]
  inline std::istream& readLine(std::istream&, G4bool skipWhite = true);
  
  /// @brief Deprecated function
  /// @deprecated Use `std::string::erase` instead
  [[deprecated("Use std::string::erase instead")]]
  inline G4String& remove(size_type);

  /// @brief Deprecated function
  /// @deprecated Use `G4StrUtil::contains` instead
  [[deprecated("Use G4StrUtil::contains instead")]]
  inline G4bool contains(const std::string&) const;

  /// @brief Deprecated function
  /// @deprecated Use `G4StrUtil::contains` instead 
  [[deprecated("Use G4StrUtil::contains instead")]]
  inline G4bool contains(char) const;

  /// @brief Deprecated function
  /// @deprecated Use `G4StrUtil` functions instead
  [[deprecated("Use G4StrUtil::{lstrip,rstrip,strip}_copy instead")]]
  [[nodiscard]] inline G4String strip(stripType strip_Type = trailing, char ch = ' ');

  /// @brief Deprecated function
  /// @deprecated Use `G4StrUtil` functions instead
  [[deprecated("Use G4StrUtil::to_lower/to_lower_copy instead")]]
  inline void toLower();

  /// @brief Deprecated function
  /// @deprecated Use `G4StrUtil` functions instead
  [[deprecated("Use G4StrUtil::to_upper/to_upper_copy instead")]]
  inline void toUpper();
};

/// @brief Query and manipulation functions for G4String
//
/// Additional free functions that are not part of the `std::string` interface as
/// of the minimum C++ standard supported by Geant4 (currently C++17).
///
/// @see `G4String`
namespace G4StrUtil
{
  /// @brief Convert string to lowercase
  /// @param[in, out] str the string to lowercase
  inline void to_lower(G4String& str);

  /// @brief Return lowercased copy of string
  /// @param[in] str the string to lowercase
  /// @return lowercased copy of `str`
  inline G4String to_lower_copy(G4String str);

  /// @brief Convert string to uppercase
  /// @param[in, out] str the string to uppercase
  inline void to_upper(G4String& str);

  /// @brief Return uppercase copy of string
  /// @param[in] str the string to upper case
  /// @return uppercased copy of `str`
  inline G4String to_upper_copy(G4String str);

  /// @brief Remove leading characters from string
  /// @param[in,out] str string to strip
  /// @param[in] ch character to remove
  /// @post `str` has any leading sequence of `ch` removed
  void lstrip(G4String& str, char ch = ' ');

  /// @brief Remove trailing characters from string
  /// @param[in,out] str string to strip
  /// @param[in] ch character to remove
  /// @post `str` has any trailing sequence of `ch` removed
  void rstrip(G4String& str, char ch = ' ');

  /// @brief Remove leading and trailing characters from string
  /// @param[in,out] str string to strip
  /// @param[in] ch character to remove
  /// @post `str` has any leading and trailing sequence of `ch` removed
  void strip(G4String& str, char ch = ' ');

  /// @brief Return copy of string with leading characters removed
  /// @param[in] str string to copy and strip
  /// @param[in] ch character to remove
  /// @return copy of `str` with any leading sequence of `ch` removed
  G4String lstrip_copy(G4String str, char ch = ' ');

  /// @brief Return copy of string with trailing characters removed
  /// @param[in] str string to copy and strip
  /// @param[in] ch character to remove
  /// @return copy of `str` with any trailing sequence of `ch` removed
  G4String rstrip_copy(G4String str, char ch = ' ');

  /// @brief Return copy of string with leading and trailing characters removed
  /// @param[in] str string to copy and strip
  /// @param[in] ch character to remove
  /// @return copy of `str` with any leading and trailing sequence of `ch` removed
  G4String strip_copy(G4String str, char ch = ' ');

  /// @brief Check if a string contains a given substring
  /// @param[in] str string to be checked
  /// @param[in] ss substring to check for
  /// @retval true if `str` contains `ss`
  /// @retval false otherwise
  inline G4bool contains(const G4String& str, std::string_view ss);

  /// @overload
  inline G4bool contains(const G4String& str, char ss);

  /// @overload
  inline G4bool contains(const G4String& str, const char* ss);

  /// @overload
  /// @note this overload is required to resolve ambiguity between the
  ///   signatures taking `std::string_view` and `const char*` substring
  ///   arguments. G4String currently provides implicit conversion to
  ///   `const char*`, which makes the calls ambiguous due to `std::string`'s
  ///   implicit conversion to `std::string_view`
  inline G4bool contains(const G4String& str, const G4String& ss);

  /// @brief Case insensitive comparison of two strings
  ///
  /// Converts both input arguments to lower case and returns the
  /// result of `std::string::compare` with these values.
  ///
  /// @param[in] lhs the first string in the comparison
  /// @param[in] rhs the second string in the comparison
  /// @return negative(positive) `G4int` if lowercased `lhs` appears
  ///   before(after) lowercased `rhs` in lexicographical order, zero if both
  ///   compare equivalent after lowercasing.
  inline G4int icompare(std::string_view lhs, std::string_view rhs);

  /// @brief Return true if a string starts with a given prefix
  /// @param[in] str string to be checked
  /// @param[in] ss prefix to check for
  /// @retval true if `str` starts with `ss`
  /// @retval false otherwise
  inline bool starts_with(const G4String& str, std::string_view ss);

  /// @overload
  inline bool starts_with(const G4String& str, G4String::value_type ss);

  /// @overload
  inline bool starts_with(const G4String& str, const char* ss);

  /// @overload
  /// @note this overload is required to resolve ambiguity between the
  ///   signatures taking `std::string_view` and `const char*` substring
  ///   arguments. G4String currently provides implicit conversion to
  ///   `const char*`, which makes the calls ambiguous due to `std::string`'s
  ///   implicit conversion to `std::string_view`
  inline bool starts_with(const G4String& str, const G4String& ss);

  /// @brief Return true if a string ends with a given suffix
  /// @param[in] str string to be checked
  /// @param[in] ss suffix to check for
  /// @retval true if `str` ends with `ss`
  /// @retval false otherwise
  inline bool ends_with(const G4String& str, std::string_view ss);

  /// @overload
  inline bool ends_with(const G4String& str, G4String::value_type ss);

  /// @overload
  inline bool ends_with(const G4String& str, const char* ss);

  /// @overload
  inline bool ends_with(const G4String& str, const G4String& ss);

  /// @brief Remove specified in-range characters from string
  /// 
  /// Equivalent to `std::string::erase(index, count)` with erasure only occuring
  /// if `index <= size()`. When `index > size()` the string is left unmodified.
  ///
  /// @deprecated It is strongly recommended to use `std::string::erase` if the
  ///   start index is already checked, or otherwise known, to be in range. Otherwise,
  ///   implement the index-size comparison instead of using this function.
  ///
  /// @param[in,out] str string to erase characters from
  /// @param[in] index position to start removal from
  /// @param[in] count number of characters to remove
  /// @post `str` is unchanged if `index > str.size()`, otherwise `str` has
  ///   `min(count, str.size() - index)` characters removed starting at `index`
  inline void safe_erase(G4String& str, G4String::size_type index = 0,
                         G4String::size_type count = G4String::npos);

  /// @brief Read characters into a G4String from an input stream until end-of-line
  ///
  /// @deprecated It is strongly recommended to use `std::getline` instead of this
  ///   function, plus `G4StrUtil::lstrip` if leading whitespace removal is required. 
  /// 
  /// @param[in] is input stream to read from
  /// @param[in,out] str string to read into
  /// @param[in] skipWhite if true, discard leading whitespace from `is`
  /// @return `is`
  inline std::istream& readline(std::istream& is, G4String& str, G4bool skipWhite = true);
}  // namespace G4StrUtil

#include "G4String.icc"

#endif
