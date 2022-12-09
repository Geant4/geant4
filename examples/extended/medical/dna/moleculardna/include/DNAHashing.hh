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
//
//
//  MolecularDNAHashing.hh
//  Geant4
//
//  Created by Mathieu Karamitros
//
//

#ifndef MolecularDNAHashing_h
#define MolecularDNAHashing_h

#include <string>
#include <type_traits>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace G4::hashing
{
  namespace crc32
  {
    //--------------------------------------------------------------------------
    // Generate CRC lookup table
    //--------------------------------------------------------------------------

    template <unsigned long long c, int k = 8>
    struct f : f<((c & 1) ? 0xedb88320 : 0) ^ (c >> 1), k - 1>
    {};
    template <unsigned long long c>
    struct f<c, 0>
    {
      enum
      {
        value = c
      };
    };

#define _A_(x) _B_(x) _B_(x + 128llu)
#define _B_(x) _C_(x) _C_(x + 64llu)
#define _C_(x) _D_(x) _D_(x + 32llu)
#define _D_(x) _E_(x) _E_(x + 16llu)
#define _E_(x) _F_(x) _F_(x + 8llu)
#define _F_(x) _G_(x) _G_(x + 4llu)
#define _G_(x) _H_(x) _H_(x + 2llu)
#define _H_(x) _I_(x) _I_(x + 1llu)
#define _I_(x) f<x>::value,

    static constexpr unsigned fCrc_table[] = { _A_(0llu) };

    //--------------------------------------------------------------------------
    // CRC32 ALGO
    //--------------------------------------------------------------------------

    //--------------------------------
    // RUN TIME ALGO
    uint32_t Hash(const char* str, size_t len);

    //--------------------------------------------------------------------------
    // Details for compilation time
    namespace detail
    {
      // CRC32 Table (zlib polynomial)
      template <size_t idx>
      constexpr uint32_t Combine_crc32(const char* str, uint32_t part)
      {
        return (part >> 8) ^ fCrc_table[(part ^ str[idx]) & 0x000000FF];
      }

      // recursion
      template <size_t idx>
      constexpr uint32_t Crc32(const char* str)
      {
        return Combine_crc32<idx>(str, Crc32<idx - 1>(str));
      }

      // stop-recursion
      template <>
      constexpr uint32_t Crc32<size_t(-1)>(const char* /*str*/)
      {
        return 0xFFFFFFFF;
      }
    }  // namespace detail

    //--------------------------------------------------------------------------
    // use COMPILATION TIME CRC32 as a hash function
    template <size_t len>
    constexpr uint32_t Hash(const char (&str)[len])
    {
      return detail::Crc32<len - 2>(str) ^ 0xFFFFFFFF;
    }

    //--------------------------------------------------------------------------
    // std::string hasher for convenience
    uint32_t Hash(const std::string& str);
    // in window
    //warning C4717: 'G4::hashing::crc32::crc32_hasher<G4String>::operator()': recursive on all control paths, function will cause runtime stack overflow
    /*
    template <typename T>
    struct crc32_hasher
    {
      template <size_t len>
      uint32_t constexpr operator()(const char (&str)[len]) const
      {
        return hash(str);
      }

      uint32_t operator()(T str) const { return hash(str); }
    };

    template <typename T>
    uint32_t constexpr hash(T&& t)
    {
      return crc32_hasher<typename std::decay<T>::type>()(std::forward<T>(t));
    }
    */
  }  // namespace crc32

  //--------------------------------------------------------------------------
  // fnv ALGO
  //--------------------------------------------------------------------------

  namespace fnv
  {
    // Fowler / Noll / Vo (fnv) Hash - 1a
    // magic numbers from http://www.isthe.com/chongo/tech/comp/fnv/
#define fnv_offset_basis 0xcbf29ce484222325U
#define fnv_prime 0x100000001b3U  // 64 bits
    // 224 + 28 + 0x93

    //------------------------------------------------------------------------
    // Run time
    size_t Hash(const std::string& str);

    //--------------------------------------------------------------------------
    // Details for compilation time
    namespace detail
    {
      // recursion
      template <int idx>
      constexpr size_t Fnv(const char* str, size_t len, size_t hash)
      {
        return Fnv<idx - 1>(str, len, (hash ^ (str[len - idx])) * fnv_prime);
      }

      // stop-recursion
      template <>
      constexpr size_t Fnv<-1>(const char*, size_t, size_t hash)
      {
        return hash;
      }
    }  // namespace detail

    //------------------------------------------------------------------------
    // COMPILATION TIME
    template <size_t len>
    constexpr size_t Hash(const char (&str)[len])
    {
      return detail::Fnv<len - 2>(str, len - 2, fnv_offset_basis);
    }

    //------------------------------------------------------------------------
    template <class>
    struct fnv_hasher;

    template <typename T>
    struct fnv_hasher
    {
      template <size_t len>
      size_t constexpr operator()(const char (&str)[len]) const
      {
        return Hash(str);
      }

      size_t operator()(const T& str) const { return Hash(str); }
    };

    //      template<>
    //      size_t fnv_hasher<std::string>::operator()(const std::string& str)
    //      const{
    //          return hash(str.c_str());
    //      }

    template <typename T>
    size_t constexpr Hash(T&& t)
    {
      return fnv_hasher<typename std::decay<T>::type>()(std::forward<T>(t));
    }
  }  // namespace fnv

  //--------------------------------------------------------------------------
  // LARSON
  //--------------------------------------------------------------------------

  namespace larson
  {
    //------------------------------------------------------------------------
    // Details for compilation time
    namespace detail
    {
      // recursion
      template <int idx>
      constexpr size_t Larson(const char* str, size_t len, size_t hash)
      {
        return Larson<idx - 1>(str, len + 1, hash * 101 + str[len]);
      }

      // stop-recursion
      template <>
      constexpr size_t Larson<-1>(const char*, size_t, size_t hash)
      {
        return hash;
      }
    }  // namespace detail

    //------------------------------------------------------------------------
    // COMPILATION TIME
    template <size_t len>
    constexpr size_t Cthash(const char (&str)[len], unsigned int seed = 0)
    {
      return detail::Larson<len - 2>(str, 0, seed);
    }

    //------------------------------------------------------------------------
    // RUN TIME
    size_t Hash(const char* str, unsigned int seed = 0);

    size_t Hash(std::string str, unsigned int seed = 0);
  }  // namespace larson

  //--------------------------------------------------------------------------
  // DEFAULT STRING HASHER
  //--------------------------------------------------------------------------

  template <typename T>
  struct hasher
  {
    std::size_t constexpr operator()(char const* input) const
    {
      return *input ? static_cast<std::size_t>(*input) + 33 * (*this)(input + 1)
                    : 5381;
    }

    std::size_t operator()(const std::string& str) const
    {
      return (*this)(str.c_str());
    }
  };

  template <typename T>
  std::size_t constexpr Hash(T&& t)
  {
    return hasher<typename std::decay<T>::type>()(std::forward<T>(t));
  }

  // not necessary
  template <>
  std::size_t constexpr Hash<const char*>(const char*&& str)
  {
    return hasher<std::string>()(str);
  }

  inline namespace literals
  {
    std::size_t constexpr operator"" _hash(const char* s, size_t)
    {
      return hasher<std::string>()(s);
    }
  }  // namespace literals
}  // namespace G4::hashing

using namespace G4::hashing::literals;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
