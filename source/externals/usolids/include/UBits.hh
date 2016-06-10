//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UBits
//
// Class description:
//
// Container of bits
//
// This class provides a simple container of bits.
// Each bit can be set and tested via the functions SetBitNumber and
// TestBitNumber.
// The default value of all bits is false.
// The size of the container is automatically extended when a bit
// number is either set or tested.  To reduce the memory size of the
// container use the Compact function, this will discard the memory
// occupied by the upper bits that are 0.
//
// Created for UTessellatedSolid
//
// 19.10.12 Marek Gayer
//          Created from original implementation in ROOT (TBits)
// --------------------------------------------------------------------

#ifndef UBits_HH
#define UBits_HH

#include <cstring>
#include <ostream>

class UBits
{

  public:
    UBits(unsigned int nbits = 0);
    UBits(const UBits&);
    UBits& operator=(const UBits& rhs);
    virtual ~UBits();

    //----- bit manipulation
    //----- (note the difference with TObject's bit manipulations)
    void ResetAllBits(bool value = false); // if value=1 set all bits to 1
    void ResetBitNumber(unsigned int bitnumber);
    void SetBitNumber(unsigned int bitnumber, bool value = true);
    bool TestBitNumber(unsigned int bitnumber) const;

    //----- Accessors and operator
    bool operator[](unsigned int bitnumber) const;

    /*
    UBits& operator&=(const UBits& rhs) { DoAndEqual(rhs); return *this; }
    UBits& operator|=(const UBits& rhs) {  DoOrEqual(rhs); return *this; }
    UBits& operator^=(const UBits& rhs) { DoXorEqual(rhs); return *this; }
    UBits& operator<<=(unsigned int rhs) { DoLeftShift(rhs); return *this; }
    UBits& operator>>=(unsigned int rhs) { DoRightShift(rhs); return *this; }
    UBits  operator<<(unsigned int rhs) { return UBits(*this)<<= rhs; }
    UBits  operator>>(unsigned int rhs) { return UBits(*this)>>= rhs; }
    UBits  operator~() { UBits res(*this); res.DoFlip(); return res; }
    */

    //----- Optimized setters
    // Each of these will replace the contents of the receiver with the bitvector
    // in the parameter array.  The number of bits is changed to nbits.  If nbits
    // is smaller than fNBits, the receiver will NOT be compacted.

    void   Set(unsigned int nbits, const char* array);
//  void   Set(unsigned int nbits, const unsigned char *array) { Set(nbits, (const char*)array); }
//  void   Set(unsigned int nbits, const short *array);
    //void   Set(unsigned int nbits, const unsigned short *array) { Set(nbits, (const short*)array); }
    void   Set(unsigned int nbits, const int* array);
//  void   Set(unsigned int nbits, const unsigned int *array) { Set(nbits, (const int*)array); }

    //----- Optimized getters
    // Each of these will replace the contents of the parameter array with the
    // bits in the receiver.  The parameter array must be large enough to hold
    // all of the bits in the receiver.
    // Note on semantics: any bits in the parameter array that go beyond the
    // number of the bits in the receiver will have an unspecified value.  For
    // example, if you call Get(Int*) with an array of one integer and the UBits
    // object has less than 32 bits, then the remaining bits in the integer will
    // have an unspecified value.
    void   Get(char* array) const;
//  void   Get(unsigned char *array) const { Get((char*)array); }
//  void   Get(short *array) const;
//  void   Get(unsigned short *array) const { Get((short*)array); }
    void   Get(int* array) const;
//  void   Get(unsigned int *array) const { Get((int*)array); }

    //----- Utilities
    void    Clear();
    void    Compact();               // Reduce the space used.


    unsigned int  GetNbits()      const
    {
      return fNBits;
    }
    unsigned int  GetNbytes()     const
    {
      return fNBytes;
    }

    /*
    unsigned int  CounUBits(unsigned int startBit=0)     const ;  // return number of bits set to 1
    unsigned int  FirstNullBit(unsigned int startBit=0)  const;
    unsigned int  FirstSetBit(unsigned int startBit=0)   const;
    */

//  bool  operator==(const UBits &other) const;
//  bool  operator!=(const UBits &other) const { return !(*this==other); }

    void    Print() const;  // to show the list of active bits
    void    Output(std::ostream&) const;

  protected:
    void ReserveBytes(unsigned int nbytes);

    /*
    void DoAndEqual(const UBits& rhs);
    void DoOrEqual (const UBits& rhs);
    void DoXorEqual(const UBits& rhs);
    void DoLeftShift(unsigned int shift);
    void DoRightShift(unsigned int shift);
    void DoFlip();
    */

  protected:
    unsigned int   fNBits;         // Highest bit set + 1
    unsigned int   fNBytes;        // Number of UChars in fAllBits

  public:
    unsigned char* fAllBits;       //[fNBytes] array of UChars

};

/*
inline UBits operator&(const UBits& lhs, const UBits& rhs)
{
  UBits result(lhs);
  result &= rhs;
  return result;
}

inline UBits operator|(const UBits& lhs, const UBits& rhs)
{
  UBits result(lhs);
  result |= rhs;
  return result;
}

inline UBits operator^(const UBits& lhs, const UBits& rhs)
{
  UBits result(lhs);
  result ^= rhs;
  return result;
}

inline std::ostream &operator<<(std::ostream& os, const UBits& rhs)
{
  rhs.Output(os); return os;
}
*/

// inline functions...

inline void UBits::SetBitNumber(unsigned int bitnumber, bool value)
{
  // Set bit number 'bitnumber' to be value
  if (bitnumber >= fNBits)
  {
    unsigned int new_size = (bitnumber / 8) + 1;
    if (new_size > fNBytes)
    {
      if (new_size < 100 * 1024 * 1024)
        new_size *= 2;
      unsigned char* old_location = fAllBits;
      fAllBits = new unsigned char[new_size];
      std::memcpy(fAllBits, old_location, fNBytes);
      std::memset(fAllBits + fNBytes , 0, new_size - fNBytes);
      fNBytes = new_size;
      delete [] old_location;
    }
    fNBits = bitnumber + 1;
  }
  unsigned int  loc = bitnumber / 8;
  unsigned char bit = bitnumber % 8;
  if (value)
    fAllBits[loc] |= (1 << bit);
  else
    fAllBits[loc] &= (0xFF ^ (1 << bit));
}

inline bool UBits::TestBitNumber(unsigned int bitnumber) const
{
  // Return the current value of the bit

  if (bitnumber >= fNBits) return false;
  unsigned int  loc = bitnumber / 8;
  unsigned char value = fAllBits[loc];
  unsigned char bit = bitnumber % 8;
  bool  result = (value & (1 << bit)) != 0;
  return result;
  // short: return 0 != (fAllBits[bitnumber/8] & (1<< (bitnumber%8)));
}

inline void UBits::ResetBitNumber(unsigned int bitnumber)
{
  SetBitNumber(bitnumber, false);
}

inline bool UBits::operator[](unsigned int bitnumber) const
{
  return TestBitNumber(bitnumber);
}

#endif
