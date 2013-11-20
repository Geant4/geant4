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
// 19.10.12 Marek Gayer
//          Created from original implementation in ROOT (TBits)
// --------------------------------------------------------------------

#include "UBits.hh"
#include <stdio.h>
//______________________________________________________________________________
UBits::UBits(unsigned int nBits) : fNBits(nBits)
{
  // UBits constructor.  All bits set to 0

  if (fNBits <= 0) fNBits = 0;
  fNBytes  = fNBits ? ((fNBits - 1) / 8) + 1 : 1;
  fAllBits = new unsigned char[fNBytes];
  // this is redundant only with libNew
  std::memset(fAllBits, 0, fNBytes);
}

//______________________________________________________________________________
UBits::UBits(const UBits& original) : fNBits(original.fNBits),
  fNBytes(original.fNBytes)
{
  // UBits copy constructor

  fAllBits = new unsigned char[fNBytes];
  std::memcpy(fAllBits, original.fAllBits, fNBytes);

}


//______________________________________________________________________________
UBits& UBits::operator=(const UBits& rhs)
{
  // UBits assignment operator
   // Check assignment to self    
 if (this == &rhs)  { return *this; }

    //      TObject::operator=(rhs);
    fNBits   = rhs.fNBits;
    fNBytes  = rhs.fNBytes;
    delete [] fAllBits;
    if (fNBytes != 0) {
      fAllBits = new unsigned char[fNBytes];
      std::memcpy(fAllBits,rhs.fAllBits,fNBytes);
    } else {
      fAllBits = 0;
    }
  return *this;
}

//______________________________________________________________________________
UBits::~UBits()
{
  // UBits destructor

  delete [] fAllBits;
}

//______________________________________________________________________________
void UBits::Clear()
{
  // Clear the value.

  delete [] fAllBits;
  fAllBits = 0;
  fNBits   = 0;
  fNBytes  = 0;
}

//______________________________________________________________________________
void UBits::Compact()
{
  // Reduce the storage used by the object to a minimun

  if (!fNBits || !fAllBits) return;
  unsigned int needed;
  for (needed = fNBytes - 1;
       needed > 0 && fAllBits[needed] == 0;)
  {
    needed--;
  };
  needed++;

  if (needed != fNBytes)
  {
    unsigned char* old_location = fAllBits;
    fAllBits = new unsigned char[needed];

    std::memcpy(fAllBits, old_location, needed);
    delete [] old_location;

    fNBytes = needed;
    fNBits = 8 * fNBytes;
  }
}

/*
//______________________________________________________________________________
unsigned int UBits::CounUBits(unsigned int startBit) const
{
  // Return number of bits set to 1 starting at bit startBit

  static const int nBitsCached[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

    unsigned int i,count = 0;
    if (startBit == 0) {
      for(i=0; i<fNBytes; i++) {
        count += nBitsCached[fAllBits[i]];
      }
      return count;
    }
    if (startBit >= fNBits) return count;
    unsigned int startByte = startBit/8;
    unsigned int ibit = startBit%8;
    if (ibit) {
      for (i=ibit;i<8;i++) {
        if (fAllBits[startByte] & (1<<ibit)) count++;
      }
      startByte++;
    }
    for(i=startByte; i<fNBytes; i++) {
      count += nBitsCached[fAllBits[i]];
    }
    return count;
}
*/

/*
//______________________________________________________________________________
void UBits::DoAndEqual(const UBits& rhs)
{
  // Execute (*this) &= rhs;
  // Extra bits in rhs are ignored
  // Missing bits in rhs are assumed to be zero.

  unsigned int min = (fNBytes<rhs.fNBytes) ? fNBytes : rhs.fNBytes;
  for(unsigned int i=0; i<min; ++i) {
    fAllBits[i] &= rhs.fAllBits[i];
  }
  if (fNBytes>min) {
    std::memset(&(fAllBits[min]),0,fNBytes-min);
  }
}
*/

/*
//______________________________________________________________________________
void UBits::DoOrEqual(const UBits& rhs)
{
  // Execute (*this) &= rhs;
  // Extra bits in rhs are ignored
  // Missing bits in rhs are assumed to be zero.

  unsigned int min = (fNBytes<rhs.fNBytes) ? fNBytes : rhs.fNBytes;
  for(unsigned int i=0; i<min; ++i) {
    fAllBits[i] |= rhs.fAllBits[i];
  }
}

//______________________________________________________________________________
void UBits::DoXorEqual(const UBits& rhs)
{
  // Execute (*this) ^= rhs;
  // Extra bits in rhs are ignored
  // Missing bits in rhs are assumed to be zero.

  unsigned int min = (fNBytes<rhs.fNBytes) ? fNBytes : rhs.fNBytes;
  for(unsigned int i=0; i<min; ++i) {
    fAllBits[i] ^= rhs.fAllBits[i];
  }
}

//______________________________________________________________________________
void UBits::DoFlip()
{
  // Execute ~(*this)

  for(unsigned int i=0; i<fNBytes; ++i) {
    fAllBits[i] = ~fAllBits[i];
  }
  // NOTE: out-of-bounds bit were also flipped!
}

//______________________________________________________________________________
void UBits::DoLeftShift(unsigned int shift)
{
  // Execute the left shift operation.

  if (shift==0) return;
  const unsigned int wordshift = shift / 8;
  const unsigned int offset = shift % 8;
  if (offset==0) {
    for(unsigned int n = fNBytes - 1; n >= wordshift; --n) {
      fAllBits[n] = fAllBits[ n - wordshift ];
    }
  } else {
    const unsigned int sub_offset = 8 - offset;
    for(unsigned int n = fNBytes - 1; n > wordshift; --n) {
      fAllBits[n] = (fAllBits[n - wordshift] << offset) |
        (fAllBits[n - wordshift - 1] >> sub_offset);
    }
    fAllBits[wordshift] = fAllBits[0] << offset;
  }
  std::memset(fAllBits,0,wordshift);
}

//______________________________________________________________________________
void UBits::DoRightShift(unsigned int shift)
{
  // Execute the left shift operation.

  if (shift==0) return;
  const unsigned int wordshift = shift / 8;
  const unsigned int offset = shift % 8;
  const unsigned int limit = fNBytes - wordshift - 1;

  if (offset == 0)
    for (unsigned int n = 0; n <= limit; ++n)
      fAllBits[n] = fAllBits[n + wordshift];
  else
  {
    const unsigned int sub_offset = 8 - offset;
    for (unsigned int n = 0; n < limit; ++n)
      fAllBits[n] = (fAllBits[n + wordshift] >> offset) |
      (fAllBits[n + wordshift + 1] << sub_offset);
    fAllBits[limit] = fAllBits[fNBytes-1] >> offset;
  }

  std::memset(&(fAllBits[limit + 1]),0, fNBytes - limit - 1);
}
*/

/*
//______________________________________________________________________________
unsigned int UBits::FirstNullBit(unsigned int startBit) const
{
  // Return position of first null bit (starting from position 0 and up)

  static const int fbits[256] = {
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,
    0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,8};

    unsigned int i;
    if (startBit == 0) {
      for(i=0; i<fNBytes; i++) {
        if (fAllBits[i] != 255) return 8*i + fbits[fAllBits[i]];
      }
      return fNBits;
    }
    if (startBit >= fNBits) return fNBits;
    unsigned int startByte = startBit/8;
    unsigned int ibit = startBit%8;
    if (ibit) {
      for (i=ibit;i<8;i++) {
        if ((fAllBits[startByte] & (1<<i)) == 0) return 8*startByte+i;
      }
      startByte++;
    }
    for(i=startByte; i<fNBytes; i++) {
      if (fAllBits[i] != 255) return 8*i + fbits[fAllBits[i]];
    }
    return fNBits;
}

//______________________________________________________________________________
unsigned int UBits::FirstSetBit(unsigned int startBit) const
{
  // Return position of first non null bit (starting from position 0 and up)

  static const int fbits[256] = {
    8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};

    unsigned int i;
    if (startBit == 0) {
      for(i=0; i<fNBytes; i++) {
        if (fAllBits[i] != 0) return 8*i + fbits[fAllBits[i]];
      }
      return fNBits;
    }
    if (startBit >= fNBits) return fNBits;
    unsigned int startByte = startBit/8;
    unsigned int ibit = startBit%8;
    if (ibit) {
      for (i=ibit;i<8;i++) {
        if ((fAllBits[startByte] & (1<<i)) != 0) return 8*startByte+i;
      }
      startByte++;
    }
    for(i=startByte; i<fNBytes; i++) {
      if (fAllBits[i] != 0) return 8*i + fbits[fAllBits[i]];
    }
    return fNBits;
}
*/

//______________________________________________________________________________
void UBits::Output(std::ostream& os) const
{
  // Print the value to the std::ostream
  for (unsigned int i = 0; i < fNBytes; ++i)
  {
    unsigned char val = fAllBits[fNBytes - 1 - i];
    for (unsigned int j = 0; j < 8; ++j)
    {
      os << (bool)(val & 0x80);
      val <<= 1;
    }
  }
}

//______________________________________________________________________________
void UBits::Print() const
{
  // Print the list of active bits
  int count = 0;
  for (unsigned int i = 0; i < fNBytes; ++i)
  {
    unsigned char val = fAllBits[i];
    for (unsigned int j = 0; j < 8; ++j)
    {
      if (val & 1) printf(" bit:%4d = 1\n", count);
      count++;
      val = val >> 1;
    }
  }
}

//______________________________________________________________________________
void UBits::ResetAllBits(bool value)
{
  if (fAllBits) std::memset(fAllBits, value ? 0xFF : 0, fNBytes);
}

//______________________________________________________________________________
void UBits::ReserveBytes(unsigned int nbytes)
{
  // Reverse each bytes.

  if (nbytes > fNBytes)
  {
    // do it in this order to remain exception-safe.
    unsigned char* newBits = new unsigned char[nbytes];
    delete[] fAllBits;
    fNBytes = nbytes;
    fAllBits = newBits;
  }
}

//______________________________________________________________________________
void UBits::Set(unsigned int nBits, const char* array)
{
  // Set all the bytes
  unsigned int nbytes = (nBits + 7) >> 3;

  ReserveBytes(nbytes);

  fNBits = nBits;
  std::memcpy(fAllBits, array, nbytes);
}

//______________________________________________________________________________
void UBits::Get(char* array) const
{
  // Copy all the byes.
  std::memcpy(array, fAllBits, (fNBits + 7) >> 3);
}

// If we are on a little endian machine, a bitvector represented using
// any integer type is identical to a bitvector represented using bytes. -- FP.


void UBits::Set(unsigned int nBits, const int* array)
{
  // Set all the bytes.

  Set(nBits, (const char*)array);
}

void UBits::Get(int* array) const
{
  // Get all the bytes.

  Get((char*)array);
}



/*
bool UBits::operator==(const UBits &other) const
{
  // Compare object.

  if (fNBits == other.fNBits) {
    return !memcmp(fAllBits, other.fAllBits, (fNBits+7)>>3);
  } else if (fNBits <  other.fNBits) {
    return !memcmp(fAllBits, other.fAllBits, (fNBits+7)>>3) && other.FirstSetBit(fNBits) == other.fNBits;
  } else {
    return !memcmp(fAllBits, other.fAllBits, (other.fNBits+7)>>3) && FirstSetBit(other.fNBits) == fNBits;
  }
}
*/
