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
// $Id: G4SurfBits.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4SurfBits
//
// Class description:
//
// This class provides a simple container of bits, to be used for
// optimization of tessellated surfaces (G4TessellatedSolid).
// Each bit can be set and tested via the functions SetBitNumber and
// TestBitNumber.
// The default value of all bits is false.
// The size of the container is automatically extended when a bit
// number is either set or tested.  To reduce the memory size of the
// container use the Compact function, this will discard the memory
// occupied by the upper bits that are 0.
//

// History:
// 19.10.12 Marek Gayer, created and adapted from original implementation
//                       of Root's TBits class by P.Canal
// --------------------------------------------------------------------

#ifndef G4SurfBits_HH
#define G4SurfBits_HH

#include <cstring>

#include "G4Types.hh"

class G4SurfBits
{
  public:

    G4SurfBits(unsigned int nbits = 0);
    G4SurfBits(const G4SurfBits&);
    G4SurfBits& operator=(const G4SurfBits&);
   ~G4SurfBits();

    //----- Bit manipulation
    void ResetAllBits(G4bool value=false);  // if value=1 set all bits to 1
    void ResetBitNumber(unsigned int bitnumber);
    void SetBitNumber(unsigned int bitnumber, G4bool value = true);
    G4bool TestBitNumber(unsigned int bitnumber) const;

    //----- Accessors and operator
    G4bool operator[](unsigned int bitnumber) const;

    //----- Optimized setters
    // Each of these will replace the contents of the receiver with the
    // bitvector in the parameter array. The number of bits is changed
    // to nbits. If nbits is smaller than fNBits, the receiver will NOT
    // be compacted.
    void   set(unsigned int nbits, const char *array);
    void   set(unsigned int nbits, const G4int *array);

    //----- Optimized getters
    // Each of these will replace the contents of the parameter array with the
    // bits in the receiver.  The parameter array must be large enough to hold
    // all of the bits in the receiver.
    // Note on semantics: any bits in the parameter array that go beyond the
    // number of the bits in the receiver will have an unspecified value. For
    // example, if you call Get(Int*) with an array of one integer and the
    // G4SurfBits object has less than 32 bits, then the remaining bits in the
    // integer will have an unspecified value.
    void   Get(char *array) const;
    void   Get(G4int *array) const;

    //----- Utilities
    void    Clear();
    void    Compact();               // Reduce the space used.

    unsigned int  GetNbits()      const { return fNBits; }
    unsigned int  GetNbytes()     const { return fNBytes; }

    void    Print() const;  // to show the list of active bits
    void    Output(std::ostream &) const;

  protected:

    void ReserveBytes(unsigned int nbytes);

  public:

    unsigned char *fAllBits;       // [fNBytes] array of UChars

  protected:

    unsigned int   fNBits;         // Highest bit set + 1
    unsigned int   fNBytes;        // Number of UChars in fAllBits
};

// inline functions...

inline void G4SurfBits::SetBitNumber(unsigned int bitnumber, G4bool value)
{
  // set bit number 'bitnumber' to be value
  if (bitnumber >= fNBits) {
    unsigned int new_size = (bitnumber/8) + 1;
    if (new_size > fNBytes) {
      if (new_size < 100 * 1024 * 1024)
        new_size *= 2;
      unsigned char *old_location = fAllBits;
      fAllBits = new unsigned char[new_size];
      std::memcpy(fAllBits,old_location,fNBytes);
      std::memset(fAllBits+fNBytes ,0, new_size-fNBytes);
      fNBytes = new_size;
      delete [] old_location;
    }
    fNBits = bitnumber+1;
  }
  unsigned int  loc = bitnumber/8;
  unsigned char bit = bitnumber%8;
  if (value)
    fAllBits[loc] |= (1<<bit);
  else
    fAllBits[loc] &= (0xFF ^ (1<<bit));
}

inline G4bool G4SurfBits::TestBitNumber(unsigned int bitnumber) const
{
  // Return the current value of the bit

  if (bitnumber >= fNBits) return false;
  unsigned int  loc = bitnumber/8;
  unsigned char value = fAllBits[loc];
  unsigned char bit = bitnumber%8;
  G4bool  result = (value & (1<<bit)) != 0;
  return result;
  // short: return 0 != (fAllBits[bitnumber/8] & (1<< (bitnumber%8)));
}

inline void G4SurfBits::ResetBitNumber(unsigned int bitnumber)
{
  SetBitNumber(bitnumber,false);
}

inline G4bool G4SurfBits::operator[](unsigned int bitnumber) const
{
  return TestBitNumber(bitnumber);
}

#endif
