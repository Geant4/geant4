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

// Author: Marek Gayer (CERN), 19.10.2012 - Created and adapted from original
//                                          implementation of Root/TBits class.
// --------------------------------------------------------------------
#ifndef G4SURFBITS_HH
#define G4SURFBITS_HH

#include <cstring>

#include "G4Types.hh"

/**
 * @brief G4SurfBits provides a simple container of bits, to be used for
 * optimization of tessellated surfaces. The size of the container is
 * automatically extended when a bit number is either set or tested.
 */

class G4SurfBits
{
  public:

    /**
     * Constructor given the number of bits.
     *  @param[in] nbits The number of bits.
     */
    G4SurfBits(unsigned int nbits = 0);

    /**
     * Destructor. Clears all allocated bits.
     */
   ~G4SurfBits();

    /**
     * Copy constructor and assignment operator.
     */
    G4SurfBits(const G4SurfBits&);
    G4SurfBits& operator=(const G4SurfBits&);

    /**
     * Methods for bit manipulation.
     */
    void ResetAllBits(G4bool value = false);  // if value=1 set all bits to 1
    inline void ResetBitNumber(unsigned int bitnumber);
    inline void SetBitNumber(unsigned int bitnumber, G4bool value = true);
    inline G4bool TestBitNumber(unsigned int bitnumber) const;

    /**
     * Accessor operator.
     */
    inline G4bool operator[](unsigned int bitnumber) const;

    /**
     * Optimized setters. Each of these will replace the contents of the
     * receiver with the bitvector in the parameter array. The number of bits
     * is changed to 'nbits'. If nbits is smaller than fNBits, the receiver
     * will NOT be compacted.
     */
    void set(unsigned int nbits, const char* array);
    void set(unsigned int nbits, const G4int* array);

    /**
     * Optimized getters. Each of these will replace the contents of the
     * parameter array with the bits in the receiver. The parameter array must
     * be large enough to hold all of the bits in the receiver.
     * Note on semantics: any bits in the parameter array that go beyond the
     * number of the bits in the receiver will have an unspecified value. For
     * example, if calling Get(Int*) with an array of one integer and the
     * G4SurfBits object has less than 32 bits, then the remaining bits in the
     * integer will have an unspecified value.
     */
    void Get(char* array) const;
    void Get(G4int* array) const;

    /**
     * Utilities to clear or reduce the space used.
     */
    void Clear();
    void Compact();               // Reduce the space used.

    /**
     * Accessors.
     */
    inline unsigned int GetNbits() const { return fNBits; }
    inline unsigned int GetNbytes() const { return fNBytes; }

    /**
     * Logging functions.
     */
    void Print() const;  // to show the list of active bits
    void Output(std::ostream &) const;

  public:

    unsigned char* fAllBits = nullptr; // [fNBytes] array of UChars

  private:

    void ReserveBytes(unsigned int nbytes);

  private:

    unsigned int fNBits;         // Highest bit set + 1
    unsigned int fNBytes;        // Number of UChars in fAllBits
};

// inline functions...

#include "G4SurfBits.icc"

#endif
