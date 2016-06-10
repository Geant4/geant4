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
// $Id: G4SurfBits.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4SurfBits implementation
//
// History:
// 19.10.12 Marek Gayer, created and adapted from original implementation
//                       of Root's TBits class by P.Canal
// --------------------------------------------------------------------

#include "G4SurfBits.hh"
#include "G4ios.hh"

//______________________________________________________________________________
G4SurfBits::G4SurfBits(unsigned int nBits) : fNBits(nBits)
{
  // G4SurfBits constructor.  All bits set to 0

  if (fNBits <= 0) fNBits = 0;
  fNBytes  = fNBits ? ((fNBits-1)/8) + 1 : 1;
  fAllBits = new unsigned char[fNBytes];
  // this is redundant only with libNew
  std::memset(fAllBits,0,fNBytes);
}

//______________________________________________________________________________
G4SurfBits::G4SurfBits(const G4SurfBits &original) : fNBits(original.fNBits),
  fNBytes(original.fNBytes)
{
  // G4SurfBits copy constructor

  fAllBits = new unsigned char[fNBytes];
  std::memcpy(fAllBits,original.fAllBits,fNBytes);
}

//______________________________________________________________________________
G4SurfBits& G4SurfBits::operator=(const G4SurfBits& rhs)
{
  // G4SurfBits assignment operator
  if (this != &rhs) {
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
  }
  return *this;
}

//______________________________________________________________________________
G4SurfBits::~G4SurfBits()
{
  // G4SurfBits destructor

  delete [] fAllBits;
}

//______________________________________________________________________________
void G4SurfBits::Clear()
{
  // Clear the value.

  delete [] fAllBits;
  fAllBits = 0;
  fNBits   = 0;
  fNBytes  = 0;
}

//______________________________________________________________________________
void G4SurfBits::Compact()
{
  // Reduce the storage used by the object to a minimun

  if (!fNBits || !fAllBits) return;
  unsigned int needed;
  for(needed=fNBytes-1;
    needed > 0 && fAllBits[needed]==0; ) { needed--; };
    needed++;

  if (needed!=fNBytes) {
    unsigned char *old_location = fAllBits;
    fAllBits = new unsigned char[needed];

    std::memcpy(fAllBits,old_location,needed);
    delete [] old_location;

    fNBytes = needed;
    fNBits = 8*fNBytes;
  }
}

//______________________________________________________________________________
void G4SurfBits::Output(std::ostream &os) const
{
  // Print the value to the std::ostream
  for(unsigned int i=0; i<fNBytes; ++i) {
    unsigned char val = fAllBits[fNBytes - 1 - i];
    for (unsigned int j=0; j<8; ++j) {
      os << (G4bool)(val&0x80);
      val <<= 1;
    }
  }
}

//______________________________________________________________________________
void G4SurfBits::Print() const
{
  // Print the list of active bits
  G4int count = 0;
  for(unsigned int i=0; i<fNBytes; ++i) {
    unsigned char val = fAllBits[i];
    for (unsigned int j=0; j<8; ++j) {
      if (val & 1) G4cout << " bit:" << count << " = 1" << G4endl;
      count++;
      val = val >> 1;
    }
  }
}

//______________________________________________________________________________
void G4SurfBits::ResetAllBits(G4bool value)
{
  if (fAllBits) std::memset(fAllBits, value ? 0xFF : 0,fNBytes);
}

//______________________________________________________________________________
void G4SurfBits::ReserveBytes(unsigned int nbytes)
{
  // Reverse each bytes.

  if (nbytes > fNBytes) {
    // do it in this order to remain exception-safe.
    unsigned char *newBits=new unsigned char[nbytes];
    delete[] fAllBits;
    fNBytes=nbytes;
    fAllBits=newBits;
  }
}

//______________________________________________________________________________
void G4SurfBits::set(unsigned int nBits, const char *array)
{
  // set all the bytes
  unsigned int nbytes=(nBits+7)>>3;

  ReserveBytes(nbytes);

  fNBits=nBits;
  std::memcpy(fAllBits, array, nbytes);
}

//______________________________________________________________________________
void G4SurfBits::Get(char *array) const
{
  // Copy all the byes.
  std::memcpy(array, fAllBits, (fNBits+7)>>3);
}

// If we are on a little endian machine, a bitvector represented using
// any integer type is identical to a bitvector represented using bytes.

//______________________________________________________________________________
void G4SurfBits::set(unsigned int nBits, const G4int *array)
{
  // set all the bytes.

  set(nBits, (const char*)array);
}

//______________________________________________________________________________
void G4SurfBits::Get(G4int *array) const
{
  // Get all the bytes.

  Get((char*)array);
}
