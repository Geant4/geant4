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
//

// class description:
//
//  This class represents a line of JPEG code. This class must be used exclusively
// by G4RTJpegCoder.
//

#ifndef G4RTOutBitStream_H
#define G4RTOutBitStream_H 1
#include "G4RTJpeg.hh"

static const u_char BitFullMaskT[8] = {0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f
, 0xff};

class G4OutBitStream
{
  public:
    G4OutBitStream(int size);
    ~G4OutBitStream();
        void SetBits(int v, int numBits);
        void SetByte(u_char dat);
        void SetWord(u_int dat);
        void CopyByte(const char* src, int n);

        u_char* GetStreamAddress(void){return mHeadOfBuf;};
        int GetStreamSize(void){return int(mBuf - mHeadOfBuf);};


  protected:
        u_char* mHeadOfBuf;
        u_char* mBuf;
        u_char* mEndOfBuf;
        int mBitPos;
        int mWriteFlag;

        void IncBuf(void);
        void FullBit(void);
        void Set8Bits(u_char v, int numBits);
        void SetFewBits(u_char v, int numBits);
        void SetBits2Byte(u_char v, int numBits);
};

#endif
