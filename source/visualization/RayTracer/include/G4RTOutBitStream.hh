// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTOutBitStream.hh,v 1.5 2000-03-09 17:38:32 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
        void SetBits(int v, int numBits);
        void SetByte(u_char dat);
        void SetWord(u_int dat);
        void CopyByte(const char* src, int n);

        u_char* GetStreamAddress(void){return mHeadOfBuf;};
        int GetStreamSize(void){return mBuf - mHeadOfBuf;};


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
