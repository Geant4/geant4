// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTJpegCoder.hh,v 1.3 2000-03-09 15:36:33 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4RTJpegCoder_H
#define G4RTJpegCoder_H 1
#include "G4RTJpeg.hh"

//HuffmanTable
struct
G4HuffmanCodeTable
{
        int     numOfElement;
        int*    SizeT;
        int*    CodeT;
};

class G4OutBitStream;

class G4JpegCoder
{
  public:
    G4JpegCoder(u_char* colorR,u_char* colorG,u_char* colorB);
    ~G4JpegCoder(void);

    void GetJpegData(char** aJpegData,int& size);

    void SetJpegProperty(const G4JpegProperty& aProperty);

    int DoCoding(void);


  protected:
        u_char* mRgb[3];
        int mYBlock[4][64];
        int mCbBlock[64];
        int mCrBlock[64];
        double mCosT[8][8];
        int mDCTData[64];
        int mPreDC[3];

        G4JpegProperty mProperty;
        int mNumVUnits;
        int mNumHUnits;

        G4OutBitStream *mOBSP;

        void CodeMCU();

        void makeYCC(int ux,int uy);

        void CodeHuffman(int cs);

        void ForwardDCT(int* picData);

        void Quantization(int cs);

        void WriteHeader(void);
        void WriteEOI(void);
};

const u_int JFIFLength = 16;
const u_int JFIFVersion = 0x0102;     //JFIF Ver 1.02
const u_char YSampleF = 0x22;
const u_char CSampleF = 0x11;

#endif
