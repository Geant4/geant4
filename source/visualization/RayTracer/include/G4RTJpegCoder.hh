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
// $Id: G4RTJpegCoder.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
//

// class description:
//
//  This class converts 8 bit unsighned ints array to JPEG code.
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
