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
//

#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "G4RTJpeg.hh"
#include "G4RTOutBitStream.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTJpegCoder.hh"
#include "G4RTJpegCoderTables.hh"


G4JpegCoder::G4JpegCoder(u_char* colorR,u_char* colorG,u_char* colorB)
{
  mRgb[0] = colorR;
  mRgb[1] = colorG;
  mRgb[2] = colorB;

  mPreDC[0] = mPreDC[1] = mPreDC[2] = 0;
  mOBSP = 0;

  for(int n=0; n<8; n++)
                for(int im=0; im<8; im++)
                                mCosT[n][im] = std::cos((2 * im + 1) * n * PaiDiv16);
}

G4JpegCoder::~G4JpegCoder(void)
{}

void
G4JpegCoder::GetJpegData(char** aJpegData, int& size)
{
  if (mOBSP != 0){
    *aJpegData = (char*)mOBSP->GetStreamAddress();
        size = mOBSP->GetStreamSize();
        }
        else{
          *aJpegData = 0;
          size = 0;
        }

}

int
G4JpegCoder::DoCoding(void)
{
  mNumVUnits = (mProperty.nRow / 16) + ((mProperty.nRow % 16) ? 1 : 0);
  mNumHUnits = (mProperty.nColumn / 16) + ((mProperty.nColumn % 16) ? 1 : 0);

  int size = mProperty.nColumn * mProperty.nRow * 3;
  if(size < 10240)
    size = 10240;

  try{
       mOBSP = new G4OutBitStream(size);
           WriteHeader();
           for(int yu=0; yu<mNumVUnits; yu++){
                   for(int xu=0; xu<mNumHUnits; xu++){
                           makeYCC(xu, yu);

  //mRgb->YCrCb
  #ifdef GRAY
  for(int i=0; i<64; i++)
        mCbBlock[i] = mCrBlock[i] = 0;
  #endif
         CodeMCU();
         }
       }
       WriteEOI();
           return M_NoError;
     }

  catch(G4MemoryError &me){
                return M_RuntimeError;
  }
  catch(G4BufferError &be){
                return M_RuntimeError;
  }
  catch(G4IndexError &ie){
                return M_RuntimeError;
  }
}

//MCU
void
G4JpegCoder::CodeMCU(void)
{
  for(int n=0; n<4; n++){
                ForwardDCT(mYBlock[n]);
                Quantization(0);
                CodeHuffman(0);
  }
  ForwardDCT(mCbBlock);
  Quantization(1);
  CodeHuffman(1);

  ForwardDCT(mCrBlock);
  Quantization(2);
  CodeHuffman(2);
}

void
G4JpegCoder::makeYCC(int ux, int uy)
{
  u_char  rv, gv, bv;
  int tCrBlock[4][64];
  int tCbBlock[4][64];

  for(int u=0; u<4; u++){
        int *yp = mYBlock[u];
        int *cbp = tCbBlock[u];
        int *crp = tCrBlock[u];

        int sx = ux * 16 + ((u&1) ? 8 : 0);
        int ex = sx + 8;
        int sy = uy * 16 + ((u>1) ? 8 : 0);
        int ey = sy + 8;

     for(int iv=sy; iv<ey; iv++){
          int ii = iv < mProperty.nRow ? iv : mProperty.nRow - 1;
          for(int ih=sx; ih<ex; ih++){
            int jj = ih < mProperty.nColumn ? ih : mProperty.nColumn - 1;
                int index = ii * mProperty.nColumn + jj;
                rv = mRgb[0][index];
                gv = mRgb[1][index];
                bv = mRgb[2][index];

                *yp++ = int((0.2990 * rv) + (0.5870 * gv) + (0.1140 * bv) - 128)
;
                *cbp++ = int(-(0.1687 * rv) - (0.3313 * gv) + (0.5000 * bv));
                *crp++ = int((0.5000 * rv) - (0.4187 * gv) - (0.0813 * bv));
                                }       // ih
                           }    //iv
  }     //u

 int   n = 0;
  for(int b=0; b<4; b++){
        switch(b){
                case 0:         n=0;    break;
                case 1:         n=4;    break;
                case 2:         n=32;   break;
                case 3:         n=36;
        }
        for(int y=0; y<8; y+=2){
                for(int x=0; x<8; x+=2){
                        int idx = y * 8 + x;
                        mCrBlock[n] = tCrBlock[b][idx];
                        mCbBlock[n] = tCbBlock[b][idx];
                        n++;
                }
                n += 4;
        }
  }
}

void
G4JpegCoder::CodeHuffman(int cs)
{
  const G4HuffmanCodeTable& dcT = cs ? CDcHuffmanT : YDcHuffmanT;
  const G4HuffmanCodeTable& acT = cs ? CAcHuffmanT : YAcHuffmanT;
  const int eobIdx = cs ? CEOBidx : YEOBidx;
  const int zrlIdx = cs ? CZRLidx : YZRLidx;

  int diff = mDCTData[0] - mPreDC[cs];
  mPreDC[cs] = mDCTData[0];
  int absDiff = std::abs(diff);
  int dIdx = 0;

  while(absDiff > 0){
        absDiff >>= 1;
        dIdx++;
  }
  if(dIdx > dcT.numOfElement)
        throw(G4IndexError(dcT.numOfElement, dIdx, "CodeHuffman:DC"));
  mOBSP->SetBits((dcT.CodeT)[dIdx], (dcT.SizeT)[dIdx]);

  if(dIdx){
        if(diff < 0)
                diff--;
        mOBSP->SetBits(diff, dIdx);
  }

  int run = 0;
  for(int n=1; n<64; n++){
        int absCoefficient = std::abs( mDCTData[ Zigzag[n] ] );
        if( absCoefficient ){
                while( run > 15 ){
                mOBSP->SetBits((acT.CodeT)[zrlIdx], (acT.SizeT)[zrlIdx]);
                run -= 16;
                }
                int is = 0;
                while( absCoefficient > 0 ){
                        absCoefficient >>= 1;
                        is++;
                }
                int     aIdx = run * 10 + is + (run == 15);
                if( aIdx >= acT.numOfElement ) 
                  throw( G4IndexError( acT.numOfElement, aIdx, "CodeHuffman:AC" )
 );
  mOBSP->SetBits( (acT.CodeT)[aIdx], (acT.SizeT)[aIdx] );
                int     v = mDCTData[ Zigzag[n] ];
                if( v < 0 )
                  v--;
                mOBSP->SetBits( v, is );
                run = 0;
                }
                else{
                  if(n == 63)
                  mOBSP->SetBits( (acT.CodeT)[eobIdx], (acT.SizeT)[eobIdx] );
                else
                run++;
                }
        }
}


void
G4JpegCoder::Quantization(int cs)
{
  int* qt = (int*)(cs ? CQuantumT : YQuantumT);
  for( int i=0; i<64; i++ ){
    mDCTData[i] /= qt[i];
  }
}


void
G4JpegCoder::ForwardDCT(int* picData)
{
  for( int v=0; v<8; v++ ){
    double cv = v ? 1.0 : DisSqrt2;
    for( int u=0; u<8; u++ ){
      double cu = u ? 1.0 : DisSqrt2;
      double sum = 0;
      for( int y=0; y<8; y++ )
        for( int x=0; x<8; x++ )
          sum += picData[ y * 8 + x ] * mCosT[u][x] * mCosT[v][y];
      mDCTData[ v * 8 + u ] = int( sum * cu * cv / 4 );
    }
  }
}


void
G4JpegCoder::WriteHeader( void )
{
  int i = 0;    //counter
  //SOI
  mOBSP->SetByte( M_Marker );   //FF
  mOBSP->SetByte( M_SOI );      //SOI

  //APP0(JFIF Header)
  mOBSP->SetByte( M_Marker );   //FF
  mOBSP->SetByte( M_APP0 );     //APP0
  mOBSP->SetWord( JFIFLength );        //parameter
  mOBSP->CopyByte( (char*)JFIF, 5 );   //"JFIF\0"
  mOBSP->SetWord( JFIFVersion );       //Version
  mOBSP->SetByte( mProperty.Units );
  mOBSP->SetWord( mProperty.HDensity );
  mOBSP->SetWord( mProperty.VDensity );
  mOBSP->SetByte( 0 );
  mOBSP->SetByte( 0 );

 //comment
  if( mProperty.Comment != 0 ){
    mOBSP->SetByte( M_Marker ); //FF
    mOBSP->SetByte( M_COM );    //comment
    int length = (int)strlen( mProperty.Comment ) + 1;
    mOBSP->SetWord( length + 2 );
    mOBSP->CopyByte( mProperty.Comment, length );
  }

  //DQT
  mOBSP->SetByte( M_Marker );
  mOBSP->SetByte( M_DQT );
  mOBSP->SetWord( 67 );
  mOBSP->SetByte( 0 );
  for( i=0; i<64; i++ )
        mOBSP->SetByte( u_char( YQuantumT[Zigzag[i]] ) );
  mOBSP->SetByte( M_Marker );
  mOBSP->SetByte( M_DQT );
  mOBSP->SetWord( 67 );
  mOBSP->SetByte( 1 );
  for( i=0; i<64; i++ )
        mOBSP->SetByte( u_char( CQuantumT[Zigzag[i]] ) );
   // DHT
  mOBSP->CopyByte( (char*)YDcDht, DcDhtLength );
  mOBSP->CopyByte( (char*)CDcDht, DcDhtLength );
  mOBSP->CopyByte( (char*)YAcDht, AcDhtLength );
  mOBSP->CopyByte( (char*)CAcDht, AcDhtLength );

  // Frame Header
  mOBSP->SetByte( M_Marker );   // FF
  mOBSP->SetByte( M_SOF0 );
  mOBSP->SetWord( 3 * mProperty.Dimension + 8 );
  mOBSP->SetByte( mProperty.SamplePrecision );
  mOBSP->SetWord( mProperty.nRow );
  mOBSP->SetWord( mProperty.nColumn );
  mOBSP->SetByte( mProperty.Dimension );

  mOBSP->SetByte( 0 );
  mOBSP->SetByte( YSampleF );
  mOBSP->SetByte( 0 );

  mOBSP->SetByte( 1 );
  mOBSP->SetByte( CSampleF );

  mOBSP->SetByte( 1 );
  mOBSP->SetByte( 2 );
  mOBSP->SetByte( CSampleF );
  mOBSP->SetByte( 1 );

  //Scan Header
  mOBSP->SetByte( M_Marker );
  mOBSP->SetByte( M_SOS );
  mOBSP->SetWord( 2 * mProperty.Dimension + 6 );
  mOBSP->SetByte( mProperty.Dimension );
  for( i=0; i<mProperty.Dimension; i++ ){
        mOBSP->SetByte( i );
        mOBSP->SetByte( i==0 ? 0 : 0x11 );
  }
  mOBSP->SetByte( 0 );  //Ss
  mOBSP->SetByte( 63 ); //Se
  mOBSP->SetByte( 0 );  //Ah,Al
}

//EOI
void
G4JpegCoder::WriteEOI( void )
{
  mOBSP->SetByte( M_Marker );
  mOBSP->SetByte( M_EOI );
}

//SetJpegProperty
void
G4JpegCoder::SetJpegProperty(const G4JpegProperty& aProperty )
{
  mProperty = aProperty;
  mProperty.Dimension = 3;
  mProperty.SamplePrecision = 8;
  mProperty.Format = 1;
  mProperty.MajorRevisions = 1;
  mProperty.MinorRevisions = 2;
  mProperty.HThumbnail = 0;
  mProperty.VThumbnail = 0;
}
