// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTJpeg.hh,v 1.4 2000-03-09 17:38:31 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

// class description:
//
//  This header file defines some static constant variables and error classes
// used internally by G4JpegMaker and related classes
//

#ifndef G4RTJpeg_H
#define G4RTJpeg_H 1

typedef	unsigned char	u_char;
typedef unsigned int	u_int;

const char      JFIF[] = "JFIF";
const char      JFXX[] = "JFXX";

const double    Sqrt2 = 1.41421356;
const double    DisSqrt2 = 1.0 / Sqrt2;
const double    PaiDiv16 = 3.14159265 / 16;

//Zigzag
static const int Zigzag[64] = {
                 0,  1,  8, 16,  9,  2,  3, 10,
                17, 24, 32, 25, 18, 11,  4,  5,
                12, 19, 26, 33, 40, 48, 41, 34,
                27, 20, 13,  6,  7, 14, 21, 28,
                35, 42, 49, 56, 57, 50, 43, 36,
                29, 22, 15, 23, 30, 37, 44, 51,
                58, 59, 52, 45, 38, 31, 39, 46,
                53, 60, 61, 54, 47, 55, 62, 63
};

//ProcessResult
enum
jProcessResult{
  M_NoError = 0,
  M_RuntimeError = -1,
  M_DataError = -2
};

// JpegMarkerCode
enum
jMarker{

        M_SOF0  = 0xc0,
        M_SOF1  = 0xc1,
        M_SOF2  = 0xc2,
        M_SOF3  = 0xc3,

        M_SOF5  = 0xc5,
        M_SOF6  = 0xc6,
        M_SOF7  = 0xc7,

        M_JPG   = 0xc8,
        M_SOF9  = 0xc9,
        M_SOF10 = 0xca,
        M_SOF11 = 0xcb,

        M_SOF13 = 0xcd,
        M_SOF14 = 0xce,
        M_SOF15 = 0xcf,

        M_DHT   = 0xc4,

        M_DAC   = 0xcc,

        M_RST0  = 0xd0,         M_RST1  = 0xd1,
        M_RST2  = 0xd2,         M_RST3  = 0xd3,
        M_RST4  = 0xd4,         M_RST5  = 0xd5,
        M_RST6  = 0xd6,         M_RST7  = 0xd7,

        M_SOI   = 0xd8,
        M_EOI   = 0xd9,
        M_SOS   = 0xda,
        M_DQT   = 0xdb,
        M_DNL   = 0xdc,
        M_DRI   = 0xdd,
        M_DHP   = 0xde,
        M_EXP   = 0xdf,
        M_COM   = 0xfe,

        M_APP0  = 0xe0,         M_APP1  = 0xe1,
        M_APP2  = 0xe2,         M_APP3  = 0xe3,
        M_APP4  = 0xe4,         M_APP5  = 0xe5,
        M_APP6  = 0xe6,         M_APP7  = 0xe7,
        M_APP8  = 0xe8,         M_APP9  = 0xe9,
        M_APP10 = 0xea,         M_APP11 = 0xeb,
        M_APP12 = 0xec,         M_APP13 = 0xed,
        M_APP14 = 0xee,         M_APP15 = 0xef,


        M_JPG0  = 0xf0,         M_JPG1  = 0xf1,
        M_JPG2  = 0xf2,         M_JPG3  = 0xf3,
        M_JPG4  = 0xf4,         M_JPG5  = 0xf5,
        M_JPG6  = 0xf6,         M_JPG7  = 0xf7,
        M_JPG8  = 0xf8,         M_JPG9  = 0xf9,
        M_JPG10 = 0xfa,         M_JPG11 = 0xfb,
        M_JPG12 = 0xfc,         M_JPG13 = 0xfd,


        M_TEM   = 0x01,
        M_RESst  = 0x02,
        M_RESend = 0xbf,

        M_Error  = 0xff,
        M_Marker  = 0xff
};

//JpegProperty
struct
G4JpegProperty{
  int nRow;
  int nColumn;
  int Dimension;
  int SamplePrecision;
  const char * Comment;
  int Format;
  u_char MajorRevisions;
  u_char MinorRevisions;
  int Units;
  int HDensity;
  int VDensity;
  int HThumbnail;
  int VThumbnail;
  int ExtensionCode;
};


//MemoryError
class G4MemoryError
{
  public:
    G4MemoryError(int size, const char* message)
        {mSize = size;          mMessage = message;};
        int mSize;
        const char* mMessage;
};

//IndexError
class G4IndexError
{
  public:
    G4IndexError(int maxIndex, int errorIndex, const char* mes)
        {mMaxIndex = maxIndex;  mErrorIndex = errorIndex;  mMessage = mes;};
        int mMaxIndex;
        int mErrorIndex;
        const char* mMessage;
};

//BufferError
class G4BufferError
{
  public:
    G4BufferError(const char* mes)
        {mMessage = mes;};
        const char* mMessage;
};

//DataFormatError
class G4DataFormatError
{
  public:
    G4DataFormatError(void* address, const char* message)
        {mAddress = address;            mMessage = message;};
        void* mAddress;
        const char* mMessage;
};


//NotSupported
class G4NotSupported
{
  public:
    G4NotSupported(jMarker aMark, void* address)
        {mMark = aMark;         mAddress = address;};
        jMarker mMark;
        void* mAddress;
};

#endif
