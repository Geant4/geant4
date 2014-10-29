#ifndef G4LEPTSDiffXS_h
#define G4LEPTSDiffXS_h 1


#include <string>
using namespace std;

class G4LEPTSDiffXS {

public:

  G4LEPTSDiffXS( string);   // Constructor

  void readDXS();    // Read file
  void BuildCDXS();
  void BuildCDXS(G4double, G4double);
  void NormalizeCDXS();
  void InterpolateCDXS();
  void PrintDXS(int);

  G4double SampleAngle(G4double);
  G4double SampleAngleMT(G4double, G4double);
  G4double SampleAngleEthylene(G4double, G4double);
  G4bool IsFileFound() const {
    return bFileFound;
  }

private:
  string fileName;
  int NumAng;
  int INumAng;
  int NumEn;
  char DXSTypeName[8];
  int DXSType;
  G4double Eb[100];
  //  G4double DXS[100][190], CDXS[100][190], IDXS[100][19000], ICDXS[100][19000];
  G4double DXS[100][190], CDXS[100][190], ICDXS[100][19000];
  //  G4double KT[100][190],  CKT[100][190],  IKT[100][19000];
  G4double KT[100][190],  IKT[100][19000];

  G4bool bFileFound;
};

#endif
