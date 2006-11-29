//Id:  G4tgrFileOut.h
//   ostream class for handling the output
// 
//   History: v1.0 
//   Pedro Arce 3/98

#ifndef FILEOUT_H
#define FILEOUT_H

#include <fstream>
#include <iostream>
#include <vector>


class G4tgrFileOut : public ofstream 
{
public:
  G4tgrFileOut(){};
  G4tgrFileOut( const G4String& name ): ofstream(), theName(name){};
  ~G4tgrFileOut(){};

  // get the instance of file with name filename
  static G4tgrFileOut& GetInstance( const G4String& filename );

 // Access data members
  const G4String GetName() { return theName; }

// private DATA MEMEBERS
private:
  // Class only instance
  static vector<G4tgrFileOut*> theInstances;

  /// Name of file
  G4String theName; 
};

#endif 

