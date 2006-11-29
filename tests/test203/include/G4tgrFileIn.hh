#ifndef G4tgrFileIn_HH
#define G4tgrFileIn_HH

#include "globals.hh"
#include <vector>


class G4tgrFileIn 
{
public:
  G4tgrFileIn(){};
  ~G4tgrFileIn(){}
private:
  G4tgrFileIn( const G4String& name ): theName(name){}
  
 public:
  // Get the only instance opening the file
  static G4tgrFileIn& GetInstance( const G4String& name ); 
  
  // Get the only instance when file should be already opened
  static G4tgrFileIn& GetInstanceOpened( const G4String& name ); 
  
  // Read a line and transform it to a vector of words 
  G4int GetWordsInLine( std::vector<G4String>& wl );

  // Print out an error message indicating the line being read
  void ErrorInLine();
  
  // Access data members
  const G4int Nline() { return theLineNo[theCurrentFile]; }
  
  const G4String& GetName() { return theName; }
  
  G4bool EndOfFile();
  void Close();
  void DumpException( const G4String& sent );
  
  // private:
  void OpenNewFile( const char* filename );
  
  // private:
  std::vector< std::ifstream* > theFiles;
  // Number of line being read
  std::vector<G4int> theLineNo;
  std::vector<G4String> theNames;
  int theCurrentFile;  // index of file being read in theFiles
  
  // private DATA MEMEBERS
  // Vector of class instances (each one identified by its name)
  static std::vector<G4tgrFileIn*> theInstances;
  
  // Name of file
  G4String theName;
  
};

#endif 
