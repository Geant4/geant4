// $Id: G4RootTransManager.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootTransManager.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef ROOT_TRANS_MANAGER_HH
#define ROOT_TRANS_MANAGER_HH 1

#include <string>
#include <map>

class G4RootIOManager;
class TFile;
class TTree;

typedef std::map<std::string,TFile*,std::less<std::string> > RootFileCatalog;

// Class inherited:
#include "G4VTransactionManager.hh"

// Class Description:
//   Implementation of transaction manager with HepODBMS persistency package

class G4RootTransManager
 : public G4VTransactionManager
{
    public: // With description
      G4RootTransManager();
      // Constructor

      virtual ~G4RootTransManager() {};
      // Destructor

    public: // With description
      bool SelectReadFile(std::string obj, std::string file);
      // Select the input file for type "obj" and sets the container ref.

      bool SelectWriteFile(std::string obj, std::string file);
      // Select the output file for type "obj" and sets the container ref.

      bool StartUpdate() {return true;};
      // Start the database update transaction.  Dummy in ROOT I/O.

      bool StartRead() {return true;};
      // Start the database read transaction.  Dummy in ROOT I/O.

      void Commit();
      // Commit the changes to the ROOT I/O files.

      void Abort() {};
      // Abort the transaction.  Dummy in ROOT I/O.

      TFile* LastReadFile() { return f_lastreadfile; };
      // returns the pointer of the last read file

      TFile* LastWriteFile() { return f_lastwritefile; };
      // returns the pointer of the last read file

      TTree* LastReadTree() { return f_lastreadtree; };
      // returns the pointer of the last read tree

      TTree* LastWriteTree() { return f_lastwritetree; };
      // returns the pointer of the last read tree

      void Initialize();
      // Initialize the ROOT I/O file map.

      int BufferSize();
      // Returns the ROOT I/O buffer size.

      int SplitLevel() { return f_split; };
      // Returns the ROOT I/O split level.

      int CompressionLevel() { return f_comp; };
      // Returns the ROOT I/O compression level.

      void SetBufferSize(int buf) { f_bufsize = buf; };
      // Sets the ROOT I/O buffer size.

      void SetSplitLevel(int s) { f_split = s; };
      // Sets the ROOT I/O split level.

      void SetCompressionLevel(int c) { f_comp = c; };
      // Sets the ROOT I/O compression level.

    private:
      TFile* LookUp(string file, string mode);
      // Returns the Root file pointer to the object type and file name.
      // Opens a new file if the access is for the first time.

    private:
      RootFileCatalog f_catalog;
      TFile*          f_lastreadfile;
      TFile*          f_lastwritefile;
      TTree*          f_lastreadtree;
      TTree*          f_lastwritetree;
      std::string     f_treeName;
      std::string     f_treeDesc;
      int f_bufsize;
      int f_split;
      int f_comp;

}; // End of class G4RootTransManager

#endif

