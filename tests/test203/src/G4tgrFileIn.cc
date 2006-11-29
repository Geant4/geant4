#include "globals.hh"
#include "G4tgrFileIn.hh"
#include "G4tgrMessenger.hh"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <strstream>
////#include <strstream.h>

std::vector<G4tgrFileIn*> G4tgrFileIn::theInstances;

//-----------------------------------------------------------------------
G4tgrFileIn& G4tgrFileIn::GetInstance( const G4String& filename )
{
  std::vector<G4tgrFileIn*>::const_iterator vfcite;
  for( vfcite = theInstances.begin(); vfcite != theInstances.end(); vfcite++) {
    if( (*vfcite)->GetName() == filename) {
      return *(*vfcite);
    }
  }

  G4tgrFileIn* instance = 0;
  if( vfcite == theInstances.end() ) {
    instance = new G4tgrFileIn( filename );
    
    instance->theCurrentFile = -1;
    instance->OpenNewFile( filename.c_str() );

    theInstances.push_back( instance );
  }

  return *instance;
  
}

//-----------------------------------------------------------------------
void G4tgrFileIn::OpenNewFile( const char* filename )
{ 
  theCurrentFile++;
  std::ifstream* fin = new std::ifstream(filename);
  theFiles.push_back(fin);

  //-  G4int lineno = new G4int;
  //-  G4int lineno = 0;
  theLineNo.push_back( 0 );

  theNames.push_back( filename );

#ifndef OS_SUN_4_2
  if( !fin->is_open()) {
    G4cerr << "!!!! Input file does not exist: " << filename << G4endl;
    exit(1);
  }
#endif

}



//-----------------------------------------------------------------------
G4tgrFileIn& G4tgrFileIn::GetInstanceOpened( const G4String& filename )
{

  G4tgrFileIn& filein = G4tgrFileIn::GetInstance(filename);
  if (filein.GetName() != filename ) {
    G4cerr << "Error: file not opened yet " << filename << G4endl; 
    exit(0); 
  } else {
    return filein;
  }

}


//----------------------------------------------------------------------- 
G4int G4tgrFileIn::GetWordsInLine( std::vector<G4String>& wordlist)
{
  G4int isok = 1;

  //---------- Read a line of file:
  //@@@@--- Cannot be read with a istream_iterator, becasuse it uses G4cout, and then doesn't read '\n'
  //----- Clear wordlist
  G4int wsiz = wordlist.size();
  G4int ii;
  for (ii = 0; ii < wsiz; ii++) {
    wordlist.pop_back();
  } 

  //---------- Loop lines while there is an ending '\' or line is blank   
  const G4int NMAXLIN = 1000;
  char ltemp[NMAXLIN]; //there won't be lines longer than NMAXLIN characters
  for (;;) {
    (theLineNo[theCurrentFile])++;
    for( ii = 0; ii < NMAXLIN; ii++) ltemp[ii] = ' ';
    theFiles[theCurrentFile]->getline( ltemp, NMAXLIN ); 
    //---------- Check for lines longer than NMAXLIN character
    G4int ii;
    for ( ii=0; ii < NMAXLIN; ii++) {
      if ( ltemp[ii] == '\0' ) break;
    }
    if ( ii == NMAXLIN-1 ) {
      ErrorInLine();
      G4cerr << "!!!! line longer than " << NMAXLIN << " characters" << 
	G4endl << " please split it putting a '\\' at the end of line" << G4endl;
      exit(0);
    }
    
    //---------- End of file
    //-    if ( theFiles[theCurrentFile]->eof() ) {
    if ( EndOfFile() ) {
      //t          exit(0);
      return 0;
    }
    
    //---------- Convert line read to istrstream to split it in words 
    std::istrstream istr_line(ltemp);
     
    //--------- count how many words are there in ltemp (this sohuld not be needed, but sun compiler has problems) !! this has to be nvestigated...
    G4int NoWords = 0;
    char* tt = ltemp;
    G4String stemp(ltemp);
    do{ 
      if( *tt != ' ' && *(tt) != '\0' ) {
	if( tt == ltemp) {
	  NoWords++;
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) G4cout << "G4tgrFileIn NoWords" << NoWords << ltemp << G4endl;
#endif
	} else if( *(tt-1) == ' ' ||  *(tt-1) == '\015' ||  *(tt-1) == '\t') {
	  NoWords++; 
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	      G4cout << "G4tgrFileIn NoWords" << NoWords << ltemp << G4endl;
#endif
	}
      }
      tt++;
    }while(*tt != '\0' & stemp.length()!=0);
    G4String stempt (ltemp);
    if(stempt.length() == 0) NoWords = 0;
    
    //--------- Read words from istr_line and write them into wordlist
    for( ii=0; ii < NoWords; ii++) {
      G4String stemp = "";
      istr_line >> stemp;   //?? gives warning in Insure++
      if ( stemp.length() == 0 ) break;
      G4int comment = stemp.find(G4String("//") );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	G4cout << "!!!COMMENT" << comment << stemp.c_str() << G4endl;
#endif
      if ( comment == 0 ) {
	break; 
      } else if ( comment > 0 ) {
	stemp = stemp.substr( 0, comment );
	wordlist.push_back(stemp);
	break;
      } 
      wordlist.push_back(stemp);
    }
    
    //These two algorithms should be the more STL-like way, but they don't work for files whose lines end without '\015'=TAB (STL problem: doesn't find end of string??)
    // istream_iterator<G4String, ptrdiff_t> G4String_iter(istr_line);
    // istream_iterator<G4String, ptrdiff_t> eosl;
    // copy(G4String_iter, eosl, back_inserter(wordlist));
    // typedef istream_iterator<G4String, ptrdiff_t> G4String_iter;
    // copy(G4String_iter(istr_line), G4String_iter(), back_inserter(wordlist));
    
    if ( wordlist.size() != 0 ) {
      if( (*(wordlist.end()-1)).compare("\\") == 0 ) {   //use '\' to mark continuing line  
	wordlist.pop_back();
      } else {
	break;
      }
    }
  }
  
  //--------- A pair of double quotes delimits a word, therefore, look for the case where there is more than one word between two double quotes
  std::vector<G4String> wordlist2;
  G4String wordq;
  uint imerge = 0;
  for( ii = 0; ii < wordlist.size(); ii++) {
    if( wordlist[ii].substr(0,1) == "\"" ) {
      //-      cout << " l quote found " << wordlist[ii] << endl;
      imerge = 1;
    } 
    if( wordlist[ii][ wordlist[ii].size()-1 ] == '\"' ) {
      if( imerge != 1 ) {
	DumpException(" word with trailing '\"' while there is no previous word with leading '\"' in line " );
      }
      //-      cout << " t quote found " << wordlist[ii] << endl;
      imerge = 2;
    }
    if( imerge == 0 ) {
      wordlist2.push_back( wordlist[ii] );
    } else if( imerge == 1 ) {
      wordq.append( wordlist[ii].substr(1,wordlist[ii].size()) );
      wordq.append(" ");
    } else if( imerge == 2 ) {
      wordq.append( wordlist[ii].substr(0,wordlist[ii].size()-1));
      wordlist2.push_back( wordq );
      wordq = "";
      imerge = 0;
    }
  }
  if( imerge == 1 ) {
    DumpException(" word with leading '\"' in line while there is no later word with trailing '\"' in line " );
  }

  wordlist = wordlist2;

  //or why not like this?:
  //typedef istream_iterator<G4String, ptrdiff_t> string_iter;
  //copy(string_iter(istr_line), string_iter(), back_inserter(wordlist));
  
  // check if including a new file
  if( wordlist[0] == "#include" ) {
    if( wordlist.size() != 2 ) {
      ErrorInLine();
      std::cerr << "'#include' should have as second argument the filename " << G4endl;
      exit(0);
    }

#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4cout << "include found " << G4endl;
#endif
    OpenNewFile( wordlist[1].c_str() );
    isok = GetWordsInLine( wordlist);

  }

  return isok;  

}


//-----------------------------------------------------------------------
void G4tgrFileIn::ErrorInLine()
{
  std::cerr << "!! EXITING: ERROR IN LINE No " << theLineNo[theCurrentFile] << " file: " << theNames[theCurrentFile] << " : ";

}


//-----------------------------------------------------------------------
G4bool G4tgrFileIn::EndOfFile()
{
  G4bool isok = theFiles[theCurrentFile]->eof();
  if( isok ) {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4cout << " eof theCurrentFile " << theCurrentFile << G4endl;
#endif
    theCurrentFile--;
    if( theCurrentFile != -1 ) Close();  // last file will be closed by the user
  }
  //only real closing if all files are closed
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " eof " << isok << " " << theCurrentFile << G4endl;
#endif
  if( theCurrentFile != -1 ) { 
    return 0;
  } else {
    return isok;
  }
}


//-----------------------------------------------------------------------
void G4tgrFileIn::Close()
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "G4tgrFileIn close " << theCurrentFile << " size " << theFiles.size() << G4endl;
#endif
  //  if( theCurrentFile+1 != 0 ) {
  //  ErrorInLine();
  //  std::cerr << "trying to close file while reading other files included in it " << theCurrentFile+1 << G4endl;
  //  //    exit(0);
  //  } else { 
    theFiles[theCurrentFile+1]->close();
    theFiles.pop_back();
    //  }
}

//-----------------------------------------------------------------------
void G4tgrFileIn::DumpException( const G4String& sent )
{
  std::cerr << "!!! EXITING: " << sent << " in file " << theName << " line No " << theLineNo[theCurrentFile] << G4endl;
  exit(1); 

}

