#include <iomanip.h>
#include <fstream.h>
#include <strstream.h>

#define maxLineLength 80

template<class t>
ostream& operator<<(ostream& o,ReadList<t>& L) 
{
  for (vector<t>::iterator X=L.begin(); X != L.end(); X++) {
    o << *X << endl;
  }
  return o;
}


template<class t>
ReadList<t>::ReadList(char* fileName) 
  : vector<t>(),N(0),file(new ifstream(fileName) ) {}

template<class t>
void ReadList<t>::readIn()
{
  char Line[120];
  int noLines = 0;
  char c;
  while ( !file->eof() && file->get(Line,maxLineLength,'\n') ) {
    file->get(c);
    ++noLines;
    int j=0;
    while ( Line[j] == ' ' ) j++;
    if ( Line[j] != '#' && Line[j] ) {
      istrstream inputLine(Line,strlen(Line));
      try {
	t h(inputLine);
	bookIn(h);
	++N;
      }
      catch (char*) { 
	cerr << "Read Error in Line " << noLines << ":\n";
	throw;
      }
    }
  }
}

template<class t>
ReadList<t>::~ReadList()
{
  delete file;
  erase(begin(),end());
}

template<class t>
int ReadList<t>::getIndex(const t& h)
{
  for (int i=0; i<N; i++) 
    if ( (*this)[i] == h ) 
      return i;
  return -1;
};

// ----------------------------------------------------------------------

template<class t>
ostream& operator<<(ostream& o,ReadList_P<t>& L) 
{
  for (vector<t>::iterator X=L.begin(); X != L.end(); X++) {
    o << *(*X) << endl;
  }
  return o;
}


template<class t>
ReadList_P<t>::ReadList_P(char* fileName) 
  : vector<t*>(),N(0),file(new ifstream(fileName) ) {}

template<class t>
void ReadList_P<t>::readIn()
{
  char Line[120];
  int noLines = 0;
  char c;
  while ( !file->eof() && file->get(Line,maxLineLength,'\n') ) {
    file->get(c);
    ++noLines;
    int j=0;
    while ( Line[j] == ' ' ) j++;
    if ( Line[j] != '#' && Line[j] ) {
      istrstream inputLine(Line,strlen(Line));
      try {
	t* h = new t(inputLine);
	bookIn(*h);
	++N;
      }
      catch (char*) { 
	cerr << "Read Error in Line " << noLines << ":\n";
	throw;
      }
    }
  }
}

template<class t>
ReadList_P<t>::~ReadList_P()
{
  delete file;
  for (int i=0; i<size(); i++)
    delete (*this)[i];
  erase(begin(),end());
}

template<class t>
int ReadList_P<t>::getIndex(const t& h)
{
  for (int i=0; i<N; i++) 
    if ( *(*this)[i] == h ) 
      return i;
  return -1;
};


