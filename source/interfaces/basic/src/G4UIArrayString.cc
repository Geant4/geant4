// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIArrayString.cc,v 1.2 2000-06-14 03:18:59 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/iomanip"
#include "G4UIArrayString.hh"

static const char strESC= '\033';

////////////////////////////////////////////////////////
G4UIArrayString::G4UIArrayString(const G4String& stream)
////////////////////////////////////////////////////////
{
  nElement=0;
  nColumn=5;  // temporal assignment

  G4String tmpstr= stream;  // G4String::strip() CONST !!
  G4String astream= tmpstr.strip(G4String::both);

  // tokenize...
  G4int indx=0;
  while(1) {
    G4int jc= astream.index(" ", indx);
    nElement++;
    if(jc == G4String::npos) break;
    for(G4int i=1; jc+i< astream.length(); i++ ) {  // skip continuing spaces
      if(astream[jc+i]==' ') jc++;
      else break;
    }
    indx= jc+1;
  }

  // allocate string array
  stringArray= new G4String[nElement];   

  // push...
  indx=0;
  for(G4int i=0; i<nElement; i++){
    G4int jc= astream.index(" ", indx);
    if(jc != G4String::npos)
      stringArray[i]= astream(indx, jc-indx);
    else {  // last token
      jc= astream.length()+1;
      stringArray[i]= astream(indx, jc-indx);
    }
    for(G4int j=1; jc+j< astream.length(); j++ ) { // skip continuing spaces
      if(astream(jc+j)==' ') jc++;
      else break;
    }
    indx= jc+1;
  }
}

///////////////////////////////////
G4UIArrayString::~G4UIArrayString()
///////////////////////////////////
{ 
  delete [] stringArray;
}

///////////////////////////////////////////////////////////////////
G4String* G4UIArrayString::GetElement(G4int icol, G4int irow) const
///////////////////////////////////////////////////////////////////
{  
  if( !(icol>=1 && irow>=1)) // offset of column/row is "1".
    G4cerr << "G4UIArrayString: overrange" << G4endl;
  if(icol>nColumn) G4cerr << "G4UIArrayString: overrange" << G4endl;

  G4int jq= (irow-1)*nColumn + icol;
  if(jq> nElement) G4cerr << "G4UIArrayString: overrange" << G4endl;

  jq--;
  return &stringArray[jq];
}

////////////////////////////////////////////
G4int G4UIArrayString::GetNRow(int icol) const
////////////////////////////////////////////
{
  G4int ni;
  if(nElement%nColumn ==0) ni= nElement/nColumn;
  else ni= nElement/nColumn + 1;

  G4int nn= nElement%nColumn;
  if(nn==0) nn= nColumn;

  if(icol<= nn) return ni;
  else return ni-1;
}

////////////////////////////////////////////////
G4int G4UIArrayString::GetNField(int icol) const
////////////////////////////////////////////////
{
  G4int maxWidth=0;
  for (G4int iy=1; iy<= GetNRow(icol); iy++) {
    G4int ilen= GetElement(icol,iy)-> length();
    // care for color code
    // if(GetElement(icol,iy)-> index(strESC,0) != G4String::npos) {
    // if(strESC == (*GetElement(icol,iy))[0] ) {
    const char tgt = (*GetElement(icol,iy))[0];
    if(strESC == tgt) {
      ilen-= 5;
      if(ilen<0) G4cout << "length(c) cal. error." << G4endl;
    }
    if(ilen> maxWidth) maxWidth= ilen;
  }

  return maxWidth;
}

/////////////////////////////////////////////////
int G4UIArrayString::CalculateColumnWidth() const
/////////////////////////////////////////////////
{
  G4int totalWidth= 0;

  for(G4int ix=1; ix<= nColumn; ix++) {
    totalWidth+= GetNField(ix);
  }

  const G4int nwSpace= 2;
  totalWidth+= (nColumn-1)*nwSpace;  // for space

  return totalWidth;
}

//////////////////////////////////////
void G4UIArrayString::Show(G4int ncol)
//////////////////////////////////////
{
  // calculate #colums in need...
  while( CalculateColumnWidth()< ncol ) {
    nColumn++;
  }
  while( CalculateColumnWidth()> ncol && nColumn>1 ) {
    nColumn--;
  }
  
  for(G4int iy=1; iy<= GetNRow(1); iy++) {
    G4int nc= nColumn;
    if(iy == GetNRow(1)) { // last row
      nc= nElement%nColumn;
      if(nc==0) nc= nColumn;
    }
    for(G4int ix=1; ix<=nc; ix++) {
      G4String word= GetElement(ix,iy)-> data();

      // care for color code
      G4String colorWord;
      //if(word.index(strESC,0) != G4String::npos) {
      //if(strESC == word[0]) {
      const char tgt = word[0];
      if(strESC == tgt) {
        colorWord= word(0,5);
        word.erase(0,5);
      }
      if(!colorWord.empty()) G4cout << colorWord << G4std::flush;

      G4cout << G4std::setiosflags(G4std::ios::left) << G4std::setw(GetNField(ix)) 
             << word << G4std::flush;
      if(ix != nc) G4cout << "  " << G4std::flush;
      else G4cout << G4endl;      
    }
  }
}

