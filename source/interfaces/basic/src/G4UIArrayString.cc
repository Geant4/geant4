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

#include <iomanip>
#include "G4UIArrayString.hh"

static const char strESC= '\033';

////////////////////////////////////////////////////////
G4UIArrayString::G4UIArrayString(const G4String& stream)
////////////////////////////////////////////////////////
{
  nElement=0;
  nColumn=5;  // temporal assignment

  G4String astream = G4StrUtil::strip_copy(stream);

  // tokenize...
  std::size_t indx=0;
  while(1) {
    std::size_t jc= astream.find(" ", indx);
    nElement++;
    if(jc == G4String::npos) break;
    jc++; // fix a tiny mistake...
    for(; jc< astream.length(); ) {  // skip continuing spaces
      if(astream[(G4int)jc]==' ') jc++;
      else break;
    }
    indx= jc;
  }

  // allocate string array
  stringArray= new G4String[nElement];   

  // push...
  indx=0;
  for(G4int i=0; i<nElement; ++i){
    std::size_t jc= astream.find(" ", indx);
    if(jc != G4String::npos)
      stringArray[i]= astream.substr(indx, jc-indx);
    else {  // last token
      jc= astream.length()+1;
      stringArray[i]= astream.substr(indx, jc-indx);
    }
    for(std::size_t j=1; jc+j< astream.length(); ++j ) { // skip continuing spaces
      if(astream[G4int(jc+j)]==' ') jc++;
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
  std::size_t maxWidth=0;
  for (G4int iy=1; iy<= GetNRow(icol); iy++) {
    std::size_t ilen= GetElement(icol,iy)->length();
    // care for color code
    // if(GetElement(icol,iy)-> index(strESC,0) != G4String::npos) {
    // if(strESC == (*GetElement(icol,iy))[0] ) {
    const char tgt = (*GetElement(icol,iy))[(std::size_t)0];
    if(strESC == tgt) {
      ilen-= 5;
    }
    if(ilen> maxWidth) maxWidth= ilen;
  }

  return (G4int)maxWidth;
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
      const char tgt = word[(std::size_t)0];
      if(strESC == tgt) {
        colorWord= word.substr(0,5);
        word.erase(0,5);
      }
      if(!colorWord.empty()) G4cout << colorWord << std::flush;

      G4cout << std::setiosflags(std::ios::left) << std::setw(GetNField(ix)) 
             << word.c_str() << std::flush; 
                // against problem w/ g++ iostream
      if(ix != nc) G4cout << "  " << std::flush;
      else G4cout << G4endl;      
    }
  }
}

