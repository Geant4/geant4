// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: clparse.cc,v 1.8 1999-12-05 17:50:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I.Hrivnacova
// added G3SensVol

#include "globals.hh"
#include "g4std/fstream"
#include "g4rw/rstream.h"
#include "g4rw/ctoken.h"
#include "G3toG4.hh"
#include "G3EleTable.hh"
#include "G3VolTable.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"
#include "G3RotTable.hh"
#include "G3PartTable.hh"
#include "G3DetTable.hh"
#include "G3SensVolVector.hh"

ofstream ofile;

extern "C" {
#include <stdlib.h>
}

extern ofstream ofile;

G3VolTable G3Vol;
G3MatTable G3Mat; // material G3 ID <-> G4 pointer table
G3MedTable G3Med; // trk media G3 ID <-> G4 pointer table
G3RotTable G3Rot; // rotation ID <-> G4 transform object table
G3PartTable G3Part; // particle ID <-> ParticleDefinition pointer
G3DetTable G3Det; // sensitive detector name <-> pointer
G3EleTable G3Ele; // element names table
G3SensVolVector G3SensVol; // vector of sensitive logical volumes

G4int narray;

G4int Ipar[1000];
G4double Rpar[1000];
G4String Spar[1000];

G4int G3CLTokens(G4String *line, G4String *tokens);
void G3CLEval(G4String *tokens, char *select);

// front-end decoders for G3 routines
void PG4gsvolu(G4String *tokens);
void PG4gspos (G4String *tokens);
void PG4gsposp(G4String *tokens);
void PG4gsatt (G4String *tokens);
void PG4gsrotm(G4String *tokens);
void PG4gsdvn (G4String *tokens);
void PG4gsdvt (G4String *tokens);
void PG4gsdvx (G4String *tokens);
void PG4gsdvn2(G4String *tokens);
void PG4gsdvt2(G4String *tokens);
void PG4gsmate(G4String *tokens);
void PG4gsmixt(G4String *tokens);
void PG4gstmed(G4String *tokens);
void PG4gstpar(G4String *tokens);
void PG4gspart(G4String *tokens);
void PG4gsdk  (G4String *tokens);
void PG4gsdet (G4String *tokens);
void PG4gsdetv(G4String *tokens);
void PG4gsdeta(G4String *tokens);
void PG4gsdeth(G4String *tokens);
void PG4gsdetd(G4String *tokens);
void PG4gsdetu(G4String *tokens);
void PG4ggclos();

void G3CLRead(G4String & fname, char *select = NULL){
//
//  G3CLRead
//  Read the call List file, parse the tokens, and pass the token
//  List to the Geant4 interpreter
//
//  fname: call List filename
//
  
  G4String line;
  G4String tokens[1000];

  const char* ofname = "clparse.out";
  ofile.open(ofname);
  ofile << "Output file open\n";

  G4int count = 0;
  G4int ntokens = 0;
  ifstream istr(fname);
  G4bool _debug=false;
    
  while (line.readLine(istr) && ! istr.eof()){
      count++;
      ntokens = G3CLTokens(&line,tokens);  // tokenize the line
      for (G4int i=0; i < ntokens; i++) ofile << tokens[i] << endl;
      
          // interpret the line as a Geant call
      G3CLEval(tokens, select);
  }
}


G4int G3CLTokens(G4String *line, G4String tokens[])
//
// G3CLTokens
//
// Tokenize line, returning tokens in tokens[]. Items in ".."
// are extracted as single tokens, despite embedded spaces.
//
{
    G4Tokenizer next(*line);
    // first tokenize using " to identify strings
    G4int itok = 0;
    G4int ntokens = 0;
    G4String token1, token2;
    while (!(token1=next("\"")).isNull())
        {
            itok++;
            if (itok%2 == 0 ) // even: inside a string
                {
                    tokens[ntokens++] = token1;
                } else        // not in a quoted string: finish tokenization
                {
                    G4Tokenizer lev2(token1);
                    while (!(token2=lev2()).isNull())
                        {
                            tokens[ntokens] = token2;
                            ntokens++;
                        }
                }
        }
    return ntokens;
}

void G3CLEval(G4String tokens[], char *select)
//
// G3CLEval
//
// Evaluate the token List as a Geant3 call, and execute it as
// a Geant4 call.
//
{
    const char* context = tokens[0];
    const char* routine = tokens[1];

    // If context is selected, return unless context matches.
    if (select != NULL && select != "*") if ( strcmp(select,context) ) return;

    // Branch on Geant3 routine name
    ofile << "Do routine " << routine << " in context " << context << endl;
    
    if ( !strcmp(routine,"GSVOLU") ) { 
//      volcount++;
//      if (volcount == 1) {
//        // Special handling of the first one, assumed to be global mother
//        if ( GlobalMotherVolume == 0 ) {
//        PG4gsvolu(&tokens[2]); 
//        } else {
//          G4String gblmoth="Global mother";
//          G3Vol.PutLV(&gblmoth,GlobalMotherVolume);
//        }
//      } else {
//        PG4gsvolu(&tokens[2]); 
//      }
//      if (volcount == 2) {
//        G4String vname = tokens[2];
//        SubsystemMotherVolume = G3Vol.GetLV(&vname);
//      }
        { PG4gsvolu(&tokens[2]); return;}
    }
    if ( !strcmp(routine,"GSPOS") )  { PG4gspos (&tokens[2]); return;}
    if ( !strcmp(routine,"GSPOSP") ) { PG4gsposp(&tokens[2]); return;}
    if ( !strcmp(routine,"GSATT") )  { PG4gsatt (&tokens[2]); return;}
    if ( !strcmp(routine,"GSROTM") ) { PG4gsrotm(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDVN") )  { PG4gsdvn (&tokens[2]); return;}
    if ( !strcmp(routine,"GSDVT") )  { PG4gsdvt (&tokens[2]); return;}
    if ( !strcmp(routine,"GSDVX") )  { PG4gsdvx (&tokens[2]); return;}
    if ( !strcmp(routine,"GSDVN2") ) { PG4gsdvn2(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDVT2") ) { PG4gsdvt2(&tokens[2]); return;}
    if ( !strcmp(routine,"GSMATE") ) { PG4gsmate(&tokens[2]); return;}
    if ( !strcmp(routine,"GSMIXT") ) { PG4gsmixt(&tokens[2]); return;}
    if ( !strcmp(routine,"GSTMED") ) { PG4gstmed(&tokens[2]); return;}
    if ( !strcmp(routine,"GSTPAR") ) { PG4gstpar(&tokens[2]); return;}
    if ( !strcmp(routine,"GSPART") ) { PG4gspart(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDK") )   { PG4gsdk  (&tokens[2]); return;}
    if ( !strcmp(routine,"GSDET") )  { PG4gsdet (&tokens[2]); return;}
    if ( !strcmp(routine,"GSDETV") ) { PG4gsdetv(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDETA") ) { PG4gsdeta(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDETH") ) { PG4gsdeth(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDETD") ) { PG4gsdetd(&tokens[2]); return;}
    if ( !strcmp(routine,"GSDETU") ) { PG4gsdetu(&tokens[2]); return;}
    if ( !strcmp(routine,"GGCLOS") ) { PG4ggclos(); return;}
}

void G3fillParams(G4String *tokens, char *ptypes)
//
// G3fillParams
//
// Interpret tokens to fill call parameters, based on parameter
//   types ptypes
//
{
    // loop over ptypes
    G4int i =0, ipt = 0, k = 0;
    G4int ni =0, nr = 0, ns = 0;
    while (ptypes[i] != '\0')
        {
            switch (ptypes[i]) {
            case 'i':
                Ipar[ni] = atoi(tokens[ipt].data());
                narray = Ipar[ni];
                ni++; ipt++;
                break;
            case 'r':
                Rpar[nr] = atof(tokens[ipt].data());
                nr++; ipt++;
                break;
            case 's':
                Spar[ns] = tokens[ipt];
                ns++; ipt++;
                break;
            case 'I':
                for (k=0; k < narray; k++)
                    {
                        Ipar[ni] = atoi(tokens[ipt].data());
                        ni++; ipt++;
                    }
                break;
            case 'R':
                for (k=0; k < narray; k++) 
                    { 
                        Rpar[nr] = atof(tokens[ipt].data()); 
                        nr++; ipt++;
                    }
                break;
            case 'Q':
                // special case of reading three successive R arrays 
                // into one (used in gsmixt)
                narray = 3 * abs(narray);
                for (k=0; k < narray; k++) 
                    { 
                        Rpar[nr] = atof(tokens[ipt].data()); 
                        nr++; ipt++;
                    }
                break;
            case 'S':
                for (k=0; k < narray; k++)
                    {
                        Spar[ns] = tokens[ipt];
                        ns++; ipt++;
                    }
                break;
            default:
                ofile << "unidentified ptype '" << ptypes[i] << endl;
            };
            i++;
        }
}
