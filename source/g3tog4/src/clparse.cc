// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: clparse.cc,v 1.1 1999-01-07 16:06:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3CLParse
//
// Parse the call List of Geant3 geometry calls and execute them
// in Geant4.
//
// Torre Wenaus, LLNL  6/95
//
// To do:
// - Build a linked List container for tokens rather than fixed array.
//   Not going to be appreciably slower; time goes into I/O.
//
// Comments:
// - RWCString doesn't have toInt, toFloat methods! Had to use C's
//   atol, atof
//
#include "G4ios.hh"
#include <fstream.h>
#include <rw/rstream.h>
#include <rw/ctoken.h>
#include <rw/rwfile.h>

#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"
#include "G3RotTable.hh"
#include "G3PartTable.hh"
#include "G3DetTable.hh"

ofstream ofile;

extern "C" {
#include <stdlib.h>
}

extern ofstream ofile;

// The first volume defined on the call List file is assumed to be
// the mother. 
// If GlobalMotherVolume is non-null, it is used as the mother rather
// than the volume specified on the call List file.
// If GlobalMotherVolume is null, the volume in the call List file is
// used normally.
//
// This is used by BaBar's Bogus, which uses G3toG4 as one way to
// Build subsystems. Bogus must have a way to override the global
// mother. The G3toG4 defined geometry is installed immediately below
// the global mother.

// G4LogicalVolume* GlobalMotherVolume = 0;

// Save the second volume pointer. In Bogus usage, it is the
// containing mother for the subdetector.
// G4LogicalVolume* SubsystemMotherVolume = 0;

G3VolTable G3Vol;   // volume G3 name <-> G4 pointer tables
G3MatTable G3Mat;  // material G3 ID <-> G4 pointer table
G3MedTable G3Med;  // trk media G3 ID <-> G4 pointer table
G3RotTable G3Rot;  // rotation ID <-> G4 transform object table
G3PartTable G3Part; // particle ID <-> ParticleDefinition pointer
G3DetTable G3Det;  // sensitive detector name <-> pointer

G4int narray;

G4int Ipar[1000];
G4double Rpar[1000];
RWCString Spar[1000];

G4int G3CLTokens(RWCString *line, RWCString *tokens);
void G3CLEval(RWCString *tokens, char *select);

// front-end decoders for G3 routines
void PG4gsvolu(RWCString *tokens);
void PG4gspos (RWCString *tokens);
void PG4gsposp(RWCString *tokens);
void PG4gsatt (RWCString *tokens);
void PG4gsrotm(RWCString *tokens);
void PG4gsdvn (RWCString *tokens);
void PG4gsdvt (RWCString *tokens);
void PG4gsdvx (RWCString *tokens);
void PG4gsdvn2(RWCString *tokens);
void PG4gsdvt2(RWCString *tokens);
void PG4gsmate(RWCString *tokens);
void PG4gsmixt(RWCString *tokens);
void PG4gstmed(RWCString *tokens);
void PG4gstpar(RWCString *tokens);
void PG4gspart(RWCString *tokens);
void PG4gsdk  (RWCString *tokens);
void PG4gsdet (RWCString *tokens);
void PG4gsdetv(RWCString *tokens);
void PG4gsdeta(RWCString *tokens);
void PG4gsdeth(RWCString *tokens);
void PG4gsdetd(RWCString *tokens);
void PG4gsdetu(RWCString *tokens);
void PG4ggclos();

void G3CLRead(G4String & fname, char *select = NULL)
//
//  G3CLRead
//  Read the call List file, parse the tokens, and pass the token
//  List to the Geant4 interpreter
//
//  fname: call List filename
//  select: if non-NULL, the selected context. Only the subset of
//          the call List matching the selected context will be
//          executed. If NULL, the full call List will run.
//
{
  RWCString line;
  RWCString tokens[1000];

  const char* ofname = "clparse.out";
  ofile.open(ofname);
  ofile << "Output file open\n";

  G4int count = 0;
  G4int ntokens = 0;
  ifstream istr(fname);
  G4bool _debug=false;
    
  while (line.readLine(istr) && ! istr.eof())
  {
      count++;
      ntokens = G3CLTokens(&line,tokens);  // tokenize the line
      
          // check tokens
      if (_debug) {
          G4cout << "==== Line " << count << " Tokens " << ntokens 
           << " Nvol " << G3Vol.GetEntryCount() << endl;
      
          G4cout << line << " // " << G3Vol.GetEntryCount() << endl;
      }
      
      for (G4int i=0; i < ntokens; i++) ofile << tokens[i] << endl;
      
          // interpret the line as a Geant call
      G3CLEval(tokens, select);
  }
}


G4int G3CLTokens(RWCString *line, RWCString tokens[])
//
// G3CLTokens
//
// Tokenize line, returning tokens in tokens[]. Items in ".."
// are extracted as single tokens, despite embedded spaces.
//
{
    RWCTokenizer next(*line);
    // first tokenize using " to identify strings
    G4int itok = 0;
    G4int ntokens = 0;
    RWCString token1, token2;
    while (!(token1=next("\"")).isNull())
        {
            itok++;
            if (itok%2 == 0 ) // even: inside a string
                {
                    tokens[ntokens++] = token1;
                } else        // not in a quoted string: finish tokenization
                {
                    RWCTokenizer lev2(token1);
                    while (!(token2=lev2()).isNull())
                        {
                            tokens[ntokens] = token2;
                            ntokens++;
                        }
                }
        }
    return ntokens;
}

void G3CLEval(RWCString tokens[], char *select)
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

void G3fillParams(RWCString *tokens, char *ptypes)
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
