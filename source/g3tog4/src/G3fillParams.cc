// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3fillParams.cc,v 1.2 1999-12-05 17:50:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// extracted from clparse.cc by I. Hrivnacova, 30.7.99

#include "globals.hh"
#include "g4std/fstream"

ofstream ofile;

G4int narray;
G4int Ipar[1000];
G4double Rpar[1000];
G4String Spar[1000];

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
