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
#include <ctime>
#include "LaunchG4.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char* argv[])
{
    G4int seed =  0;
    G4bool g4session = false;

    G4String macFile;
    G4bool useCustomMacFile = false;

#ifdef G4VIS_USE
    G4bool drawing = false;
#endif

    G4double incidentEnergy = 0;
    G4String strIncidentEnergy;
    G4bool chemistryflag = true;

    if (argc>0)
    {
        static char null[1] = { "" };
        for (G4int i = 1 ; i <argc ; i++)
        {
            if(!strcmp(argv[i],"-mac"))
            {
                macFile = argv[i+1];
                useCustomMacFile = true;
                argv[i]=null;
                argv[i+1]=null;
                i++;
            }

            else if (!strcmp(argv[i],"-e"))
            {
                incidentEnergy = atoi(argv[i+1]);
                strIncidentEnergy = argv[i+1];
                argv[i]=null;
                argv[i+1]=null;
                i++;
            }

            else if (!strcmp(argv[i],"-seed"))
            {
                seed = atoi(argv[i+1]);
                argv[i]=null;
                argv[i+1]=null;
                i++;
            }

            else if(!strcmp(argv[i],"-g4session"))
            {
                g4session = true;
                argv[i]=null;
            }

            else if(!strcmp(argv[i],"-nochem"))
            {
                chemistryflag = false;
                argv[i]=null;
            }

            else if(!strcmp(argv[i],"&"))
            {
                argv[i]=null;
            }
#ifdef G4VIS_USE
            else if (!strcmp(argv[i],"-drawing"))
            {
                drawing = true;
                g4session = true;
                argv[i]=null;
            }
#endif
        }

        // remove handled arguments from argument array
        int j = 0;
        for (int i = 0; i < argc; i++)
        {
            if (strcmp(argv[i], ""))
            {
                argv[j] = argv[i];
                j++;
            }
        }
        argc = j;
    }

    if(incidentEnergy==0)
    {
        strIncidentEnergy = "750"; //keV
        incidentEnergy = atoi(strIncidentEnergy);
    }

    if (argc>0)
    {
        G4bool kill = false;
        for (G4int i = 1 ; i <argc ; i++)
        {
            if (strcmp(argv[i], ""))
            {
                kill=true;
                G4cerr<<"Unknown argument : "<< argv[i] <<"\n";
            }
        }
        if(kill) exit(1);
    }

    // Choose the Random engine
    // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

    LaunchG4 *g4 = new LaunchG4();
    g4->Initialize(incidentEnergy, chemistryflag);

    if(seed==0)
    {
        seed = ((int) time(NULL)) ;
        CLHEP::HepRandom::setTheSeed(seed);

    }
    else
    {
        CLHEP::HepRandom::setTheSeed(seed);
    }

    G4cout<<"Seed used : "<< seed <<G4endl;

    if (g4session || drawing)
    {
        g4 -> NewSession(argc,argv);
    }

    if(useCustomMacFile == false)
    {
#ifdef G4VIS_USE
    g4 -> RunSimu(drawing);
#else
    g4 -> RunSimu();
#endif
    }
    else
    {
#ifdef G4VIS_USE
    g4 -> RunSimu(macFile,drawing);
#else
    g4 -> RunSimu(macFile);
#endif
    }

    if (g4session || drawing)
    {
        g4 -> StartSession();
    }

//    system("rm output.txt");
    delete g4;
    return 0;
}

