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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2Convergence.hh"

CML2Convergence::CML2Convergence(void)
{}

CML2Convergence::CML2Convergence(G4int seed, G4int saveEvents,
              G4String FileExperimentalData, G4String FileExperimentalDataOut,
              G4bool bComp, G4int maxNumEvents, G4int nRecycling, G4int maxLoops)
    :ML2ExpVoxels(0)
{
    nGeometry = 0;
    nMaxLoops = maxLoops;
    idCurrentLoop = nMaxLoops;
    bCompareExp = bComp;
    nAccumulatedEvents = 0;
    if (bCompareExp){nMaxLoops =-1;};
    fileExperimentalData = FileExperimentalData;

    // if the flag compareExp if true and the experimental data is given, create the class CML2ExpVoxels
    if (bCompareExp && fileExperimentalData!="")
    {
        ML2ExpVoxels = new CML2ExpVoxels(bCompareExp, saveEvents, seed, FileExperimentalData, FileExperimentalDataOut);
        if (!ML2ExpVoxels->loadData())
        {
            nMaxLoops =10;
            ML2ExpVoxels=0;
            G4cout << "I don't have any convergence criteria set, I'll do " << nMaxLoops << " loop(s) for each rotation" << G4endl;
        }
        else
        {
            ML2ExpVoxels->setRecycling(nRecycling);
        }
    }
    maxNumberOfEvents=maxNumEvents;
}

CML2Convergence::~CML2Convergence(void)
{
    if (ML2ExpVoxels!=0)
    {delete ML2ExpVoxels;}

}
void CML2Convergence::add(const G4Step* aStep)
{
    // accumulate events in the CML2ExpVoxels class (if created)
    if (ML2ExpVoxels!=0)
    {
        if (aStep->GetTotalEnergyDeposit()>0.)
        {ML2ExpVoxels->add(aStep);}
    }
}
G4bool CML2Convergence::stopRun()
{
    G4bool bStopRun=false;
    if (ML2ExpVoxels!=0) // true if the experimental data file exists and is used to check the convergence
    {
        bStopRun=convergenceCriteria();
        return bStopRun;
    }
    else // true if no experiemental data file is used. In this case it runs "nMaxLoops" loops.
    {
        idCurrentLoop--;
        if (idCurrentLoop==0)
        {
            bStopRun=true;
        }
    }
    return bStopRun;
}
G4bool CML2Convergence::convergenceCriteria()
{
    G4bool bStopRun=true;
    G4int nEventsAccumulated=0;
    if (bCompareExp)
    {
        nEventsAccumulated=ML2ExpVoxels->getMaxNumberOfEvents();
        // It checks if the maximum number of events is reached at least in one voxel. Having more rotations the limits is incremented each rotation
        if (ML2ExpVoxels->getMaxNumberOfEvents()>= maxNumberOfEvents)
        {bStopRun = true; ML2ExpVoxels->resetNEventsInVoxels();}
        else
        {bStopRun = false;}
    }
    G4cout <<"\n ++++++++++++++++++++ "  << G4endl;
    G4cout <<"current geometry: " << nGeometry;
    G4cout << "\nNumber of events accumulated in the current geometry:"<<
                 nEventsAccumulated<<"\nNumber of events to be accumulated:" <<
                maxNumberOfEvents<< "\n -------------------------\n" << G4endl;
    return bStopRun;
}
