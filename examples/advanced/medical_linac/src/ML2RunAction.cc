//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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


#include "ML2RunAction.hh"

CML2RunAction::CML2RunAction(CML2Convergence *convergence, G4int nBeam, G4bool bOnlyVisio)
{
	this->bRotationTranslationFilesNames=true;
	this->convergence=convergence;
	this->nBeam=nBeam;
	this->bOnlyVisio=bOnlyVisio;
}

CML2RunAction::~CML2RunAction(void)
{
}
void CML2RunAction::BeginOfRunAction(const G4Run *)
{
	G4String fullName;
	if (this->bRotationTranslationFilesNames)
	{fullName=CML2AcceleratorConstruction::GetInstance()->getCurrentRotationString()+CML2PhantomConstruction::GetInstance()->getCurrentTranslationString();}
	else
	{fullName="";}
	CML2PhantomConstruction::GetInstance()->setNewName(fullName);

	CML2AcceleratorConstruction::GetInstance()->writeInfo();
	CML2PhantomConstruction::GetInstance()->writeInfo();

	std::cout<<"*********************************************"<<'\n';
	if (this->convergence->getNMaxLoops()<0 || this->bOnlyVisio)
	{	
		std::cout << "loop n. "<<++this->nLoop <<'\n';
	}
	else
	{
		std::cout << "loop n. "<<++this->nLoop<<"/" <<this->convergence->getNMaxLoops() <<'\n';
	}
	if (!this->bOnlyVisio)
	{std::cout << "Launched "<< this->nBeam <<" random primary particles" << '\n';}
	std::cout<<"*********************************************"<<'\n';
	this->MyTime.Start();
}
void CML2RunAction::EndOfRunAction(const G4Run *)
{
	CML2WorldConstruction::GetInstance()->savePhantomData();
	CML2WorldConstruction::GetInstance()->savePhaseSpaceData();
	this->convergence->saveResults();

	this->MyTime.Stop();
	this->loopElapsedTime=MyTime.GetUserElapsed();
	std::cout << "loop elapsed time [s] : "<< this->loopElapsedTime << '\n';
	std::cout <<'\n';
}
