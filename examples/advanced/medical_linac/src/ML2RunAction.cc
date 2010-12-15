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
