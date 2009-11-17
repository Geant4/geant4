#include "ML2Convergence.h"

CML2Convergence::CML2Convergence(G4int seed, G4int saving_in_Selected_Voxels_every_events, G4String FileExperimentalData, G4bool bCompareExp, G4int minNumberOfEvents)
:ML2ExpVoxels(0)
{
	this->bCompareExp=bCompareExp;
	this->fileExperimentalData=FileExperimentalData;

// if the flag compareExp if true and the experimental data is given create the class CML2ExpVoxels
	if (this->bCompareExp && this->fileExperimentalData!="")
	{
		this->ML2ExpVoxels=new CML2ExpVoxels(this->bCompareExp, saving_in_Selected_Voxels_every_events, seed, FileExperimentalData);
		if (!this->ML2ExpVoxels->loadData())
		{
			this->ML2ExpVoxels=0;
		}
	}
	this->minNumberOfEvents=minNumberOfEvents;
}

CML2Convergence::~CML2Convergence(void)
{
	if (this->ML2ExpVoxels!=0)
	{delete this->ML2ExpVoxels;}

}
void CML2Convergence::add(const G4Step* aStep)
{
// accumulate events in the CML2ExpVoxels class (if created)
	if (this->ML2ExpVoxels!=0)
	{
		G4double energyDep=aStep->GetTotalEnergyDeposit();
		if (energyDep>0.)
		{
			G4double density=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();
			this->ML2ExpVoxels->add(aStep->GetPreStepPoint()->GetPosition(), energyDep, density);
		}
	}
}
G4bool CML2Convergence::runAgain()
{
	G4bool bAgain=true;
	if (this->ML2ExpVoxels!=0)
	{
		bAgain=!this->convergenceCriteria();
		return bAgain;
	}
	return bAgain;
}
G4bool CML2Convergence::convergenceCriteria()
{
	if (this->bCompareExp)
	{
		if (this->ML2ExpVoxels->getMinNumberOfEvents() > this->minNumberOfEvents)
		{return false;}
		else
		{return true;}
	}
	return false;
}
