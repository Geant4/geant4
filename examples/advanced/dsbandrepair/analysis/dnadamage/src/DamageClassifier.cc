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
/// \file DamageClassifier.cc
/// \brief Implementation of the DamageClassifier class

#include "DamageClassifier.hh"

#include <iostream>
#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<ClassifiedDamage> DamageClassifier::MakeCluster(
	std::vector<Damage>& pListDamage,unsigned int pDSBLength, bool pBase)
{
	unsigned long int copyNb;
	unsigned long int lastCopyNb = -pDSBLength-1;

	// sort the list of damage by ascending bp
	std::sort(pListDamage.begin(), pListDamage.end(),
		[](const Damage& a, const Damage& b) {
			return a.GetCopyNb() < b.GetCopyNb();
		});

	std::vector<ClassifiedDamage> listClassifiedDamage;
	ClassifiedDamage classDamage;

	for(auto it=pListDamage.begin();it!=pListDamage.end();it++)
	{
		classDamage.SetIncludeBase(pBase);
		Damage tempDamage = (*it);

		if(tempDamage.GetDamageType()==Damage::DamageType::fBackbone)
		{

			copyNb=tempDamage.GetCopyNb();

			if(classDamage.GetNumDamage()<=0)
			{
				classDamage.AddDamage(tempDamage);
				lastCopyNb = copyNb;
			}
			else
			{
				// New Damage
				if(copyNb>lastCopyNb+pDSBLength)
				{
					classDamage.ComputeBp();
					classDamage.ComputeType();
					classDamage.ComputeComplexity();
					listClassifiedDamage.push_back(classDamage);
					classDamage.Reset();
				}

				classDamage.AddDamage(tempDamage);
				lastCopyNb = copyNb;
			}
		}
	}

	if(classDamage.GetNumDamage()>0)
	{
		classDamage.ComputeBp();
		classDamage.ComputeType();
		classDamage.ComputeComplexity();
		listClassifiedDamage.push_back(classDamage);
	}

	// TODO: include base damage if include base (pBase) is true

	return listClassifiedDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumSSB(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numSSB = 0;
	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fSSB)
		{
			numSSB++;
		}
	}
	return numSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumDSB(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numDSB = 0;
	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fDSB)
		{
			numDSB++;
		}
	}
	return numDSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumComplexDSB(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numComplexDSB = 0;

	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(((it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fDSB))&&(it->GetComplexity()>1))
		{
			numComplexDSB++;
		}
	}
	return numComplexDSB;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > DamageClassifier::SortDamageByEvent(
	const std::vector<Damage>& pListDamage)
{
	std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > mapDamage;
	for(auto it=pListDamage.begin();it!=pListDamage.end();it++)
	{
		mapDamage[it->GetEvt()][it->GetChromo()].push_back((*it));
	}
	return mapDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > DamageClassifier::SortDamageByChromo(
	const std::vector<Damage>& pListDamage)
{
	std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > mapDamage;
	for(auto it=pListDamage.begin();it!=pListDamage.end();it++)
	{
		mapDamage[it->GetChromo()][it->GetEvt()].push_back((*it));
	}
	return mapDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumDSBwithDirectDamage(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numDSBwDir = 0;
	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fDSB && 
			it->GetIsThereDirectComponentContribution())
		{
			numDSBwDir++;
		}
	}
	return numDSBwDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumDSBwithIndirectDamage(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numDSBwIn = 0;
	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fDSB && 
			it->GetIsThereIndirectComponentContribution())
		{
			numDSBwIn++;
		}
	}
	return numDSBwIn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int DamageClassifier::GetNumDSBwithBothDirectIndirectDamage(
	const std::vector<ClassifiedDamage>& pListClassifiedDamage) const
{
	unsigned int numDSBwDirIn = 0;
	for(auto it=pListClassifiedDamage.begin();it!=pListClassifiedDamage.end();it++)
	{
		if(it->GetClassifiedDamageType()==ClassifiedDamage::ClassifiedDamageType::fDSB && 
			it->GetIsThereDirectComponentContribution() && 
			it->GetIsThereIndirectComponentContribution())
		{
			numDSBwDirIn++;
		}
	}
	return numDSBwDirIn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......