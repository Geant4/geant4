//
// SBTVisManager
//
// Definition of visualization manager for SBT
//

#ifndef SBTVisManager_H
#define SBTVisManager_H

#include "G4VisManager.hh"
#include "G4NullModel.hh"

class SBTFakeModel : public G4NullModel {
	public:
	SBTFakeModel(const G4ModelingParameters* pMP) : G4NullModel(pMP) { ; }
	~SBTFakeModel() {;}
	
	inline virtual void DescribeYourselfTo( G4VGraphicsScene &scene ) { ; }
};
		

class SBTVisManager : public G4VisManager
{
	public:
	SBTVisManager() {}
	~SBTVisManager() {}
	
	G4int BuildFakeWorld() const;
	
	private:
	void RegisterGraphicsSystems();
}; 

#endif
