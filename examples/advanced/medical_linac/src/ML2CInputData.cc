#include "ML2CInputData.h"

CML2CInputData::CML2CInputData(void)
{
	this->bOnlyVisio=false;

// instantiate the messenger for the general and convergence data
	this->ML2MainMessenger=new CML2MainMessenger(this);
}

CML2CInputData::~CML2CInputData(void)
{
}
