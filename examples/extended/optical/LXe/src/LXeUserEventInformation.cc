#include "LXeUserEventInformation.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeUserEventInformation::LXeUserEventInformation()
  :hitCount(0),photonCount_Scint(0),photonCount_Ceren(0),absorptionCount(0),
   boundaryAbsorptionCount(0),totE(0.),eWeightPos(0.),reconPos(0.),convPos(0.),
   convPosSet(false),posMax(0.),pmtsAboveThreshold(0)
{
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeUserEventInformation::~LXeUserEventInformation()
{
}

