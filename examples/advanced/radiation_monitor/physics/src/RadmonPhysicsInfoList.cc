//
// File name:     RadmonPhysicsInfoList.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfoList.cc,v 1.1 2005-11-10 08:14:10 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsInfoList.hh"

                                                RadmonPhysicsInfoList :: RadmonPhysicsInfoList(const RadmonPhysicsInfoList & copy)
:
 infoVector(copy.infoVector)
{
}





RadmonPhysicsInfoList &                         RadmonPhysicsInfoList :: operator=(const RadmonPhysicsInfoList & copy)
{
 infoVector=copy.infoVector;
 
 return (*this);
}





G4bool                                          RadmonPhysicsInfoList :: CollidesWith(const RadmonPhysicsInfoList & other) const
{
 InfoVector::const_iterator j;
 InfoVector::const_iterator i(infoVector.begin());
 const InfoVector::const_iterator endJ(other.infoVector.end());
 const InfoVector::const_iterator endI(infoVector.end());
 
 while (i!=endI)
 {
  j=other.infoVector.begin();
  
  while (j!=endJ)
  {
   if (i->CollidesWith(*j))
    return true;
    
   j++;
  }
  
  i++;
 }
 
 return false;
}





void                                            RadmonPhysicsInfoList :: InsertPhysicsInfo(const RadmonPhysicsInfo & info)
{
 infoVector.push_back(info);
}





std::ostream &                                  operator<<(std::ostream & out, const RadmonPhysicsInfoList & infoList)
{
 const G4int n(infoList.GetNPhysicsInfos());
 
 for (G4int i(0); i<n; i++)
 {
  if (i>0)
   out << ", ";
   
  out << infoList.GetPhysicsInfo(i);
 }
 
 return out;
}
