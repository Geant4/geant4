#include "G4HadTmpUtil.hh"
#include <vector>
#include <iostream>
using namespace std;

int main()
{
  vector<int *> it;
  for(int i=10; i>0; i--)   it.push_back(new int(i));
  for(size_t i=0; i<it.size(); i++)
  {
    cout << *it[i]<<endl;
  }
  G4PtrSort<int>(&it);
  cout << "And then..."<<endl;
  for(size_t i=0; i<it.size(); i++)
  {
    cout << *it[i]<<endl;
  }
}
