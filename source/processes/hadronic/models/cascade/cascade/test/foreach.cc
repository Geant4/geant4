#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

void initRadius(int *elem){
  *elem = 1;
}

void printItem(int *elem){
  cout << *elem << endl;
}

int main(){
  vector<int*> test;
  for(int i=0; i < 3; i++) test.push_back(new int(i)); 

  for_each(test.begin(), test.end(), initRadius);
  for_each(test.begin(), test.end(), printItem);

  return 0;
}
