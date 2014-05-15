#ifdef WIN32
#include <unordered_map>
#else
#include <tr1/unordered_map>
#define unordered_map tr1::unordered_map
#endif 

int main()
{
  std::unordered_map<int,int> uom(100);
  return 0;
}
