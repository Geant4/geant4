#ifdef WIN32
#include <unordered_map>
#define UOM "UOM"
#elif (defined(__GNUC__) && ((__GNUC__==4 && __GNUC_MINOR__>=1) || __GNUC__>4 ))
#include <tr1/unordered_map>
#define unordered_map tr1::unordered_map
#define UOM "UOM"
#else
#include <map>
#define unordered_map map
#define UOM "MAP"
#endif 

#include <iostream>

void Print(const char* p, int i )
{
  std::cout << p << " = " << i << std::endl;
}

int main()
{
  std::unordered_map<int,int> uom;

#ifdef __GNUC__
  Print("__GNUC__",__GNUC__);
  Print("__GNUC_MINOR__",__GNUC_MINOR__);
#endif

  std::cout << UOM << std::endl;
  return 0;
}
