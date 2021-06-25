// Check that trivial use of <filesystem> works
// - https://en.cppreference.com/w/cpp/filesystem
#include <iostream>
#if __has_include( <filesystem>)
#  include <filesystem>
namespace fs = std::filesystem;
#elif __has_include( <experimental/filesystem>)
#  include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

int main( int argc, char* argv[] ) {
  fs::path p{ argv[0] };
  std::cout << p << ", " << fs::canonical( p ) << std::endl;
}