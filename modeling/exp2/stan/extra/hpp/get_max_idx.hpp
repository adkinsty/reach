#include <algorithm>
using namespace std; 
template <typename T0__>
int get_max_idx(const std::vector<T0__>& EU, std::ostream* pstream__) {
  return std::max_element(EU.begin(),EU.end()) - EU.begin() + 1;
}