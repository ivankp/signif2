#ifndef IVANP_ARRAY_OPS
#define IVANP_ARRAY_OPS

#include <array>

template <typename T, size_t N, typename X>
inline std::array<T,N>& operator*=(std::array<T,N>& a, X x)
noexcept(noexcept(std::get<0>(a) *= std::declval<X>()))
{
  for (size_t i=0; i<N; ++i) a[i] *= x;
  return a;
}

template <typename T, typename X>
inline bool in(X x, const std::array<T,2>& a)
noexcept(
  noexcept(std::declval<T>() < std::declval<X>()) &&
  noexcept(std::declval<X>() < std::declval<T>()) )
{
  return (a[0] < x && x < a[1]);
}

#ifdef IVANP_ARRAY_BOOST_PO

#include "runtime_exception.hh"

namespace std {
template <typename T, size_t N>
void validate(boost::any& v,
              const std::vector<std::string>& values,
              std::array<T,N>*, int)
{
  namespace po = boost::program_options;
  po::validators::check_first_occurrence(v);
  const std::string& s = po::validators::get_single_string(values);

  std::array<T,N> arr;

  unsigned prev = 0, n = 0;
  for (unsigned i=0; i<s.size(); ++i) {
    if (s[i]==':') {
      arr[n++] = boost::lexical_cast<T>(s.data()+prev,i-prev);
      prev = i+1;
    }
  }
  arr[n++] = boost::lexical_cast<T>(s.data()+prev,s.size()-prev);

  if (n!=N) throw rte(s," does not have ",N," : separated fields");

  v = boost::any(arr);
}
}

#endif

#endif
