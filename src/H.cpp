#include <cpp11.hpp>
using namespace cpp11;

[[cpp11::register]]
double Hcpp(int n, cpp11::doubles P, cpp11::doubles h)
{
  if (n > 1)
    return P[n - 1] * h[n - 1] + (1 - P[n - 1]) * Hcpp(n - 1, P, h);
  else
    return h[0];
}
