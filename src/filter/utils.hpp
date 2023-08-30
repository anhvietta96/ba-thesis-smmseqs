#ifndef FILTER_UTILS_HPP
#define FILTER_UTILS_HPP

constexpr size_t constexpr_pow(const size_t base, const size_t exponent){
  return (exponent == 1) ? base : base * constexpr_pow(base,exponent-1);
}

constexpr size_t binom(size_t a,size_t b)
{
  if (b == 0 or a == b)
  {
    return size_t(1);
  }
  if (2*b > a and a != b)
  {
    return binom(a,a-b);
  }
  return binom(a-1,b-1) + binom(a-1,b);
}


#endif