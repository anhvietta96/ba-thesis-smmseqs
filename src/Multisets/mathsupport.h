#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include "minmax.h"

static inline bool gt_float_isclose(float a, float b,float rel_tol,
                                    float abs_tol)
{
  const float max_value = rel_tol * MAX(fabsf(a),fabsf(b));
  return fabsf(a - b) <= MAX(max_value,abs_tol) ? true : false;
}


static inline unsigned long gt_required_bits2(unsigned long value)
{
#if __has_builtin(__builtin_clz)
  return value == 0 ? 0 : (sizeof(value) * CHAR_BIT - __builtin_clz(value));
#else
  for(unsigned long count = 0; value > 0; count++)
  {
    value >>= 1;
  }
  return count;
#endif
}

static inline unsigned long gt_required_bits(unsigned long value)
{
  return value == 0 ? 0 : ((unsigned long) ceil(log2((double) value)));
}

static inline unsigned long ulong_pow(unsigned long base,unsigned long exponent)
{
  unsigned long idx, powvalue = 1UL;

  for (idx = 0; idx < exponent; idx++)
  {
    powvalue *= base;
  }
  return powvalue;
}

static inline int code_bits_from_alphasize(size_t alphasize,
                                           unsigned long qgram_length)
{
  unsigned long num_different_qgrams = ulong_pow(alphasize,qgram_length);
  return (int) gt_required_bits(num_different_qgrams);
}
#endif
