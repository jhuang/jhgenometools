/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <limits.h>
#include "core/ensure.h"
#include "core/safearith.h"
#include "core/unused_api.h"

void gt_safe_default_overflow_handler(const char *src_file, int src_line,
                                      GT_UNUSED void *data)
{
  fprintf(stderr, "%s, l.%d: overflow in operation\n", src_file, src_line);
  exit(EXIT_FAILURE);
}

int gt_safe_abs_check_func(int j, const char *src_file, int src_line,
                      GtOverflowHandlerFunc handler_func, void *data)
{
  int rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

long gt_safe_labs_check_func(long j, const char *src_file, int src_line,
                             GtOverflowHandlerFunc handler_func, void *data)
{
  long rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

long long gt_safe_llabs_check_func(long long j, const char *src_file,
                                   int src_line,
                                   GtOverflowHandlerFunc handler_func,
                                   void *data)
{
  long long rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

uint32_t gt_safe_mult_u32_check_func(uint32_t a, uint32_t b,
                                     const char *src_file, int src_line,
                                     GtOverflowHandlerFunc handler_func,
                                     void *data)
{
  unsigned long long x = (unsigned long long) a * b;
  if (x > 0xffffffff) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return x;
}

uint64_t gt_safe_mult_u64_check_func(uint64_t a, uint64_t b,
                                     const char *src_file, int src_line,
                                     GtOverflowHandlerFunc handler_func,
                                     void *data)
{
  uint32_t a_hi = a >> 32,
           a_lo = a & 0xffffffff,
           b_hi = b >> 32,
           b_lo = b & 0xffffffff;
  /*
     x = 0xffffffff
     a = a_hi*x + a_lo
     b = b_hi*x + b_lo
     a * b = (a_hi*x + a_lo) * (b_hi*x + b_lo)
           = a_hi*x*b_hi*x + a_hi*x*b_lo + a_lo*b_hi*x + a_lo*b_lo
  */
  if (a_hi && b_hi) {    /* overflow */
    handler_func(src_file, src_line, data);
  }
  a = (uint64_t)(a_hi) * b_lo + (uint64_t)(a_lo) * b_hi;
  if (a > 0xffffffff) {  /* overflow */
    handler_func(src_file, src_line, data);
  }
  return (a << 32) + (uint64_t)(a_lo) * b_lo;
}

unsigned long gt_safe_mult_ulong_check_func(unsigned long a, unsigned long b,
                                            const char *src_file, int src_line,
                                            GtOverflowHandlerFunc handler_func,
                                            void *data)
{
  gt_assert(sizeof (unsigned long) == 4 || sizeof (unsigned long) == 8);
  if (sizeof (unsigned long) == 4) {
    return gt_safe_mult_u32_check_func(a, b, src_file, src_line, handler_func,
                                       data);
  } else { /* sizeof (unsigned long) == 8 */
    return gt_safe_mult_u64_check_func(a, b, src_file, src_line, handler_func,
                                       data);
  }
}

long gt_safe_cast2long_check_func(unsigned long value, const char *src_file,
                                  int src_line,
                                  GtOverflowHandlerFunc handler_func,
                                  void *data)
{
  if (value > (~0UL >> 1)) {
    handler_func(src_file, src_line, data);
  }
  return value;
}

unsigned long gt_safe_cast2ulong_check_func(long value, const char *src_file,
                                            int src_line,
                                            GtOverflowHandlerFunc handler_func,
                                            void *data)
{
  if (value < 0) {
    handler_func(src_file, src_line, data);
  }
  return value;
}

int gt_safearith_example(GT_UNUSED GtError *err)
{
  unsigned long ulong;
  long slong;
  unsigned int a, b, c;
  int dest, src;
  gt_error_check(err);

  /* safe assignments */
  slong = 256;
  gt_safe_assign(ulong, slong);
  gt_assert(ulong = 256);

  /* safe additions */
  a = 256;
  b = 1;
  gt_safe_add(c, a, b);
  gt_assert(c == 257);

  /* safe subtractions */
  a = 256;
  b = 1;
  gt_safe_sub(c, a, b);
  gt_assert(c == 255);

  /* safe absolutes */
  src = -256;
  dest = gt_safe_abs(src);
  gt_assert(dest == 256);

  return 0;
}

int gt_safearith_unit_test(GtError *err)
{
  int had_err = 0;
  gt_error_check(err);

  {
    ensure(had_err, __MIN(char) == -128);
    ensure(had_err, __MAX(char) == 127);
    ensure(had_err, __MIN(unsigned char) == 0);
    ensure(had_err, __MAX(unsigned char) == 255);

    ensure(had_err, __MIN(short) == SHRT_MIN);
    ensure(had_err, __MAX(short) == SHRT_MAX);
    ensure(had_err, __MIN(unsigned short) == 0);
    ensure(had_err, __MAX(unsigned short) == USHRT_MAX);

    ensure(had_err, __MIN(int) == INT_MIN);
    ensure(had_err, __MAX(int) == INT_MAX);
    ensure(had_err, __MIN(unsigned int) == 0);
    ensure(had_err, __MAX(unsigned int) == UINT_MAX);

    ensure(had_err, __MIN(long) == LONG_MIN);
    ensure(had_err, __MAX(long) == LONG_MAX);
    ensure(had_err, __MIN(unsigned long) == 0);
    ensure(had_err, __MAX(unsigned long) == ULONG_MAX);

#ifdef LLONG_MIN
    ensure(had_err, __MIN(long long) == LLONG_MIN);
    ensure(had_err, __MAX(long long) == LLONG_MAX);
    ensure(had_err, __MIN(unsigned long long) == 0);
    ensure(had_err, __MAX(unsigned long long) == ULLONG_MAX);
#endif
  }

  {
    unsigned long ulong;
    long slong;

    slong = -1;
    ensure(had_err, assign(ulong, slong));

    ulong = 0;
    slong = LONG_MAX;
    ensure(had_err, !assign(ulong, slong) && ulong == LONG_MAX);

    ulong = ULONG_MAX;
    ensure(had_err, assign(slong, ulong));

    slong = 0;
    ulong = LONG_MAX;
    ensure(had_err, !assign(slong, ulong) && slong == LONG_MAX);
  }

  {
    int x;
    ensure(had_err, add_of(x, INT_MAX, 1));
    ensure(had_err, add_of(x, INT_MAX, 256));
    ensure(had_err, add_of(x, INT_MAX, INT_MAX));

    x = 0; ensure(had_err, !add_of(x, INT_MAX - 1, 1) && x == INT_MAX);
    x = 0; ensure(had_err, !add_of(x, INT_MAX - 256, 256) && x == INT_MAX);
    x = 0; ensure(had_err, !add_of(x, INT_MAX, 0) && x == INT_MAX);

    ensure(had_err, add_of(x, 0x100000000ll, 0x100000000ll));
    x = INT_MAX;
    ensure(had_err, !add_of(x, 0x100000000ll, -0x100000000ll) && x == 0);

    ensure(had_err, sub_of(x, INT_MIN, 1));
    ensure(had_err, sub_of(x, INT_MIN, 256));
    ensure(had_err, sub_of(x, INT_MIN, INT_MAX));

    x = 0; ensure(had_err, !sub_of(x, INT_MIN + 1, 1) && x == INT_MIN);
    x = 0; ensure(had_err, !sub_of(x, INT_MIN + 256, 256) && x == INT_MIN);
    x = 0; ensure(had_err, !sub_of(x, INT_MIN, 0) && x == INT_MIN);
  }

  {
    unsigned int x;
    ensure(had_err, add_of(x, UINT_MAX, 1));
    ensure(had_err, add_of(x, UINT_MAX, 256));
    ensure(had_err, add_of(x, UINT_MAX, UINT_MAX));

    x = 0; ensure(had_err, !add_of(x, UINT_MAX - 1, 1) && x == UINT_MAX);
    x = 0; ensure(had_err, !add_of(x, UINT_MAX - 256, 256) && x == UINT_MAX);
    x = 0; ensure(had_err, !add_of(x, UINT_MAX, 0) && x == UINT_MAX);

    ensure(had_err, add_of(x, 0x100000000ll, 0x100000000ll));
    x = UINT_MAX;
    ensure(had_err, !add_of(x, 0x100000000ll, -0x100000000ll) && x == 0);

    ensure(had_err, sub_of(x, 0, 1));
    ensure(had_err, sub_of(x, 0, 256));
    ensure(had_err, sub_of(x, 0, UINT_MAX));

    x = 0; ensure(had_err, !sub_of(x, 1, 1) && x == 0);
    x = 0; ensure(had_err, !sub_of(x, 256, 256) && x == 0);
    x = 0; ensure(had_err, !sub_of(x, 0, 0) && x == 0);
  }

  {
    int i;
    long l;
    long long ll;

    i = gt_safe_abs(0);
    ensure(had_err, i == 0);

    i = gt_safe_abs(-1);
    ensure(had_err, i == 1);

    i = gt_safe_abs(INT_MIN + 1);
    ensure(had_err, i == INT_MAX);

    l = gt_safe_labs(0);
    ensure(had_err, l == 0);

    l = gt_safe_labs(-1);
    ensure(had_err, l == 1);

    l = gt_safe_labs(LONG_MIN + 1);
    ensure(had_err, l == LONG_MAX);

    ll = gt_safe_llabs(0);
    ensure(had_err, ll == 0);

    ll = gt_safe_llabs(-1);
    ensure(had_err, ll == 1);

#ifdef LLONG_MIN
    ll = gt_safe_llabs(LLONG_MIN + 1);
    ensure(had_err, ll == LLONG_MAX);
#endif
  }

  return had_err;
}
