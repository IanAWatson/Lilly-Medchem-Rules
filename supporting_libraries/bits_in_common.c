/**************************************************************************

    Copyright (C) 2011  Eli Lilly and Company

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*#define CHECK_RESULTS*/

#if defined (CHECK_RESULTS)

static const unsigned char one_bit_8[8] = {128, 64, 32, 16, 8, 4, 2, 1};

static const unsigned int one_bit_32[32] = {
  0x80000000,
  0x40000000,
  0x20000000,
  0x10000000,
  0x08000000,
  0x04000000,
  0x02000000,
  0x01000000,
  0x00800000,
  0x00400000,
  0x00200000,
  0x00100000,
  0x00080000,
  0x00040000,
  0x00020000,
  0x00010000,
  0x00008000,
  0x00004000,
  0x00002000,
  0x00001000,
  0x00000800,
  0x00000400,
  0x00000200,
  0x00000100,
  0x00000080,
  0x00000040,
  0x00000020,
  0x00000010,
  0x00000008,
  0x00000004,
  0x00000002,
  0x00000001 };


static int
check_bits_in_common (unsigned int b1, unsigned int b2,
                      int other_result)
{
  int bic = 0;
  unsigned int band = b1 & b2;
  int i;

  for (i = 0; i < 32; i++)
  {
    if (one_bit_32[i] & band)
      bic++;
  }

  if (bic == other_result)
    return 1;

  fprintf (stderr, "Warning, mismatch %u (%x) %u (%x) and %x, bic = %d, other = %d\n",
           b1, b1, b2, b2, band, bic, other_result);

  return 0;
}

#endif

/*
   Sun   METHOD1    sun compiler
   Linux METHOD4    Intel compiler
   SGI   METHOD4    SGI compiler
*/

#define BIC_METHOD4

#if defined (BIC_METHOD1)

#include "precompbit8.h"

#define IW_LAST_BYTE 0x000000ff

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  unsigned int c;

  for (i = 0; i < nwords; i++)
  {
#if defined(CHECK_RESULTS)
    int tmp;
#endif

    c = *b1 & *b2;

#if defined(CHECK_RESULTS)
    tmp = eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
    check_bits_in_common (*b1, *b2, tmp);
    rc += tmp;
#else
    rc += eight_bit_count[c & IW_LAST_BYTE] + eight_bit_count[(c >> 16) & IW_LAST_BYTE] + eight_bit_count[(c >> 8) & IW_LAST_BYTE] + eight_bit_count[c >> 24];
#endif

    b1++;
    b2++;
  }

  return rc;
}

#elif defined (BIC_METHOD2)

#include "precompbit.h"

static unsigned int andtmp[512];

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int nsw = nwords * 2;
  int rc = 0;
  const unsigned short * s = (const unsigned short *) &andtmp;

  for (i = 0; i < nwords; i++)
  {
    andtmp[i] = (b1[i]) & (b2[i]);
  }

#ifdef CHECK_AND
  for (i = 0; i < nwords; i++)
  {
    unsigned int j = b1[i] & b2[i];
    if (j != andtmp[i])
    {
      fprintf (stderr, "Vectorisation failure, i = %d, vector %u actually %u\n", i, andtmp[i], j);
    }
  }
#endif

  for (i = 0; i < nsw; i++)
  {
    rc += sixteen_bit_count[s[i]];
  }

  return rc;
}

#elif defined (BIC_METHOD3)

#include "precompbit.h"

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  int ns = nwords * 2;     /* 2 shorts per word */
  const unsigned short * sb1 = (const unsigned short *) b1;
  const unsigned short * sb2 = (const unsigned short *) b2;
  unsigned short c;

  for (i = 0; i < ns; i++)
  {
    c = *sb1 & *sb2;
    sb1++;
    sb2++;
    rc += sixteen_bit_count[c];
  }

  return rc;
}

#elif defined (BIC_METHOD4)

#include "precompbit.h"

struct Uc2
{
  unsigned short s0;
  unsigned short s1;
};

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int j;
  int rc = 0;

  union And_Target {
    unsigned int uint;
    struct Uc2 uc2;
  } and_target;

  for (j = 0; j < nwords; j++)
  {
#if defined (CHECK_RESULTS)
    int tmp;
#endif

    and_target.uint = (*b1 & *b2);
    b1++;
    b2++;

#ifndef __i386__
    /* found this slowed things down on Intel */
    if (0 == and_target.uint)
      continue;
#endif

#if defined(CHECK_RESULTS)
    tmp = sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
    check_bits_in_common (*(b1 - 1), *(b2 - 1), tmp);
    rc += tmp;
#else
    rc += sixteen_bit_count[and_target.uc2.s0] + 
          sixteen_bit_count[and_target.uc2.s1];
#endif
  }

  return rc;
}

#elif defined (BIC_METHOD5)

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int i;
  int rc = 0;
  unsigned int x;

  for (i = 0; i < nwords; i++)
  {
    x = b1[i] & b2[i];

    x -= (x >>1) & 0x55555555;
    x  = ((x >> 2) & 0x33333333) + (x & 0x33333333);
    x  = ((x >> 4) + x) & 0x0f0f0f0f0f;
    x *= 0x01010101;

    rc += x >> 24;
  }

  return rc;
}


#elif defined (BIC_METHOD6)

#include "precompbit8.h"

struct Uc4
{
  unsigned char c0;
  unsigned char c1;
  unsigned char c2;
  unsigned char c3;
};

int
bits_in_common (const unsigned int * b1, const unsigned int * b2,
                const int nwords)
{
  int j;
  int rc = 0;

  union And_Target {
    unsigned int uint;
    struct Uc4 uc4;
  } and_target;

  for (j = 0; j < nwords; j++)
  {
#if defined (CHECK_RESULTS)
    int tmp;
#endif

    and_target.uint = (*b1 & *b2);
    b1++;
    b2++;

#ifndef __i386__
    /* found this slowed things down on Intel */
    if (0 == and_target.uint)
      continue;
#endif

#if defined(CHECK_RESULTS)
    tmp = eight_bit_count[and_target.uc4.c0] + 
          eight_bit_count[and_target.uc4.c1] +
          eight_bit_count[and_target.uc4.c2] +
          eight_bit_count[and_target.uc4.c3];
    check_bits_in_common (*(b1 - 1), *(b2 - 1), tmp);
    rc += tmp;
#else
    rc += eight_bit_count[and_target.uc4.c0] + 
          eight_bit_count[and_target.uc4.c1] +
          eight_bit_count[and_target.uc4.c2] +
          eight_bit_count[and_target.uc4.c3];
#endif
  }

  return rc;
}

#else

#error "Must define a BIC_METHOD symbol"

#endif
