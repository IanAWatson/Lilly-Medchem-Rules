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

/*
  hash functions for hash maps
*/

#include "iw_stl_hash_map.h"

/*
  This seems to be the fastest one I've found so far
*/

#define VERSION_FROM_SDBM
#ifdef VERSION_FROM_SDBM

size_t
IWStringHash::operator () (const IWString & s) const
{
  size_t rc = 5381;

  const unsigned char * q = (const unsigned char *) s.rawchars ();    // beware alignment

  int n = s.length();

  for (int i = 0; i < n; i++)
  {
    rc = ((rc << 5) + rc) + q[i];
  }

  return rc;
}
#endif

/*
  Found at
  http://www.azillionmonkeys.com/qed/hash.html
  Written by
  Paul Hsieh 

  I played this with stl hash maps, and found that it was always
  slower than the hashes above - not necessarily computing the
  hash function itself, but the overall programme performance.
*/

#ifdef ZILLIONMONKEY

#include <stdint.h> /* Replace with <stdint.h> if appropriate */
#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t
SuperFastHash (const char * data, int len) 
{
  uint32_t hash = len, tmp;
  int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

size_t
IWStringHash::operator () (const IWString & s) const
{
  return SuperFastHash(s.rawchars(), s.length());
}
#endif

#if defined (IW_INTEL_COMPILER)

bool
IWStringHash::operator () (const IWString & s1, const IWString & s2) const
{
  if (s1.length () < s2.length ())
    return true;

  if (s1.length () > s2.length ())
    return false;

// same length, use strcmp

  if (s1.strncmp (s2, s2.length ()) < 0)
    return true;

  return false;
}

#endif
