/**************************************************************************

    Copyright (C) 2012  Eli Lilly and Company

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
/*
  Random number stuff
*/

#include <stdlib.h>
#include <sys/types.h>
#ifdef __WIN32__
#include <process.h>
#else
#include <unistd.h>
#endif

#include <time.h>
#include <math.h>
#include <stdio.h>

#include <assert.h>

#include <iostream>
#include <limits>
//#include <iomanip>

using namespace std;

#include "iwrandom.h"

#include "iwbits.h"

/*
  We initialise the last two bytes with an arbitrary value
  We rely on memory allocation to give a random value for the
  rest. No real reason why we couldn't just rely on memory 
  allocation for the whole thing....
*/

void
Random_Number_Working_Storage::_default_values ()
{
  _initialised = 0;

  _buffer[4] = 0x33;
  _buffer[5] = 0x0E;

  return;
}

Random_Number_Working_Storage::Random_Number_Working_Storage ()
{
  _default_values ();

  return;
}

/*
*/

Random_Number_Working_Storage::Random_Number_Working_Storage (random_number_seed_t s)
{
  _default_values ();

  set_seed (s);

  return;
}

int
Random_Number_Working_Storage::debug_print (ostream & os) const
{
  os << "Random number storage object, buffer : ";
  for (int i = 0; i < 6; i++)
    os << hex << (unsigned int) _buffer[i];

  os << dec << endl;

  return os.good ();
}

/*
  Note that this is somewhat inelegant, and could break if different
  word sizes are encountered.
  With 32 bit arithmetic, _buffer[4] should not change, so we check
  it either size of the call. If ever this becomes a problem, change
  it to do it right.
*/

void
Random_Number_Working_Storage::set_seed (random_number_seed_t s)
{
  char b4 = _buffer[4];

  random_number_seed_t * tmp = (random_number_seed_t *) & _buffer;

  *tmp = s;

  assert (_buffer[4] == b4);

  _initialised = 1;

  return;
}

random_number_seed_t
Random_Number_Working_Storage::choose_random_seed ()
{
  random_number_seed_t tt = time (0);
  tt += getpid ();

  set_seed (tt);

  return tt;
}

/*
  This is not great. Conventional wisdom says that you should
  never just test one bit of the stream from a random number
  generator, but that is what we are doing here.
*/

int
Random_Number_Working_Storage::random_one_or_zero ()
{
  long tmp = jrand48 ( (unsigned short *) & _buffer);

  return tmp > 0;
}

int
Random_Number_Working_Storage::intbtwij (int low, int high)
{
  return low + (int) (random_number () * (double) (high - low + 1));
}

random_number_t
iwrandom ()
{
  return drand48 ();
}

random_number_t
iwrandom (Random_Number_Working_Storage & rnws)
{
  return erand48 (rnws.buffer ());
}

int
intbtwij (int low, int high)
{
  return low + (int) (iwrandom() * (double) (high - low + 1));
}

int
random_one_or_zero (Random_Number_Working_Storage & rnws)
{
  long int i = nrand48 (rnws.buffer ());

  if (i < numeric_limits<long int>::max() / 2)
    return 0;

  return 1;
}

void
iw_set_rnum_seed (random_number_seed_t seed)
{
  srand48 (seed);

  return;
}

static void
jumble_bytes (char * b, int nb,
              int bfrom, int bto)
{
  char bsave = b[bto];

  b[bto] = b[bfrom];
  b[bfrom] = bsave;

  return;
}

random_number_seed_t
iw_random_seed (void)
{
  unsigned int tt = static_cast<unsigned int> (time (0));

  jumble_bytes (reinterpret_cast<char *> (&tt), sizeof(tt), 0, 3);

  int pid = static_cast<int> (getpid ());

  jumble_bytes (reinterpret_cast<char *> (&pid), sizeof(pid), 2, 3);

  random_number_seed_t seed = static_cast<random_number_seed_t>(tt | pid);

  cerr << "Using seed " << seed << endl;
  srand48(seed);

  return seed;
}
