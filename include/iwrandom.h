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
#ifndef IWRANDOM_H
#define IWRANDOM_H

#include <stdlib.h>
#include <iostream>

using namespace std;

typedef unsigned long int random_number_seed;
typedef unsigned long int random_number_seed_t;
typedef double random_number;
typedef double random_number_t;

/*
  On an SGI, the random number stuff is done in 48 bit
  arithmetic
*/

class Random_Number_Working_Storage
{
  private:
    char _buffer[6];
    int  _initialised;

// private functions
  
  void _default_values ();

  public:
    Random_Number_Working_Storage ();
    Random_Number_Working_Storage (random_number_seed_t);

    int initialised () const { return _initialised;}

    random_number_seed_t choose_random_seed ();

    int debug_print (ostream & = cerr) const;

    void set_seed (random_number_seed_t);

    unsigned short * buffer () const { return (unsigned short *) & _buffer;}

    int random_one_or_zero ();
    random_number_t random_number () { return ::erand48 ((unsigned short *) &_buffer);}
    int intbtwij (int, int);
};

extern random_number iwrandom (void);
extern int intbtwij (int, int);
extern void iw_set_rnum_seed (random_number_seed_t);
extern random_number_seed_t iw_random_seed ();

template <typename T> T random_number_between (const T & low, const T & high);

#ifdef RANDOM_NUMBER_BETWEEN_IMPLEMENTATION

template <typename T>
T
random_number_between (const T & low, const T & high)
{
  double range = static_cast<double> (high - low + 1);

  T delta = static_cast<T> (range * iwrandom ());

  return low + delta;
}

#endif

#endif
