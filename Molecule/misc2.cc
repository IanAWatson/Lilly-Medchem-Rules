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

#ifdef UNIX
#include <unistd.h>
#include <stdio.h>
#endif

#include "misc2.h"
#include "iwconfig.h"

/*
  This file contains miscelaneous functions which cannot go in misc1.cc
  because of misc1 instantiates templates.
*/

#include "iwstring_data_source.h"
int
iw_rename (const char * old_name, const char * new_name)
{
  assert (NULL != old_name);
  assert (NULL != new_name);

  assert (strlen (old_name));
  assert (strlen (new_name));

  return rename (old_name, new_name);
}

int
iw_getpid ()
{
  return IW_GETPID ();
}

/*
  An int comparitor for use with qsort.
*/

int
int_comparitor_larger (const int * pi1, const int * pi2)
{
  if (*pi1 < *pi2)
    return -1;
  else if (*pi1 == *pi2)
    return 0;
  else   
    return 1;
}

int
uint64_comparitor_smaller (const iw_uint64_t * pi1, const iw_uint64_t * pi2)
{
  if (*pi1 < *pi2)
    return 1;
  else if (*pi1 == *pi2)
    return 0;
  else   
    return -1;
}

void
iwxor (const int * i1,
       int * i2,
       int n)
{
  for (int i = 0; i < n; i++)
  {
    i2[i] = (i2[i] ^ i1[i]);
  }

  return;
}

#define DEBUG_SHORTENING

/*
  This function is used in the ring finding stuff based on bonds.

  We have two sets of bonds. See if the first, I1, can be XOR'd with
  I2 to shorten I2. If so, do the XOR, and return 1
*/



/*
  A frequent operation when parsing smiles/smarts is to extract one or
  more digits.  We return the number of characters we fully parse.
*/

int
fetch_numeric (const const_IWSubstring & zstring, int & result, int max_chars)
{
  int rc = 0;
  result = 0;

  int characters_to_search = zstring.length ();
  if (max_chars > 0 && max_chars < characters_to_search)
    characters_to_search = max_chars;

  for (int i = 0; i < characters_to_search; i++)
  {
    int tmp = zstring[i] - '0';
    if (tmp < 0 || tmp > 9)
      return rc;

    result = result * 10 + tmp;
    rc++;
  }

  return rc;
}

void
iwabort ()
{
  abort ();
}


/*
  The structure reading functions need the ability to skip over the rest
  of a bad connection table
*/

int
skip_to_string (iwstring_data_source & input,
                const char * target,
                int report_discard)
{
  if (input.at_eof ())
    return 0;

  int records_discarded = 0;
  int start_record = input.lines_read ();

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if (buffer.starts_with (target))
    {
      if (records_discarded && report_discard)
        cerr << "skip to string: '" << target << "' discarded " <<
                 records_discarded << " records starting at " << start_record << endl;
      return 1;
    }

    records_discarded++;
  }

  if (report_discard)
    cerr << "skip to string: reached EOF, no '" << target << "' found\n";

  return 0;
}
