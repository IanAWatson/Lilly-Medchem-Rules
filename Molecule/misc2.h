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
#ifndef IW_MISC2_H
#define IW_MISC2_H

#include "iwmtypes.h"

extern void iwabort ();

extern int iw_rename (const char *, const char *);

extern int iw_getpid ();

extern int int_comparitor_larger (const int *, const int *);

extern int uint64_comparitor_smaller (const iw_uint64_t *, const iw_uint64_t *);

extern void iwxor (const int *, int *, int);

class const_IWSubstring;

extern int fetch_numeric (const const_IWSubstring & string, int & value, int max_chars = 0);

/*
  Sometimes we need to compute combinatorial permutations and we
  may be dealing with numbers larger than can be held in an int
*/

extern iw_uint64_t iw_combinatorial_combinations (int n, int k);

class iwstring_data_source;

extern int skip_to_string (iwstring_data_source & input,
                const char * target,
                int report_discard);
#endif

