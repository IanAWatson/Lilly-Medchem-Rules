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
#ifndef MDL_FUNCTIONS_H
#define MDL_FUNCTIONS_H

#include "iwstring.h"
#include "iwmtypes.h"

#include "atom_alias.h"

/*
  A couple of supporting functions used by the various mdl routines
*/

extern int int3d (const const_IWSubstring &, int &, int &, int * = NULL);

#define MDL_RADICAL -998

extern int convert_from_mdl_charge (int);

extern int convert_from_mdl_number_to_bond_type (int int_rep, bond_type_t & bt);

/*
  When reading the M lines in MDL files, we have pairs of atom numbers and atom properties.
  This class describes such a grouping.
*/

struct Aprop
{
  int _atom_number;
  int _property;
};

typedef struct Aprop Aprop;

extern int write_v30_record (IWString & buffer, ostream & output);

/*
  According to the documentation, there can be a max of MAX_PAIRS of these pairs on a record
*/

#define MAX_PAIRS 10

extern int fill_atom_property_array (const IWString & buffer, int &, Aprop * atom_properties);

extern int parse_bond_record (const_IWSubstring & buffer,
                   int na,
                   atom_number_t & a1, atom_number_t & a2,
                   int & bond_type_read_in,
                   int & directionality);

class Atom;

extern Atom * create_mdl_atom (const const_IWSubstring & ss,
                 int msdif,
                 int chg,
                 int is_radical);

void delete_digits_objects_in_mdl_file ();

#endif
