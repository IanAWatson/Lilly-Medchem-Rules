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
#ifndef IW_SET_OF_ATOMS_H
#define IW_SET_OF_ATOMS_H

#include "iwmtypes.h"
#include "iwaray.h"

class Set_of_Atoms : public resizable_array<atom_number_t>
{
  private:
  public:
    Set_of_Atoms ();
    Set_of_Atoms (int);
    Set_of_Atoms (const Set_of_Atoms &);

    Set_of_Atoms & operator = (const Set_of_Atoms &);

    int write (ostream &) const;

    int increment_vector (int *, int = 1) const;

    int set_vector (int *, int) const;     // change to member template sometime
    int set_vector (float *, float) const;

    int any_members_set_in_array (const int * haystack, int needle) const;
    int any_members_set_in_array (const int *) const;    // looks for non-zero items

    int count_members_set_in_array (const int * haystack, int needle) const;

    int offset_atom_numbers (int);
    ostream & write_atom_numbers (ostream &, int = 0, int = 0) const;

    int adjust_for_loss_of_atom (atom_number_t, int = 0);
    int all_members_set_in_array (const int * v, int target) const;

    int all_members_non_zero_in_array (const int * v) const;
    int number_members_non_zero_in_array (const int * v) const;

    int any_members_in_common (const Set_of_Atoms &) const;
    atom_number_t first_member_in_common (const Set_of_Atoms &) const;
    int members_in_common (const Set_of_Atoms &) const;

    int write_as_mdl_v30_collection_block (const const_IWSubstring & zname, const const_IWSubstring & subname, ostream &) const;

    template <typename T> void each (Molecule &, T &) const;
};

template <typename T>
void
Set_of_Atoms::each (Molecule & m, T & o) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    o (m, _things[i]);
  }

  return;
}

ostream &
operator << (ostream &, const Set_of_Atoms &);
ostream &
operator << (ostream &, const Set_of_Atoms *);

#endif
