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
#ifndef IW_RING_CLOSURE_BONDS_H
#define IW_RING_CLOSURE_BONDS_H

#include "bond.h"
#include "iwaray.h"

/*
  We need to keep track of whether or not given pairs of atoms are
  stored. Some trickery with the _present array does that if both
  atoms are hit just once. Otherwise, we store integers
  (a1 * _atoms_in_molecule + a2) in the resizable_array<int>
*/

class Ring_Closure_Bonds : public resizable_array<int>
{
  private:
    int _atoms_in_molecule;

    int * _present;

//  private functions

    int _form_corresponding_integer (atom_number_t a1, atom_number_t a2) const;

  public:
    Ring_Closure_Bonds ();
    Ring_Closure_Bonds (const Ring_Closure_Bonds &);
    ~Ring_Closure_Bonds ();

    Ring_Closure_Bonds & operator= (const Ring_Closure_Bonds &);

    int ok () const;

    int write_bonds (ostream & output) const;

    int reset ();

    void invalidate ();

    int activate (int);

    int add (atom_number_t, atom_number_t);
    int contains (atom_number_t, atom_number_t) const;

//  I didn't call this operator== because it the order of the integers in the resizable_array<int> may be different

    int is_the_same (const Ring_Closure_Bonds &) const;

    int report_differences (const Ring_Closure_Bonds &, ostream &) const;

    int is_subset_of (const Ring_Closure_Bonds &) const;
};

#endif
