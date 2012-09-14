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
#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "iwaray.h"

#include "bond.h"

class Beep;

#define BOND_LIST_MAGIC -134211

class Bond_list : public resizable_array_p<Bond>
{
  private:
    magic_number_t _magic;

// private functions

    int _maximum_connectivity (int *, int) const;

  public:
    Bond_list ();
    ~Bond_list ();

    int ok () const;

    int debug_print (ostream &) const;

    int nbonds () const { return _number_elements;}

    int     which_bond (atom_number_t, atom_number_t) const;

    int     remove_bonds_to_atom (atom_number_t, int = 0);
    int     remove_bond_between_atoms (atom_number_t, atom_number_t);

    Bond *  bond_between_atoms (atom_number_t, atom_number_t) const;

    int     swap_atoms   (int, int);
    int     move_atom_to_end_of_atom_list (atom_number_t, int);

//  In case someone wants rapid access to the bond types

    int     copy_bond_types (bond_type_t *) const;

    int     set_modified ();

#ifdef BONDS_KNOW_RING_MEMBERSHIP
    int     invalidate_ring_info ();
    int     assign_ring_membership_to_bonds (const resizable_array_p<Beep> & beeps);
#endif

    int     assign_bond_numbers (int istart);
    int     assign_bond_numbers_to_bonds_if_needed();
    void    invalidate_bond_numbers ();

    int     set_all_bond_types (bond_type_t);

    int     unset_all_permanent_aromatic_bonds ();

    int     cis_trans_bonds_present() const;
};

#define OK_BOND_LIST(b) ( NULL != (b) && b->ok () )

#endif
