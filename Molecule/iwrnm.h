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
#ifndef IW_RING_NUMBER_MANAGER_H
#define IW_RING_NUMBER_MANAGER_H

/*
  This class is used when constructing a smiles.
  When placing a ring opening digit, we must keep track of
    (a) The atom which completes the ring number
    (b) The type of the bond
    (c) The atom which created the ring number
*/

#include "chiral_centre.h"

class Ring_Number_Manager
{
  private:
    int _nr;
    int * _ring_id;
//  bond_type_t * _bt;
    const Bond ** _bond;
    atom_number_t * _from_atom;

    int _include_aromaticity_in_smiles;    // initialised during constructor and never changed
    int _include_cis_trans_in_smiles;      // initialised during constructor and never changed

//  private functions

    void _append_ring_closure_digits (IWString & smiles,
                            int ring_closure_number,
                            const Bond * b,
                            atom_number_t ato) const;

    int _process_ring (IWString & smiles, int ring, atom_number_t afrom);

    int _place_ring_closure (IWString & smiles,
                   atom_number_t a,
                   atom_number_t afrom);
    int _append_ring_closures_for_chiral_atom (IWString & smiles,
                    atom_number_t a,
                    const resizable_array<atom_number_t> & ring_closures_found);

    void _default_values ();

  public:
    Ring_Number_Manager ();
    Ring_Number_Manager (int);
    ~Ring_Number_Manager ();

    int debug_print (ostream &) const;

    int ok () const;

    int activate (int);

    int store_ring (IWString &, const Bond *, atom_number_t);

    int append_ring_closures_for_atom (IWString &,
                atom_number_t, 
                const resizable_array<atom_number_t> &,
                const Chiral_Centre *);
};

#endif
