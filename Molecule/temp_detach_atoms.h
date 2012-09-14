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
#ifndef TMP_DETACH_ATOMS
#define TMP_DETACH_ATOMS

/*
  When dealing with explicit hydrogens, both the donor_acceptor and
  charge_assigner objects may need to temporarily remove explicit
  Hydrogens in order for the queries to work
*/

#include "molecule.h"

class Temp_Detach_Atoms
{
  private:

    int _active;

    int _verbose;

// When detaching-reattaching explicit Hydrogens, by default we remove
// any Hydrogen no longer needed (on a Carboxyllic acid for example). 
// We ran into cases where we needed to preserve the number of atoms
// in the molecule. We break the bond, but just leave the atom.

    int _remove_hydrogens_no_longer_needed;

//  As a small consistency check, we record the number of atoms in
//  the molecule we most recently broke apart. Woe to anyone who
//  passes different molecules to different calles to detach_atoms and reattach_atoms!!

    int _matoms;

//  if we didn't detach any atoms, we don't need to reattach any

    int _need_to_reattach;

    int * _connection;

    bond_type_t _bt;

  public:
    Temp_Detach_Atoms ();
    ~Temp_Detach_Atoms ();

    void set_verbose (int v) { _verbose = v;}

//  The object recognised directives like "nodetach" and "noremove"

    int recognised_directive (const const_IWSubstring &);

//  We may no longer want a specific atom re-attached

    void do_not_reattach_to_atom (atom_number_t);

    int active () const { return _active;}
    int natoms () const { return _matoms;}

    int detach_atoms (Molecule &, atomic_number_t z = 1);
    int reattach_atoms (Molecule &);
};

#endif
