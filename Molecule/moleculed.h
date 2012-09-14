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
#ifndef COMPILING_MOLECULED
  THIS FILE SHOULD ONLY BE INCLUDED IN MOLECULED.CC
#else

    void _compute_distance_matrix ();
    int _initialise_distance_matrix ();
    int _bonds_between (atom_number_t, atom_number_t);
    int _recompute_distance_matrix (int (Molecule::*identify_first_atom) (const int *, atom_number_t &),
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &));
    void _compute_row_of_distance_matrix (int * row_of_distance_matrix,
                                 atom_number_t current_atom,
                                 int distance);
    void _compute_row_of_distance_matrix (CRDM_args & crdm,
               int (Molecule::*identify_next_atom) (const int *, atom_number_t, atom_number_t &));
    void _compute_row_of_distance_matrix (int * row_of_distance_matrix,
                                 int & distance,
                                 int * atom_stack,
                                 int stack_ptr,
                                 int * ring_atom);

    int _atoms_between (atom_number_t a1,
                          atom_number_t a2,
                          int d,
                          Set_of_Atoms & s);

#endif
