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
#ifndef MOLECULE_MAIN_H
#define MOLECULE_MAIN_H

  private:
    int _convert_set_of_atoms_to_bond_numbers (const Set_of_Atoms & s, int * barray) const;

    int  _ok_ring_info () const;
    int  _invalidate_ring_info ();
    int  _invalidate_ring_aromaticity_info ();

    void _compute_element_count (int * element_count, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;
    void _compute_element_count (int * element_count, const int * include_atom, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;
    void _compute_element_count (int * element_count, const int * atom_flag, int flag, int & highest_atomic_number, int & isotopes_present, int & non_periodic_table_elements_present) const;

    int  _remove_atom (atom_number_t);

    int _invalidate_for_changed_isotope ();
    int _exact_mass (const int * element_count, int highest_atomic_number,
                     int non_periodic_table_elements_present,
                     exact_mass_t & result) const;

    int _set_bond_length (atom_number_t a1, atom_number_t a2,
                            distance_t d, int * either_side);

    int _set_isotope_zero(atom_number_t zatom);

#endif
