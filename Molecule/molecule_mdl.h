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
#ifndef COMPILING_MDL_CC
  THIS FILE SHOULD ONLY BE INCLUDED WHEN COMPILING MDL.CC
#endif

  private:

    int  _read_molecule_mdl_ds (iwstring_data_source &, int);

//  int  _read_possibly_aromatic_molecule_mdl_ds (iwstring_data_source &, int, int);
//  int  _read_possibly_aromatic_molecule_mdl_ds (iwstring_data_source &, int);

    int  _read_mdl_atom_connection_table (iwstring_data_source &, int &, int &);
    int  _read_mdl_bond_list (iwstring_data_source &, int, int *, int *);
//  int  _read_possibly_aromatic_mdl_bond_list (iwstring_data_source &, int, int *);
    int  _read_molecule_mdl_trailing_records (iwstring_data_source &, int);

    int _write_M_RGP_records (ostream &) const;
  protected:
    int  _write_m_chg_records (ostream & os, int) const;
    int  _write_m_iso_records (ostream & os, int) const;
    int _common_parse_M_record (const const_IWSubstring & buffer,
                                int & fatal);

    int _mdl_atom_is_chiral_centre (atom_number_t zatom, int cfg);

  private:

    int _has_delocalised_neighbours (atom_number_t zatom,
                                       const int * aromatic_atoms,
                                       const int * aromatic_bonds,
                                       Set_of_Atoms & s) const;
    int _unset_aromatic_bond_entry(atom_number_t a1, 
                                     atom_number_t a2,
                                     int * aromatic_bonds) const;

    int _more_than_one_aromatic_bond(atom_number_t zatom,
                                       const int * aromatic_bond) const;

    int _read_molecule_rdf_ds (iwstring_data_source & input, IWString & possible_name);

//  May 98, stuff for Version 3 sd files

    int _read_mdl_V3 (iwstring_data_source &);
    int _write_mdl_V3 (iwstring_data_source &);
    int _read_mdl_atom_connection_table_v30 (iwstring_data_source & input, int & nb);
  protected:
    int _parse_v30_atom_record (const IWString & buffer, int = 1);
  private:
    int _read_v30_bond_list (iwstring_data_source & input, int nb, int *, int *);

    int _write_molecule_atom_list_v30 (ostream & os) const;
    int _write_molecule_bond_list_v30 (ostream & os) const;

//  Sept 2000, stuff for discerning chirality from wedge bonds

    int _discern_chirality_from_wedge_bond (atom_number_t a1, atom_number_t a2,
                          int direction);
    int _discern_chirality_from_wedge_bond_4 (atom_number_t zatom, atom_number_t a2,
                          int direction);
    int _create_chiral_centre (atom_number_t zatom,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 atom_number_t a3,
                                 atom_number_t a4,
                                 int direction);
    int _create_unspecified_chirality_object (atom_number_t zatom);

    int _mdl_atom_stereo_value (atom_number_t a) const;

    int _multiple_wedge_bonds_to_atom (atom_number_t a) const;
    int _discern_chirality_from_multiple_wedge_bonds (int bstart);

    int _parse_M_RGP_record (const const_IWSubstring & buffer);

    int _read_mdl_data_following_tag (iwstring_data_source & input);

    int _set_elements_based_on_atom_aliases (const resizable_array_p<Atom_Alias> & a);

  public:
    int  mdl_add_m_formal_charge  (int, const Aprop * atom_properties);
    int  mdl_add_m_radical (int, const Aprop * atom_properties);
    int  mdl_add_m_isotope (int, const Aprop * atom_properties);

//  Shared with the function that writes ISIS reaction files

    int _write_mdl_atom_record_element_charge_and_chirality (atom_number_t, IWString & output_buffer) const;
    int _mdl_write_atoms_and_bonds_record (ostream &, int nfc, int iat, int isis_standard_records) const;

    int _process_mdl_g_record (const IWString &, const const_IWSubstring & buffer);
    int _assign_strange_atomic_symbol_to_atom (atom_number_t zatom, const_IWSubstring s);

  private:
