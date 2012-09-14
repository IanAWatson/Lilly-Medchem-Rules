/**************************************************************************

    Copyright (C) 2012  Eli Lilly and Company

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
#ifndef IW_SMILES_H
#define IW_SMILES_H 1

#include "molecule.h"
#include "iwrandom.h"

/*
  This contains functions used internally by the smiles routines.
*/


extern boolean include_atom_in_smiles (const Molecule *, atom_number_t);

/*
  Usually, aromaticity is never written. 
  Aromaticity can be written if these functions are used.
*/

extern void set_include_aromaticity_in_smiles (int);    // smiles only
extern int  get_include_aromaticity_in_smiles ();       // smiles only

extern void set_include_cis_trans_in_smiles (int);
extern int  include_cis_trans_in_smiles ();

extern void set_include_chiral_info_in_smiles (int);
extern int  include_chiral_info_in_smiles ();

/*
  Sets include_chiral_info_in_smiles to the constructor's value and
  restores initial value on destruction
*/

class Temporarily_Set_Include_Chiral_Info_in_Smiles
{
  private:
    int _initial_value;

  public:
    Temporarily_Set_Include_Chiral_Info_in_Smiles(int);
    ~Temporarily_Set_Include_Chiral_Info_in_Smiles();
};

extern void set_smiles_random_number_seed (random_number_seed_t);
extern random_number_seed_t set_smiles_random_number_seed_random ();

extern int smiles_process_atom (Molecule *, IWString &, atom_number_t, bond_type_t,
                                 atom_number_t, chiral_type_t = NON_CHIRAL);

class Command_Line;

extern int display_standard_smiles_options (ostream &);

extern int process_standard_smiles_options (Command_Line &, int = 0, const char = 'K');

extern void set_smiles_reuse_ring_closure_numbers (int);
extern int  smiles_reuse_ring_closure_numbers ();

extern void set_append_coordinates_after_each_atom (int);

extern void set_smiles_native_ordering (int);

extern int set_datatype_name_for_structure_in_tdt_files (const char *);

extern void set_tdt_append_dataitem_content (int s);

extern int smiles_error_message (const char *, int, int, const char *);

extern int set_default_unique_smiles_aromaticity (int a);

extern void set_unset_implicit_hydrogens_known_if_possible (int s);
extern void set_unset_all_implicit_hydrogens_known_attributes (int s);

extern void set_include_isotopic_information_in_unique_smiles (int s);
extern int  include_isotopic_information_in_unique_smiles ();

extern void set_include_directional_bonding_information_in_unique_smiles (int s);

extern void set_include_implicit_hydrogens_on_aromatic_n_and_p (int);

extern void set_smiles_ring_number_offset (int s);

extern void set_display_smiles_interpretation_error_messages(int s);

extern int display_smiles_interpretation_error_messages();

#define CANON_IMPH_CONSIDER_ALL 86
#define CANON_IMPH_CONSIDER_JUST_HETEROATOMS 73
#define CANON_IMPH_CONSIDER_NONE 61

extern void set_consider_implicit_hydrogens_in_unique_smiles(int s);
extern void set_consider_implicit_hydrogens_known_in_unique_smiles(int s);

// Parse_smiles_token is external because it is also used by the smarts routines

extern int
parse_smiles_token (const char * smiles,
                    int characters_to_process,
                    const Element * &    e,
                    aromaticity_type_t & aromatic,
                    formal_charge_t &    fc,
                    int             &    hcount,
                    int             &    chiral_count,
                    int             &    atomic_mass);
int
parse_smiles_token (const char * smiles,
                    int characters_to_process,
                    const Element * &    e,
                    int & aromatic);

/*
  To avoid passing around a lot of arguments when generating a smiles, we bundle
  all the information into an object
*/

#include "iwrnm.h"

class Smiles_Formation_Info
{
  private:
    int _natoms;

    int * _already_done;

    Ring_Number_Manager _rnm;

    atom_number_t _previous_atom;

    atom_number_t _zatom;

    const int * _include_atom;

    int _write_smiles;    // are we writing smiles or smarts

//  The per-atom create embedding information comes from
//  the Smiles_Information object

    const int * _make_smarts_embedding;

//  Similarly, any user specified per-atom smarts information comes from
//  the Smiles_Information object

    const IWString * _user_specified_atomic_smarts;

  public:
    Smiles_Formation_Info (int na, int nr);
    ~Smiles_Formation_Info ();

    void set_already_done (int * s) { _already_done = s;}
    int * already_done () { return _already_done;}

    Ring_Number_Manager & rnm () { return _rnm;}

    int ok () const;

    void set_zatom (atom_number_t s) { _previous_atom = _zatom; _zatom = s;}
    void set_previous_atom (atom_number_t s) { _previous_atom = s;}

    atom_number_t zatom () { return _zatom;}
    atom_number_t previous_atom () { return _previous_atom;}

    const int * include_atom () const { return _include_atom;}

    void set_include_atom (const int * s) { _include_atom = s;}

    void set_write_smiles (int s) { _write_smiles = s;}
    int  write_smiles () const { return _write_smiles;}

    void set_make_smarts_embedding (const int * s) { _make_smarts_embedding = s;}
    int make_smarts_embedding (atom_number_t) const;

    void set_user_specified_atomic_smarts(const IWString * const s) { _user_specified_atomic_smarts = s;}
    const IWString * user_specified_atomic_smarts() const { return _user_specified_atomic_smarts;}
};

extern void set_write_smiles_with_smarts_atoms (int);
extern int  write_smiles_with_smarts_atoms ();
extern void set_max_ring_digits_following_percent (int);

extern void set_resolve_unique_smiles_ties_by_geometry(int);

extern void set_consider_isotopes_as_zero_and_non_zero (int s);
extern int  consider_isotopes_as_zero_and_non_zero (int s);

extern int  add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens();
extern void set_add_implicit_hydrogens_to_isotopic_atoms_needing_hydrogens(int);

extern void set_write_single_bonds_in_smiles(int);
extern int  write_single_bonds_in_smiles ();

extern void set_include_hcount_in_smiles(int);

extern void set_write_smiles_aromatic_bonds_as_colons(int);
extern int  write_smiles_aromatic_bonds_as_colons ();

extern  void set_display_unusual_hcount_warning_messages(int);
extern  int  display_unusual_hcount_warning_messages ();

extern  void set_include_directionality_in_ring_closure_bonds (int);

extern void reset_smi_file_scope_variables ();
extern void reset_smiles_support_file_scope_variables ();

/*
  Used for parsing leading numeric qualifiers
*/

extern int smarts_fetch_numeric (const char * string, int & value, int & qualifier);

#endif

/* arch-tag: ebd60b09-90c5-49a4-9862-1fbafd6b939c

*/
