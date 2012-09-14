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
#ifndef IW_AROMATIC_H
#define IW_AROMATIC_H

#include <iostream>
using namespace std;

/*
  Aromaticity can be either by Daylight rules, or Pearlman's
*/

#define Simple_4n_plus_2 1
#define Daylight 2
#define Pearlman 3

/*
  Dec 98. The Wang Fu Lai clogp paper. They use rules which look mostly like
  Pearlman rules, but they aromatise furan and the sulphur analogue
*/

#define WangFuLai 4

#define Vijay_Gombar 5

#define EVERYTHING_HAS_A_PI_ELECTRON 6

#define Pipeline_Pilot 7

#define PUBCHEM_AROMATICITY 8

extern int set_global_aromaticity_type (int);
extern int global_aromaticity_type ();

extern int display_standard_aromaticity_options (ostream &);

class Command_Line;

extern int process_standard_aromaticity_options (Command_Line &, int = 0, char = 'A');

extern int input_aromatic_structures ();
extern void set_input_aromatic_structures (int);

extern int allow_input_without_valid_kekule_form  ();
extern void set_allow_input_without_valid_kekule_form (int);

extern int  allow_delocalised_carbonyl_bonds ();
extern void set_allow_delocalised_carbonyl_bonds (int);

extern int  discard_non_aromatic_kekule_input ();
extern void set_discard_non_aromatic_kekule_input (int);

extern int convert_chain_aromatic_bonds ();
extern void set_convert_chain_aromatic_bonds (int);

extern void set_aromatic_chain_bonds_are_ok (int s);

/*
  When outputting forms other than SMILES, are aromatic bonds written as
  aromatic, or as kekule forms
*/

extern void set_write_aromatic_bonds (int);
extern int write_aromatic_bonds ();

extern void set_warn_aromatic_chain_atoms (int);

extern void set_kekule_try_positive_nitrogen (int s);

extern void set_all_bonds_in_aromatic_ring_must_be_aromatic (int s);

extern void set_display_no_kekule_form_message(int s);
extern int  display_no_kekule_form_message();

extern void set_allow_pipeline_pilot_aromaticity_on_input(int s);

/*
  The largest ring that can be aromatic
*/

extern void set_max_aromatic_ring_size (int s);

extern void set_perform_kekule_perception (int s);

extern void reset_aromatic_file_scope_variables ();
extern void  reset_mdl_file_scope_variables ();

#endif

/* arch-tag: cd7a479d-ca70-4047-b1d1-1d7beb0bf194 */
