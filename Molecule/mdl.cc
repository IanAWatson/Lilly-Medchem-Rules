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
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <fstream>
#include <iomanip>
#include <memory>
using namespace std;

// Be sure to define this symbol so all the private functions get defined

#define COMPILING_MDL_CC
#define COMPILING_CTB

#include "misc.h"
#include "iwstring_data_source.h"
#include "iwcrex.h"

#include "mdl.h"
#include "atom_alias.h"
#include "molecule.h"
#include "misc2.h"
#include "chiral_centre.h"
#include "rwmolecule.h"
#include "aromatic.h"
#include "mdl_atom_record.h"

//#define USE_IWMALLOC
#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

/*
  Some people wanted as close to the ISIS standard as possible
*/

static int isis_standard_records = 0;

void 
set_write_isis_standard (int i)
{
  isis_standard_records = i;
}

static int ignore_unrecognised_m_records = 0;

void
set_ignore_unrecognised_mdl_m_records (int i)
{
  ignore_unrecognised_m_records = i;
}

static int report_unrecognised_records = 1;

void
set_mdl_report_unrecognised_records (int i)
{
  report_unrecognised_records = i;
}

static int die_on_erroneous_m_input = 1;

void 
set_die_on_erroneous_m_input (int s)
{
  die_on_erroneous_m_input = s;
}

static int write_mdl_dollar_signs = 1;

void
set_write_mdl_dollars (int s)
{
   write_mdl_dollar_signs = s;
}

/*
  If set to a value > 1, that indicates write it even if the current
  molecule does not really need it
*/

static int write_mdl_m_end_record = 1;

void
set_write_mdl_m_end_record (int w)
{
  write_mdl_m_end_record = w;
}

static int write_v30_mdl_files = 0;

void
set_write_v30_mdl_files (int s)
{
  write_v30_mdl_files = s;
}

static int _ignore_self_bonds = 0;

void
set_mdl_ignore_self_bonds (int s)
{
  _ignore_self_bonds = s;
}

int
ignore_self_bonds()
{
  return _ignore_self_bonds;
}

static int _write_mdl_chiral_flags = 1;

void
set_write_mdl_chiral_flags (int s)
{
  _write_mdl_chiral_flags = s;
}

int
write_mdl_chiral_flags()
{
  return _write_mdl_chiral_flags;
}

static int _include_chiral_info_in_mdl_outputs = 1;

void
set_include_chiral_info_in_mdl_outputs (int s)
{
  _include_chiral_info_in_mdl_outputs = s;
}

int
include_chiral_info_in_mdl_outputs()
{
  return _include_chiral_info_in_mdl_outputs;
}

/*
  Jun 2001. Discovered bug in how we read chirality when explicit Hydrogens
  are present, and when the explicit Hydrogen isn't the lowest numbered atom
  in at the chiral centre
*/

static int _mdl_read_h_correct_chiral_centres = 1;

void
set_mdl_read_h_correct_chiral_centres (int s)
{
  _mdl_read_h_correct_chiral_centres = s;
}

int
mdl_read_h_correct_chiral_centres()
{
  return _mdl_read_h_correct_chiral_centres;
}

static int _mdl_write_h_correct_chiral_centres = 1;

void 
set_mdl_write_h_correct_chiral_centres (int s)
{
  _mdl_write_h_correct_chiral_centres = s;
}

static IWString insert_between_sdf_name_tokens(' ');

void
set_mdl_insert_between_sdf_name_tokens (const const_IWSubstring & s)
{
  insert_between_sdf_name_tokens = s;
}

/*
  A common need is to extract the identifier from a > <> record in an
  sdf file.

  We go to a little trouble to munge the user's input so that it matches
  ^>  <\S+>$

  Oct 2009. Start seeing things of the form

  > 10014  <MODEL.EXTREG>

*/

static IW_Regular_Expression sdf_identifier;

int
set_sdf_identifier (const const_IWSubstring & sdfid)
{
  IWString mysdfid = sdfid;

  if (mysdfid.starts_with('^'))
    mysdfid.remove_item(0);

  IWString tmp = "^>.*<";
  tmp += mysdfid;

  mysdfid = tmp;

  if (mysdfid.ends_with('$'))
    mysdfid.chop();

  mysdfid += ">";

  if (! sdf_identifier.set_pattern(mysdfid))
  {
    cerr << "Cannot set sdfid pattern to '" << mysdfid << "'\n";
    return 0;
  }

//cerr << "Pattern set to '" << mysdfid << "'\n";

  return 1;
}

static int extract_isis_extregno = 0;

void 
set_extract_isis_extregno (int e)
{
  extract_isis_extregno = e;
}

/*
  This is blank padded because we use a contains operation.
  Don't want to have to deal with variable number of spaces
  after the M
*/

static IWString name_in_m_tag;

void
set_mdl_name_in_m_tag(const const_IWSubstring & s)
{
  if (name_in_m_tag.length() > 0)
    name_in_m_tag.resize_keep_storage(0);

  name_in_m_tag << ' ' << s << ' ';
}

static int fetch_all_sdf_identifiers = 0;

void 
set_fetch_all_sdf_identifiers (int f)
{
  fetch_all_sdf_identifiers = f;
}

static int take_first_tag_as_name = 0;

void
set_take_first_tag_as_name(int s)
{
  take_first_tag_as_name = s;
}

/*
  by default, we prepend the SDF identifier to the data value and follow with a
  colon
*/

static int prepend_sdfid = 1;

void
set_prepend_sdf_identifier (int s)
{
  prepend_sdfid = s;
}

/*
  Jun 99. Found a case where we wanted to throw away what was in the name
  field of a molecule - the identifiers are in the fields
*/

static int discard_sdf_molecule_name = 0;

void
set_discard_sdf_molecule_name (int s)
{
  discard_sdf_molecule_name = s;
}

/*
  The convention is that following a tag
  >  <TAG>

  there is an arbitrary number of records of "data" followed by a blank line.
  The blank line may or may not be present. 
*/

static int _multi_record_tag_data_present = 1;

void
set_multi_record_tag_data_present (int s)
{
  _multi_record_tag_data_present = s;
}

int
multi_record_tag_data_present()
{
  return _multi_record_tag_data_present;
}

/*
  For Michal Vieth, we can optionally write aromatic bonds
*/

static int mdl_write_aromatic_bonds = 0;

void
set_mdl_write_aromatic_bonds (int s)
{
  mdl_write_aromatic_bonds = s;
}

static int mdl_write_aromatic_atoms = 0;

void
set_mdl_write_aromatic_atoms (int s)
{
  mdl_write_aromatic_atoms = s;
}

/*
  Feb 99.
  For the GER project we need to keep track of atoms which are marked
  as unspecified chiral centres - normally these are discarded.
  I don't want to make this a permanent property of all molecules, so
  we use this horrible hack.
  We also accumulate up/down and unspecified bonds
*/

static int accumulate_mdl_chirality_features = 0;

static Set_of_Atoms unspecified_chiral_atoms_last_molecule;

static int display_non_organic_chirality_messages = 1;

void
set_display_non_organic_chirality_messages (int s)
{
  display_non_organic_chirality_messages = s;
}

static int mdl_display_invalid_chiral_connectivity = 1;

void
set_mdl_display_invalid_chiral_connectivity (int s)
{
  mdl_display_invalid_chiral_connectivity = s;
}

/*
  The mdl documentation says that the first atom is the
  chiral entity in a directional bond, so we just record
  the first atom in each directional bond
*/

static Set_of_Atoms up_bonds_last_molecule, down_bonds_last_molecule, squiggle_bonds_last_molecule;

/*
  Directional flag 3 with a double bond means unspecified cis/trans
*/

static Set_of_Atoms unspecified_double_bond_atoms;

void
set_mdl_accumulate_mdl_chirality_features (int s)
{
  accumulate_mdl_chirality_features = s;
}

const Set_of_Atoms &
mdl_unspecified_chiral_atoms()
{
  return unspecified_chiral_atoms_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_up_bonds()
{
  return up_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_squiggle_bonds()
{
  return squiggle_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_atoms_with_down_bonds()
{
  return down_bonds_last_molecule;
}

const Set_of_Atoms &
mdl_unspecified_double_bond_atoms()
{
  return unspecified_double_bond_atoms;
}

static int truncate_long_symbols = 0;

static IWString change_long_symbols_to;

void
set_mdl_truncate_long_elements (int s)
{
  truncate_long_symbols = s;

  return;
}

void
set_mdl_change_long_symbols_to (const const_IWSubstring & s)
{
  assert(s.length() > 0 && s.length() <= 2);

  change_long_symbols_to = s;

  const Element * e = get_element_from_symbol_no_case_conversion(s);

  if (NULL == e)
    e = create_element_with_symbol(s);

  return;
}

static int discern_chirality_from_wedge_bonds = 0;

void
set_mdl_discern_chirality_from_wedge_bonds (int s)
{
  discern_chirality_from_wedge_bonds = s;
}

int
mdl_discern_chirality_from_wedge_bonds()
{
  return discern_chirality_from_wedge_bonds;
}

/*
  Long-time go-around with MDL on how isotopes should be represented.
  Their documentation lies
*/

static int write_isotopes_as_numbers_rather_than_differences_from_normal = 0;

void
set_mdl_write_isotopes_as_numbers_rather_than_differences_from_normal (int s)
{
  write_isotopes_as_numbers_rather_than_differences_from_normal = s;
}

static int write_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;

void
set_mdl_write_M_isotopes_as_numbers_rather_than_differences_from_normal (int s)
{
  write_M_isotopes_as_numbers_rather_than_differences_from_normal = s;
}

static int read_isotopes_as_numbers_rather_than_differences_from_normal = 0;

void
set_mdl_read_isotopes_as_numbers_rather_than_differences_from_normal (int s)
{
  read_isotopes_as_numbers_rather_than_differences_from_normal = s;
}

static int read_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;

void
set_mdl_read_M_isotopes_as_numbers_rather_than_differences_from_normal (int s)
{
  read_M_isotopes_as_numbers_rather_than_differences_from_normal = s;
}

/*
  June 2012. We have datasets that come in with various sdf tags, but we want to
  be able to replace the sdf tag with a string
*/

static IWString replace_first_sdf_tag;

void
set_replace_first_sdf_tag (const IWString & s)
{
  replace_first_sdf_tag = s;
}

/*
  This parses a string into int's, with each field being three digits long

  This would be trivial with FORTRAN, what am I missing?
*/

int
int3d (const const_IWSubstring & buffer, int & i1, int & i2, int * i3)
{
  if(buffer.length() < 6)
    return 0;

  i1 = 0;
  int tmp = -1;     // this is outside the loop to make sure that we detect blank fields.
  for (int i = 0; i < 3; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 0;

    i1 = i1 * 10 + tmp;
  }

  if (-1 == tmp)
    return 0;

  tmp = -1;

  i2 = 0;
  for (int i = 3; i < 6; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 1;

    i2 = i2 * 10 + tmp;
  }

  if (-1 == tmp)
    return 1;

  if (NULL == i3)
    return 2;

  tmp = -1;

  *i3 = 0;
  for (int i = 6; i < 9; i++)
  {
    if (' ' == buffer[i])
      continue;

    tmp = buffer[i] - '0';
    if (tmp < 0 || tmp > 9)
      return 1;

    *i3 = *i3 * 10 + tmp;
  }

  if (-1 == tmp)
    return 2;

  return 3;
}

static int allow_deuterium = 0;
static int allow_tritium = 0;

void
set_mdl_allow_deuterium (int d)
{
  allow_deuterium = d;
}

int
interpret_d_as_deuterium ()
{
  return allow_deuterium;
}

void
set_mdl_allow_tritium (int t)
{
  allow_tritium = t;
}

int
interpret_t_as_tritium ()
{
  return allow_tritium;
}

static int * input_bond_type_translation_table = NULL;

int
set_mdl_input_bond_type_translation (int zfrom, int zto)
{
  if (NULL == input_bond_type_translation_table)
  {
    input_bond_type_translation_table = new int[20];
    for (int i = 0; i < 20; i++)
    {
      input_bond_type_translation_table[i] = i;
    }
  }

  assert (zfrom >= 0 && zfrom < 20);

  input_bond_type_translation_table[zfrom] = zto;

  return 1;
}

/*
  This function converts from MDL charge codes, to charge
  codes used by this programme.

  Note that we have a special number for radicals
*/

int
convert_from_mdl_charge (int chg)
{
  if (0 == chg)     // or perhaps it really means -4!!! curse mdl
    return 0;

  if (4 == chg)
    return MDL_RADICAL;

  return 4 - chg;
}

int
convert_to_mdl_charge (int chg)
{
  switch(chg)
  {
    case 0:
      return 0;
    case 1:
      return 3;
    case 2:
      return 2;
    case 3:
      return 1;
    case 4:
      return 0;
    case -1: 
      return 5;
    case -2:
      return 6;
    case -3:
      return 7;
  }

  cerr << "convert to mdl_charge: bad charge " << chg << "\n";

  return 0;
}

static int mdl_g_records_hold_atom_symbols = 0;

void
set_mdl_g_records_hold_atom_symbols (int s)
{
  mdl_g_records_hold_atom_symbols = s;
}

/*
  A function which fills the array based on a buffer
  A typical buffer might look like

  A, A3, I4, I4, I4,

  M  CHG  4   7  -1   8  -1  11   1  12   1
  12345678901234567890123456789012345678901234567890123456789012345678901234567890
           1         2         3         4         5         6         7         8

  BEWARE, THEY LIE!!

  I've seen several cases where the count is incorrect (more tokens than
  are present). Guard against this!

  May 2002. the desktop version produces a different number of spaces than the host
  version. Change to just extracting tokens - this will presumably break if atom
  numbers go above 999. 
*/

int
fill_atom_property_array (const IWString & buffer,
                          int & npairs,
                          Aprop * atom_properties)
{
  assert (buffer.starts_with('M'));

  int nw = buffer.nwords();

  npairs = nw - 3;      // discount M  XXX c tokens

  if (npairs != npairs / 2 * 2)    // must be an even number of tokens
  {
    cerr << "mdl:must have even pairs '" << buffer << "'\n";

    if(die_on_erroneous_m_input)
      return 0;

    npairs--;    // convert to even number

    if (npairs <= 0)    // ignoring now empty record
      return 1;
  }

  npairs = npairs / 2;

  if (0 == npairs)
  {
    cerr << "mdl:very strange, M record with no info '" << buffer << "'\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);    // M
  buffer.nextword(token, i);    // XXX
  buffer.nextword(token, i);    // count

  for (int j = 0; j < npairs; j++)
  {
    buffer.nextword(token, i);
    
    if (! token.numeric_value(atom_properties[j]._atom_number))
    {
      cerr << "mdl:invalid atom number '" << buffer << "'\n";
      return 0;
    }

    buffer.nextword(token, i);

    if (! token.numeric_value(atom_properties[j]._property))
    {
      cerr << "mdl: invalid numeric value '" << buffer << "'\n";
      return 0;
    }
  }
  
  return npairs;
}

/*
  Remember that the numbers in M  ISO records are differences from the
  normal mass
*/

int
Molecule::mdl_add_m_isotope (int ntokens,
                             const Aprop * atom_properties)
{
  assert (ntokens > 0);

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    int k = atom_properties[i]._property;

//  cerr << "Isotope for atom " << j << " value " << k << "\n";
    if (0 == j && 0 == k)     // have seen this a couple of times. Hard to know what they intended
      continue;

    if (0 == k)     // have seen this one too!
      continue;

    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_isotope: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }

    Atom * a = _things[j];             // the j'th atom of this molecule

    if(read_M_isotopes_as_numbers_rather_than_differences_from_normal)
      a->set_isotope(k);
    else
      a->set_isotope(a->element()->normal_isotope() + k);
  }
  
  return 1;
}

/*
  Warning, this implementation is not complete. The definition from
  MDL says that if the M  CHG record is present, 0 charge is then
  forced on all atoms not mentioned. I'm not doing that. Change if
  this becomes a problem.
*/

int
Molecule::mdl_add_m_formal_charge (int ntokens,
                                                  const Aprop * atom_properties)
{
  assert (ntokens > 0);

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_charge: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }
    
    if (! reasonable_formal_charge_value(atom_properties[i]._property))
    {
      cerr << "mdl add formal charge: unreasonable charge value " <<
              atom_properties[i]._property << endl;
      return 0;
    }

    _things[j]->set_formal_charge(atom_properties[i]._property);
  }
  
  return 1;
}

int
Molecule::mdl_add_m_radical (int ntokens,
                                             const Aprop * atom_properties)
{
  assert (ntokens > 0);

  int matoms = _number_elements;

  for (int i = 0; i < ntokens; i++)
  {
    int j = atom_properties[i]._atom_number;
    j--;     // convert to C++ numbering
    if (j < 0 || j >= matoms)
    {
      cerr << "mdl_add_m_radical: illegal atom number " << j << ", " <<
               matoms << " atoms in the molecule\n";
      return 0;
    }

    _things[j]->set_implicit_hydrogens(0, 1);
    _things[j]->set_implicit_hydrogens_known(1);
  }
  
  return 1;
}

/*
  MDL files are somewhat strange in that there several varieties.
  Sometimes the connection table is terminated by 'M  END', sometimes by $$$$.
  Sometimes both are present.
  The way this works now is that by default, we insist on the $$$$.
  If return_on_m_end is set, then we return on 'M  END'
*/

int
Molecule::_read_molecule_mdl_ds (iwstring_data_source & input,
                                int return_on_m_end)
{
  int nb = 0;
  int v30;
  if (! _read_mdl_atom_connection_table(input, nb, v30))
  {
    skip_to_string(input, "$$$$", 1);

    return 0;
  }

  if(_bond_list.elements_allocated() < nb)
    _bond_list.resize(nb);

// Aromatic atoms are those found at the ends of aromatic bonds. Not ideal, but seems to work

  int * aromatic_atoms = NULL;
  int * aromatic_bonds = NULL;

  if (_number_elements > 0 && input_aromatic_structures())
  {
    aromatic_atoms = new_int(_number_elements);
    if (nb > 0)
      aromatic_bonds = new_int(nb);
  }

  int rc;
  if(v30)
    rc = _read_v30_bond_list(input, nb, aromatic_atoms, aromatic_bonds);
  else
    rc = _read_mdl_bond_list(input, nb, aromatic_atoms, aromatic_bonds);

  if(v30)
    skip_to_string(input, "M  END", 0);    // 0 means do it quietly

  if (NULL != aromatic_atoms)
  {
    if (! _final_processing_of_aromatic_mdl_input(aromatic_atoms, aromatic_bonds))
      rc = 0;

    delete [] aromatic_atoms;
    delete [] aromatic_bonds;
  }

  if(unconnect_covalently_bonded_non_organics_on_read())
    _do_unconnect_covalently_bonded_non_organics();

  if (v30 && return_on_m_end)    // possibly reading an rdf file
    return rc;

  if (! _read_molecule_mdl_trailing_records(input, return_on_m_end))
    return rwmolecule_error("read_molecule_mdl_ds: bad stuff at end", input);

  return rc;
}

/*
  We have an atom with > 2 connections. We need to know whether there are two
  neighbouring atoms that could be delocalised.
  The two neighbours must be set in AROMATIC_ATOMS, and must share the
  same atomic number and connectivity

  I needed this because we may have the case of a delocalised Carboxyllic
  acid adjacent to an aromatic ring. The entry in AROMATIC_ATOMS for the
  carbon in the ring will be set.

  c1ccccc1C(=O)O

  Mar 2001. Reading charged sulphonic acids from mol2 files, you can have
  three neighbours

  Oct 2001. Encountered O--S(--O)(--O)--O

  March 2006. Allow PO2 and PO3
*/

int
Molecule::_has_delocalised_neighbours (atom_number_t zatom,
                                       const int * aromatic_atoms,
                                       const int * aromatic_bonds,
                                       Set_of_Atoms & s) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();
  s.resize(acon);

  int invariant = 0;    //  a quick and dirty computation to see if two atoms are the same

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if(aromatic_atoms[j])
      ;
    else if (NULL != aromatic_bonds)   // check to see if the bond is marked aromatic
    {
      int b = which_bond(zatom, j);
      if (! aromatic_bonds[b])
        continue;
    }
    else
      continue;

    const Atom * aj = _things[j];

    int my_invariant;

    if (8 == aj->atomic_number())
      my_invariant = 10 * aj->atomic_number() + aj->ncon();   // breaks if we get more than 10 connections
    else if (7 == aj->atomic_number() && aj->ncon() < 3)   // guanidines
      my_invariant = 10 * aj->atomic_number();
    else
      continue;

//  cerr << "Atom " << j << " is an aromatic neighbour, invariant " << my_invariant << endl;

    if (0 == s.number_elements())    // first one
    {
      invariant = my_invariant;
      s.add(j);
      continue;
    }
    else if (my_invariant != invariant)    // different atom type
      continue;

    if (0 != aj->formal_charge())    // make sure anything with an existing formal charge up front
      s.insert_at_beginning(j);
    else
      s.add(j);
  }

//cerr << "Found " << s.number_elements() << " aromatic neighbours\n";

  if (0 == s.number_elements())
    return 0;

  if (NULL == aromatic_bonds)
    return(s.number_elements() > 1);

//  Implement this sometime, maybe not needed

  int nb = nedges();

  for (int i = 0; i < nb; i++)
  {
    if (0 == aromatic_bonds[i])
      continue;

    const Bond * b = _bond_list[i];

    if (! b->involves(zatom))
      continue;
  }

  return 1;
}

int
Molecule::_unset_aromatic_bond_entry(atom_number_t a1, 
                                     atom_number_t a2,
                                     int * aromatic_bonds) const
{
  if (NULL == aromatic_bonds)
    return 1;

  int b = which_bond(a1, a2);

  assert (b >= 0);

  aromatic_bonds[b] = 0;

//cerr << " Bond between " << a1 << " and " << a2 << " is bond " << b << endl;

  return 1;
}

int
Molecule::_more_than_one_aromatic_bond(atom_number_t zatom,
                                       const int * aromatic_bond) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int aromatic_bonds_found = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(zatom, i);

    int b = which_bond(zatom, j);

    if(aromatic_bond[b])
      aromatic_bonds_found++;
  }

  return aromatic_bonds_found > 1;
}

//#define DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS

/*
  Sept 97, generalise this to be any two equivalent looking 'aromatic'
  bonds bonded to the same atom.

  Mar 2001, update to include sulph* acids with more than 2 delocalised
  oxygens

  Feb 2009. ALso do Nitro groups
*/

int
Molecule::_process_delocalised_carbonyl_bonds (int * aromatic_atoms,
                                               int * aromatic_bonds)
{
  if (NULL == aromatic_bonds)    // cannot do anything
    return 1;

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
  if (NULL != aromatic_bonds)
  {
    cerr << "On entry to _process_delocalised_carbonyl_bonds\n";
    iw_write_array(aromatic_bonds, nedges(), "aromatic bonds", cerr);
  }
#endif

  for (int i = 0; i < _number_elements; i++)
  {
    if(aromatic_atoms[i])
      ;
    else if (_more_than_one_aromatic_bond(i, aromatic_bonds))
      ;
    else
      continue;

    if(is_ring_atom(i))
      continue;

    Atom * a = _things[i];

    int icon = a->ncon();

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
    cerr << "Molecule::_process_delocalised_carbonyl_bonds: atom " << i << " '" << a->atomic_symbol() << "' " << icon << " connections is non ring\n";
#endif

    if (icon < 2)
      continue;

    Set_of_Atoms s;

    if (! _has_delocalised_neighbours(i, aromatic_atoms, aromatic_bonds, s))
      continue;

    if (1 == s.number_elements())   // deal with these later
      continue;

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
    cerr << "Atoms " << s << " delocalised\n";
#endif

    aromatic_atoms[i] = 0;
    s.set_vector(aromatic_atoms, 0);

    atom_number_t a0 = s[0];

#ifdef NEED_TO_PROCESS_AROMATIC_PO_FROM_TRIPOS
// Jan 2005.  Sybyl produces P:O bonds.  No, just one bond in a PO3
// grouping is aromatic, so this isn't needed

    if (1 == s.number_elements())
    {
      if(a->implicit_hydrogens() && _things[a0]->implicit_hydrogens())
      {
        set_bond_type_between_atoms(i, a0, DOUBLE_BOND);
        a->set_modified();
        _things[a0]->set_modified();
        if (NULL != aromatic_bonds)
        {
          int b = which_bond(i, a0);
          aromatic_bonds[b] = 0;
        }
        continue;
      }

      cerr << "Molecule::_process_delocalised_carbonyl_bonds:no available bonds\n";
      continue;
    }
#endif

    _unset_aromatic_bond_entry(i, a0, aromatic_bonds);

//  Special case for the Nitro group

    if (7 == a->atomic_number() && 2 == s.number_elements() && 8 == atomic_number(a0) && 8 == atomic_number(s[1]))
    {
      set_bond_type_between_atoms(i, a0, DOUBLE_BOND);
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      _unset_aromatic_bond_entry(i, s[1], aromatic_bonds);
      continue;
    }

    set_bond_type_between_atoms(i, a0, SINGLE_BOND);

    if (-1 != _things[a0]->formal_charge())
      _things[a0]->set_formal_charge(-1);

    if (s.number_elements() > 1)
      _unset_aromatic_bond_entry(i, s[1], aromatic_bonds);

    if (2 == s.number_elements())    // the most common case, carbonyl for example
    {
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      continue;
    }

    _unset_aromatic_bond_entry(i, s[2], aromatic_bonds);

//  SO3 gets two double bonds, PO3 gets just one

    if (3 == s.number_elements())
    {
      set_bond_type_between_atoms(i, s[1], DOUBLE_BOND);
      if (16 == a->atomic_number())
        set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      else
        set_bond_type_between_atoms(i, s[2], SINGLE_BOND);
      continue;
    }

//  SO3 gets two double bonds, PO3 gets just one

    if (4 == s.number_elements())
    {
      set_bond_type_between_atoms(i, s[1], SINGLE_BOND);
      _things[s[1]]->set_formal_charge(-1);
      set_bond_type_between_atoms(i, s[2], DOUBLE_BOND);
      set_bond_type_between_atoms(i, s[3], DOUBLE_BOND);
      _unset_aromatic_bond_entry(i, s[3], aromatic_bonds);
    }
    else
      cerr << "Molecule::_process_delocalised_carbonyl_bonds: This many neighbours " << s << endl;
  }

#ifdef DEBUG_PROCESS_DELOCALISED_CARBONYL_BONDS
   cerr << "_process_delocalised_carbonyl_bonds returning\n";
#endif
  return 1;
}

//#define DEBUG_FINAL_PROCESSING_AROMATIC_MDL

int
Molecule::_final_processing_of_aromatic_mdl_input(int * aromatic_atoms,
                                int * aromatic_bonds)
{
#ifdef DEBUG_FINAL_PROCESSING_AROMATIC_MDL
  cerr << "Molecule::_final_processing_of_aromatic_mdl_input:\n";
#endif

  if(allow_delocalised_carbonyl_bonds())
  {
    if (! _process_delocalised_carbonyl_bonds(aromatic_atoms, aromatic_bonds))
    {
      cerr << "Cannot process delocalised carbonyl bonds\n";
      return 0;
    }
  }

  if (locate_item_in_array(1, _number_elements, aromatic_atoms) < 0)    // there were no aromatic atoms in the input
    return 1;

  int rc = 1;

  if (! find_kekule_form(aromatic_atoms))
  {
    cerr << "read_molecule_mdl_ds: find kekule form failed";

    if(allow_input_without_valid_kekule_form())
    {
      cerr << ", using single bonds";
      _molecule_name += " (invalid KEKULE form)";
    }
    else
      rc = 0;

    cerr << ' ' << _molecule_name << "'\n";
  }

  return rc;
}

/*
  A frequent operation is to decide whether a problem with bad chirality
  should be fatal or not
*/

static int
return_code_depending_on_ignore_incorrect_chiral_input()
{
  if(ignore_incorrect_chiral_input())
  {
    cerr << "Ignored\n";
    return 1;
  }

  return 0;
}

/*
  Main function for reading mdl files
*/

int
Molecule::read_molecule_mdl_ds (iwstring_data_source & input,
                                int return_on_m_end)
{
  assert(ok());

  resize(0);

  assert(input.good());
  if(input.eof())
    return 0;

  if (! _read_molecule_mdl_ds(input, return_on_m_end))
    return 0;

// Clean up chirality and cis-trans stuff
// Aug 2001.  If there are no chiral centres, but wedge bonds present,
// use that data - files from Afferent are like that!  Note that this
// isn't robust.  We could have some centres that are not atom marked,
// but have wedge bonds

  if(ignore_all_chiral_information_on_input())
    _chiral_centres.resize(0);
  else if(mdl_discern_chirality_from_wedge_bonds())
   (void) discern_chirality_from_wedge_bonds();
  else if (0 == _chiral_centres.number_elements() && number_up_or_down_wedge_bonds())
     (void) discern_chirality_from_wedge_bonds();
  else if (discern_chirality_from_3d_coordinates() && 3 == highest_coordinate_dimensionality())
  {
    int d = discern_chirality_from_3d_coordinates();
    if (1 == d)
      (void) discern_chirality_from_3d_structure();
    else if (_number_elements <= d)
      (void) discern_chirality_from_3d_structure();
    else
      cerr << "Molecule::_read_molecule_mdl_ds:skipped d@3d for too many atoms '" << name() << "' " << _number_elements << '\n';
  }

  if (0 == _chiral_centres.number_elements())    // none to worry about
    ;
  else if(_complete_chiral_centres_from_mdl_files())    // good
    ;
  else               // OOPS, bad chirality info
  {
    cerr << "Molecule::read_molecule_mdl_ds: erroneous chiral input\n";
    cerr << _molecule_name << endl;
    _chiral_centres.resize(0);
    if (! ignore_incorrect_chiral_input())
      return 0;
  }

// Do cis-trans bonds last because they depend on ring membership. Ran into a case, 583770, 
// where the ring perception forced a smiles ordering that was wrong because chirality hadn't
// been perceived

  if(discern_cis_trans_bonds())
   (void) discern_cis_trans_bonds_from_depiction();

  return 1;
}

/*
  Shared between the V2 and V3 programmes
*/

Atom *
create_mdl_atom (const const_IWSubstring & ss,
                 int msdif,
                 int chg,
                 int is_radical)
{
  const_IWSubstring zsymbol = ss;     // our own local copy

  if(ss.length() > 2)    // a residue or something
  {
    cerr << "Molecule::create_mdl_atom:warning: element '" << ss << "' encountered\n";

    if(change_long_symbols_to.length())
      zsymbol = change_long_symbols_to;
    else if(truncate_long_symbols)
      zsymbol.iwtruncate(2);
    else if(atomic_symbols_can_have_arbitrary_length())
      ;
    else               // no way of processing it
      return NULL;
  }

  int iso = 0;     // never used with mdl files

  const Element * e = get_element_from_symbol(zsymbol, iso);
  if (NULL == e)
  {
    if (allow_deuterium && 'D' == zsymbol)
    {
      e = get_element_from_symbol_no_case_conversion("H");
      iso = 2;
    }
    else if (allow_tritium && 'T' == zsymbol)
    {
      e = get_element_from_symbol_no_case_conversion("H");
      iso = 3;
    }
    else if (! auto_create_new_elements())
    {
      cerr << "create_mdl_atom: unrecognised element '" << zsymbol << "'\n";
      return NULL;
    }
    else
      e = create_element_with_symbol(zsymbol);

    if (NULL == e)
    {
      cerr << "Molecule::create_mdl_atom:cannot create element from '" << zsymbol << "'\n";
      return NULL;
    }
  }

  Atom * rc = new Atom(e);
  assert (NULL != rc);

  if (0 == msdif)
    ;
  else if(read_isotopes_as_numbers_rather_than_differences_from_normal)
    iso = msdif;
  else
    iso = e->normal_isotope() + msdif;

  if(iso)
    rc->set_isotope(iso);

  if(chg)
  {
    if (! reasonable_formal_charge_value(chg))
    {
      cerr << "create mdl atom: bad charge value " << chg << endl;
      return NULL;
    }
    rc->set_formal_charge(chg);
  }

  if(is_radical)
  {
    rc->set_implicit_hydrogens(0, 1);
    rc->set_implicit_hydrogens_known(1);
  }

  return rc;
}

/*static int
parse_xyz_symbol_etc (const const_IWSubstring & buffer,
                      coord_t & x,
                      coord_t & y,
                      coord_t & z,
                      const_IWSubstring & atomic_symbol,
                      int & msdiff,
                      int & chg,
                      int & astere)
{
  int nw = buffer.nwords();

  if (nw < 4)
  {
    cerr << "parse_xyz_symbol_etc:too few tokens\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(x))
  {
    cerr << "Bad X value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(y))
  {
    cerr << "Bad Y value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (token, i);
  if (! token.numeric_value(z))
  {
    cerr << "Bad Z value '" << token << "'\n";
    return 0;
  }

 (void) buffer.nextword (atomic_symbol, i);

  if (! buffer.nextword (token, i))
    return 1;

  if (! token.numeric_value(msdiff))
  {
    cerr << "Bad msdiff value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword (token, i))
    return 0;

  if (! token.numeric_value(chg))
  {
    cerr << "Bad chg value '" << token << "'\n";
    return 0;
  }

  if (! buffer.nextword (token, i))
    return 1;

  if (! token.numeric_value(astere))
  {
    cerr << "Bad astere value '" << token << "'\n";
    return 0;
  }

  return 1;
}*/

/*
  The reading of an MDL connection table is done in two functions.
  One reads the atom records, the other reads the bonds.
  The atom reading function is shared with the possibly aromatic routine.
  If it turns out we read a v30 connection table we set the variable V30
*/

int
Molecule::_read_mdl_atom_connection_table (iwstring_data_source & input, 
                                           int & nb,
                                           int & v30)
{
//  There are three header lines.
//  The line following will contain na and nb

// Note that the record with na and nb is being read here too!

  nb = 0;
  v30 = 0;

  const_IWSubstring buffer;
  for (int i = 0; i < 4; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");
    if (0 == i && ! discard_sdf_molecule_name)
      set_name(buffer);
  }

// buffer should now hold na and nb

  if (buffer.length() < 6)
  {
    cerr << "Molecule::_read_mdl_atom_connection_table: the atoms/bond record must be at least 6 chars long, line " << input.lines_read() << endl;
    cerr << buffer << endl;
    return 0;
  }

  int na;
  if (2 != int3d(buffer, na, nb))
  {
    cerr << "Molecule::_read_mdl_atom_connection_table: error from int3d '" << buffer << "'\n";
    cerr << "Line " << input.lines_read() << endl;
    return 0;
  }

//cerr << "Contains '" << na << " atoms and " << nb << " bonds\n";

  assert (na >= 0 && (nb >= 0));

  if (_elements_allocated < na)
    resize(na);

  if(accumulate_mdl_chirality_features)
  {
    unspecified_chiral_atoms_last_molecule.resize_keep_storage(0);
    down_bonds_last_molecule.resize_keep_storage(0);
    up_bonds_last_molecule.resize_keep_storage(0);
    squiggle_bonds_last_molecule.resize_keep_storage(0);
    unspecified_double_bond_atoms.resize_keep_storage(0);
  }

  if (buffer.contains(" V3000"))
  {
    v30 = 1;
    return _read_mdl_atom_connection_table_v30(input, nb);
  }

  MDL_Atom_Record mdl_atom_record;

  for (int i = 0; i < na; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

    if (! mdl_atom_record.build(buffer))
    {
      cerr << buffer << endl;
      return rwmolecule_error("read_molecule_mdl_ds:bad atom data", input);
    }

    Atom * a = mdl_atom_record.create_atom();

    if (NULL == a)
    {
      cerr << buffer << endl;
      cerr << rwmolecule_error("read_molecule_mdl_ds:bad element", input);
      return 0;
    }

    add(a);

    if (0 == mdl_atom_record.astere())     // the most common case, no chirality
      ;
    else if (! _mdl_atom_is_chiral_centre(_number_elements - 1, mdl_atom_record.astere()))
    {
      cerr << "Molecule::_read_mdl_atom_connection_table: invalid chirality on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

/*
  When reading an atom record we have encountered a chirality flag
*/

int
Molecule::_mdl_atom_is_chiral_centre (atom_number_t zatom, int cfg)
{
  const Atom * a = _things[zatom];

  if (! a->element()->organic())
  {
    if (3 == cfg)     // unspecified, who cares
      return 1;

    if (display_non_organic_chirality_messages)
      cerr << "Ignoring non-organic chirality, atom " << zatom << " type '" << a->atomic_symbol() << "'\n";
    return 1;
  }

  if (1 == cfg || 2 == cfg)
  {
    Chiral_Centre * c = create_chiral_centre(zatom, 1);
    if (NULL == c)
    {
      cerr << "Cannot create chiral center at atom " << zatom << endl;
      return 0;
    }
    c->set_chirality_known(cfg);
  }
  else if (3 == cfg)    // maybe ignore these
  {
    if(accumulate_mdl_chirality_features)
      unspecified_chiral_atoms_last_molecule.add(zatom);
  }
  else
  {
    cerr << "Molecule::_mdl_atom_is_chiral_centre: unrecognised chirality " << cfg << endl;
    return 0;
  }

  return 1;
}

/*
  We have read a bond record and DIRECTIONALITY is the 4'th number from the
  MDL bond record. Set the appropriate directionality flags in the appropriate
  bond
*/

int
Molecule::_mdl_set_bond_directionality (atom_number_t a1, atom_number_t a2,
                                        int directionality)
{
  assert (0 != directionality);

  for (int i = _bond_list.number_elements() - 1; i >= 0; i--)    // chances are this is the last bond in the list
  {
    Bond * b = _bond_list[i];
    if (! b->involves(a1, a2))
      continue;

    if(b->is_wedge_any())
    {
      cerr << "Molecule::_mdl_set_bond_directionality: bond between atoms " << a1 << " and " << a2 << " already directional\n";
      continue;
    }

    if (1 == directionality)
      b->set_wedge_up();
    else if (6 == directionality)
      b->set_wedge_down();
    else if (4 == directionality)
      b->set_wedge_either();
    else if (3 == directionality)
      b->set_cis_trans_either_double_bond();
    else
    {
      cerr << "Molecule::_mdl_set_bond_directionality: unknown directionality " << directionality << endl;
      return 0;
    }

    return 1;
  }

  cerr << "Molecule::_mdl_set_bond_directionality: HUH, no bond between atoms " << a1 << " and " << a2 << endl;
  return 0;
}

/*
  User visible numbering for bonds is
    single = 1
    double = 2
    triple = 3
    arom   = 4

  This converts those numbers into bond_type_t forms.
*/

int
convert_from_mdl_number_to_bond_type (int int_rep, bond_type_t & bt)
{
  switch(int_rep)
  {
    case 1:
      bt = SINGLE_BOND;
      return 1;

    case 2:
      bt = DOUBLE_BOND;
      return 1;

    case 3:
      bt = TRIPLE_BOND;
      return 1;

    case 4:
      bt = AROMATIC_BOND;
      return 1;

    default:
      cerr << "convert_to_bond_type: unrecognised type " << int_rep << endl;
      return 0;
  }

  assert (NULL == "Should not come here");
  return 0;
}

int
parse_bond_record (const_IWSubstring & buffer,
                   int na,
                   atom_number_t & a1, atom_number_t & a2,
                   int & bond_type_read_in,
                   int & directionality)
{
  if(buffer.length() < 9)
  {
    cerr << "Invalid bond record '" << buffer << "'\n";
    return 0;
  }

  directionality = 0;

  int ntokens = int3d(buffer, a1, a2, &bond_type_read_in);

  a1--;       // our atom numbers start at 0
  a2--;

  if ( (3 != ntokens) || a1 < 0 || a1 >= na || a2 < 0 || a2 >= na ||
       ((a1 == a2) && ! _ignore_self_bonds))
  {
    cerr << "parse_bond_record: error on bond record\n";

    if (3 != ntokens)
      cerr << "Bad token count on bond record " << ntokens << endl;
    else if (a1 < 0 || a1 > na)
      cerr << "Bad a1 value " << a1 << ", natoms = " << na << endl;
    else if (a2 < 0 || a2 > na)
      cerr << "Bad a2 value " << a2 << ", natoms = " << na << endl;
    else
      cerr << "a1 = " << a1 << " a2 = " << a2 << " natoms = " << na << endl;
    return 0;
  }

  if(buffer.length() < 12)       // no directionality present
    return 1;

  char direction = buffer[11];     // the 12'th column - probably should check columns 10 and 11....

  if ('0' == direction)            // no directionality
    return 1;

  directionality = direction - '0';    // the numeric form we return in our argument list

  if (0 == accumulate_mdl_chirality_features)
    return 1;

  if ('1' == direction)
    up_bonds_last_molecule.add(a1);
  else if ('6' == direction)
    down_bonds_last_molecule.add(a1);
  else if ('4' == direction)
    squiggle_bonds_last_molecule.add(a1);
  else if ('3' == direction)
  {
    assert (2 == bond_type_read_in);
    unspecified_double_bond_atoms.add(a1);
    unspecified_double_bond_atoms.add(a2);
  }
  else
  {
    cerr << "parse_bond_record: unrecognised bond directionality '" << direction << "'\n";
    return 0;
  }

  return 1;
}

/*
  We always attempt to read NB records.
*/

int
Molecule::_read_mdl_bond_list (iwstring_data_source & input, int nb,
                               int * aromatic_atoms,
                               int * aromatic_bonds)
{
  assert (nb >= 0);

  int na = _number_elements;

  int rc = 1;    // assume ok until proven otherwise

  const_IWSubstring buffer;
  for (int i = 0; i < nb; i++)
  {
    EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

    if (0 == rc)     // just read the records if we have already failed
      continue;

    int a1, a2;
    int bond_type_read_in;
    int directionality;

    if (! parse_bond_record(buffer, na, a1, a2, bond_type_read_in, directionality))
    {
      cerr << "Molecule::_read_mdl_bond_list: bond record " << i << " is bad '" <<
              buffer << "'\n";
      rc = 0;
      continue;
    }

    if (NULL != input_bond_type_translation_table && bond_type_read_in < 20)
      bond_type_read_in = input_bond_type_translation_table[bond_type_read_in];

    if (_ignore_self_bonds && (a1 == a2))
    {
      cerr << "Molecule::_read_mdl_bond_list: ignoring self bonds " << a1 << " to " << a2 << endl;
      continue;
    }

    bond_type_t btype = INVALID_BOND_TYPE;

    if (4 == bond_type_read_in)
    {
      if (! input_aromatic_structures())
      {
        cerr << "Molecule::_read_mdl_bond_list:aromatic input not enabled, bond between atoms " << a1 << " and " << a2 << endl;
        rc = 0;
        continue;
      }

      btype = SINGLE_BOND;
      aromatic_atoms[a1] = 1;
      aromatic_atoms[a2] = 1;
      aromatic_bonds[i] = 1;
    }
    else if (! convert_from_mdl_number_to_bond_type(bond_type_read_in, btype))
    {
      cerr << "Molecule::_read_mdl_bond_list: bad bond type " << bond_type_read_in << endl;
      rc = 0;
      continue;
    }

    add_bond(a1, a2, btype, 1);     // 1 means partially built molecule
    if (directionality)
      _mdl_set_bond_directionality(a1, a2, directionality);
  }

  if (rc && nb > 0)
  {
    check_bonding();
  }

  return rc;
}

static int
parse_m_sty_record (const const_IWSubstring & buffer,
                    int & sgroups_present)
{
  if(buffer.nwords() < 3)
  {
    cerr << "Molecule::_read_molecule_mdl_trailing_records:invalid record '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring n;
  buffer.word(2, n);

  int tmp;

  if (! n.numeric_value(tmp) || tmp < 1)
  {
    cerr << "Molecule::parse_m_sty_record:invalid count '" << buffer << "'\n";
    return 0;
  }

  sgroups_present += tmp;

  return 1;
}

static int
extract_sdf_identifier(const IWString & buffer, IWString & id)
{
  int istart = buffer.index('<');
  int istop = buffer.rindex('>');

  if (istart > istop || istart < 0 || istop < 0)
  {
    cerr << "extract_sdf_identifier: HUH: '" << buffer << "'\n";
    return 0;
  }

  buffer.from_to(istart + 1, istop - 1, id);

  if (id.contains(' '))
    id.gsub(' ', '_');

  return 1;
}

/*
  Nov 99. GER needed to flag superatoms. These are things like

OMe
  -ISIS-  10279910532D

  4  3  0  0  0  0  0  0  0  0999 V2000
    7.5837   -2.1040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6320   -3.1626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9044   -2.4178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9458   -3.4454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
A    4
MeO
M  END
>  <REGNO>(2)
2

$$$$

  which is an example from Jeff Hanson, 27 Oct 99

  Note that we also flag records starting with a G
*/

static int a_records_found = 0;
static int g_records_found = 0;

static resizable_array_p<Atom_Alias> aliases;

static int _set_elements_based_on_atom_aliases = 0;

int
set_elements_based_on_atom_aliases()
{
  return _set_elements_based_on_atom_aliases;
}

void
set_set_elements_based_on_atom_aliases (int s)
{
  _set_elements_based_on_atom_aliases = s;
}

int
a_records_found_most_recent_mdl_molecule()
{
  return a_records_found;
}

int
g_records_found_most_recent_mdl_molecule()
{
  return g_records_found;
}

const resizable_array_p<Atom_Alias> &
atom_aliases_most_recent_mdl_molecule()
{
  return aliases;
}

void
add_to_text_info (resizable_array_p<IWString> & text_info, const const_IWSubstring & zextra)
{
  IWString * s = new IWString(zextra);

  text_info.add(s);

  return;
}

/*
  Sept 2010. Turns out mdl is now putting digits in these records, so we can have. Hmmm,
  this showed up earlier, not sure why we did not find problems earlier.

> 2  <MODEL.EXTREG>

*/

#ifdef LOOKS_LIKE_SDF_TAG_OV
static int
looks_like_sdf_tag (const const_IWSubstring & buffer)
{
//cerr << "Checking for sdf tag '" << buffer << "'\n";

  if (! buffer.starts_with('>'))
    return 0;

  int n = buffer.length();

  if (n < 3)
    return 0;

  for (int i = 1; i < n; i++)
  {
    char c = buffer[i];

    if (' ' == c)
      continue;

    if ('<' != c)
      return 0;

    do
    {
      if ('>' == buffer[i])
        return 1;

      i++;
    } while (i < n);

    return 0;
  }

  return 0;
}
#endif

static int
looks_like_sdf_tag (const const_IWSubstring & buffer)
{
//cerr << "Checking for sdf tag '" << buffer << "'\n";

  if (! buffer.starts_with('>'))   // cheapest test first
    return 0;

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  if (! buffer.nextword(token, i))   // strange, no second token on line!
    return 0;

  if (token.starts_with('<'))   // the most common case
    ;
  else         // maybe a digit in there
  {
    int notused;
    if (! token.numeric_value(notused) || notused < 0)
      return 0;

    if (! buffer.nextword(token, i))   // if we have a digit, there must be a following token
      return 0;

    if (! token.starts_with('<'))   // following the digit, there must be the 
      return 0;
  }

  if (token.ends_with('>'))
    return 1;

  while (buffer.nextword(token, i))
  {
    if (token.ends_with('>'))
      return 1;
  }

  return 0;
}

/*
  By convention, data following a tag may be multi-record, but will be
  terminated by a blank line
*/

int
Molecule::_read_mdl_data_following_tag (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  for (int i = 0; input.next_record(buffer); i++)
  {
    if (looks_like_sdf_tag (buffer))
    {
      input.push_record();
      return 1;
    }

    if (read_extra_text_info())   // even if buffer is empty
      _text_info.add(new IWString(buffer));

    if (0 == buffer.length() || (1 == buffer.length() && 13 == static_cast<int>(buffer[0])))    // 13 is ^M
      return 1;

    if (0 == i && prepend_sdfid)   // on first record, no space before extra info
      _molecule_name += buffer;
    else
      _molecule_name.append_with_spacer(buffer, insert_between_sdf_name_tokens);
  }

  cerr << "Molecule::_read_mdl_data_following_tag:premature eof\n";
  return 0;
}

/*
  the MDL_Molecule class inherits from a Molecule. It needs to
  examine all the M records looking for substructural info.
  This function will process any records associated with the
  Molecule class
*/

int
Molecule::_common_parse_M_record (const const_IWSubstring & buffer,
                                  int & fatal)
{
  Aprop atom_properties[MAX_PAIRS];

#ifdef DEBUG_COMMON_PARSE_M_RECORD
  cerr << "_common_parse_M_record examining '" << buffer << "'\n";
#endif

  if (buffer.starts_with("M  CHG"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_formal_charge(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  ISO"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_isotope(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  RAD"))
  {
    int tokens;
    if (! fill_atom_property_array(buffer, tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    if (0 == tokens)
      return 1;

    if (! mdl_add_m_radical(tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with("M  RGP"))
  {
    if (! _parse_M_RGP_record(buffer))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

#ifdef DEBUG_COMMON_PARSE_M_RECORD
  cerr << "_common_parse_M_record:did not recognise '" << buffer << "'\n";
#endif

// Not recognised here

  fatal = 0;
  return 0;
}

/*
  Function for reading all the "stuff" which comes between the bond
  list and the end of the molecule
*/

int
Molecule::_read_molecule_mdl_trailing_records (iwstring_data_source & input,
                                  int return_on_m_end)
{
  if (read_extra_text_info() && 0 == _text_info.number_elements())
    _text_info.resize(10);

  a_records_found = 0;
  g_records_found = 0;
  aliases.resize_keep_storage(0);

//Aprop atom_properties[MAX_PAIRS];

// We just read and ignore the S group information

  int sgroups_present = 0;
  int msal = 0;
  int msbl = 0;
  int msmt = 0;
  int msbv = 0;

  IWString buffer;
  int trailing_lines = 0;
  int extreg_found = 0;
  int got_dollar = 0;

  while (input.next_record(buffer))
  {
    trailing_lines++;

    buffer.strip_trailing_blanks();      // do this first, get rid of any garbage out the end

    if ("$$$$" == buffer)
    {
      got_dollar = 1;
      break;
    }

    int fatal;
    if (_common_parse_M_record(buffer, fatal))   // great, recognised and good
      continue;
    else if (fatal)
    {
      cerr << "Molecule::_read_molecule_mdl_trailing_records:invalid record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    if (buffer.starts_with("M  STY"))
    {
      if (! parse_m_sty_record(buffer, sgroups_present))
      {
        if(die_on_erroneous_m_input)
          return 0;
      }

      continue;
    }

    if (buffer.starts_with("M  SAL"))
    {
      msal--;
      continue;
    }

    if (buffer.starts_with("M SBL"))
    {
      msbl--;
      continue;
    }

    if (buffer.starts_with("M  SBV"))
    {
      msbv--;
      continue;
    }

    if (buffer.starts_with("M  SMT"))
    {
      msmt--;
      continue;
    }

    if ("M  END" == buffer)
    {
      if(return_on_m_end)
        return 1;
      continue;
    }

    if (0 == extreg_found && name_in_m_tag.length() && buffer.nwords() > 2 &&
        buffer.starts_with("M ") && buffer.contains(name_in_m_tag))
    {
      _molecule_name = buffer;
      _molecule_name.remove_leading_words(2);    // M REG
      extreg_found = 1;
      continue;
    }

    if (buffer.starts_with("M  "))
    {
      if (report_unrecognised_records || ! ignore_unrecognised_m_records)
        cerr << "Unrecognised 'M  ' directive '" << buffer << "'\n";
      if (! ignore_unrecognised_m_records)
        return 0;
    }

    if (buffer.starts_with("A  "))
    {
      Atom_Alias * a = new Atom_Alias;
      if (! a->build(buffer, input))
      {
        cerr << "Invalid atom alias data, line " << input.lines_read() << endl;
        delete a;
        return 0;
      }

      aliases.add(a);

      a_records_found++;
     (void) input.next_record(buffer);

      if(buffer.length() > 0)
        input.push_record();

      continue;
    }

    if (buffer.starts_with("G  "))
    {
      g_records_found++;
      IWString g(buffer);   // needs two records, make a copy of first one
      (void) input.next_record(buffer);

      if (mdl_g_records_hold_atom_symbols)
      {
        if (! _process_mdl_g_record (g, buffer))
        {
          cerr << "Molecule::_read_molecule_mdl_trailing_records:cannot process G record '" << buffer << "'\n";
          return 0;
        }
      }

      continue;
    }

    if (read_extra_text_info())
      add_to_text_info(_text_info, buffer);

//  cerr << read_extra_text_info() << " now contains " << _text_info.number_elements() << endl;

    if (0 == buffer.length())
      continue;

//  Now all the various other identifiers possible in the file

    if (sdf_identifier.active() && sdf_identifier.matches(buffer))
    {
      IWString id;
      extract_sdf_identifier(buffer, id);

      EXTRA_STRING_RECORD (input, buffer, "read mol mdl");

      if(read_extra_text_info())
        add_to_text_info(_text_info, buffer);

      IWString tmp;
      if (replace_first_sdf_tag.length() > 0)
        tmp << replace_first_sdf_tag << ':';
      else if(prepend_sdfid)
        tmp << id << ':';

      tmp += buffer;

      _molecule_name.append_with_spacer(tmp, insert_between_sdf_name_tokens);

      continue;
    }

//  Are we fetching all the SDF data from the file.
//  We run into lots of different formats.
//  >  <IDENTIFIER> 
//  data
//
// >  <IDENTIFIER> stuff
// data

// We need to fetch IDENTIFIER

    if (fetch_all_sdf_identifiers && looks_like_sdf_tag(buffer))
    {
      IWString id;
      if (! extract_sdf_identifier(buffer, id))
        return rwmolecule_error("read_molecule_mdl_ds: cannor parse SDF identifier", input);

      if (replace_first_sdf_tag.length() > 0)
      {
        _molecule_name.append_with_spacer(replace_first_sdf_tag, insert_between_sdf_name_tokens);
        _molecule_name.add(':');
      }
      else if (prepend_sdfid)
      {
        _molecule_name.append_with_spacer(id, insert_between_sdf_name_tokens);
        _molecule_name.add(':');
      }

      if (! _read_mdl_data_following_tag(input))
        return 0;

      continue;
    }

    if (0 == _molecule_name.length() && take_first_tag_as_name &&
        looks_like_sdf_tag(buffer))
    {
      if (! _read_mdl_data_following_tag(input))
        return 0;

      continue;
    }

//  If we are storing the extra info, do so, otherwise silently ignore other
//  kinds of records.

    if (0 == extreg_found && extract_isis_extregno &&
        buffer.starts_with(">  <") && 3 == buffer.nwords())
    {
      IWString extreg;
      buffer.word(2, extreg);
      if (extreg.starts_with('(') && extreg.ends_with(')'))
      {
        extreg.remove_leading_chars(1);
        extreg.chop(1);

        if (0 == _molecule_name.length())
          _molecule_name = extreg;
        else
        {
          extreg += ' ';     // allow a space before the current name
          _molecule_name.insert(extreg, 0);
        }

        extreg_found = 1;
        continue;
      }
    }

//  If we get to here, just ignore it.

    if(report_unrecognised_records)
      cerr << "Ignoring unrecognised form, line " << input.lines_read() << " '" << buffer << "'\n";
  }

  if(set_elements_based_on_atom_aliases() && aliases.number_elements())
    _set_elements_based_on_atom_aliases(aliases);

  if(got_dollar)
    return 1;

// If we come out here, it must be EOF. That's OK.

//if (trailing_lines > 6)     May 2005, this warning doesn't seem necessary
//  cerr << "mdl_read_ds: " << trailing_lines << " lines found between bonds and EOF\n";

  cerr << "mdl_read_ds returning at EOF without $$$$\n";

  return 1;
}

/*
  When reading RDFILES, the user can specify which fields to use for the name
  These are processed in order
*/

static resizable_array_p<IWString> rdfile_identifiers;

int
add_rdfile_identifier (const IWString & new_identifier)
{
  assert(new_identifier.length() > 0);

  IWString * tmp = new IWString(new_identifier);

  return rdfile_identifiers.add(tmp);
}

static IWString rdfile_start_of_record;

int
set_rdfile_start_of_record (const const_IWSubstring & s)
{
  rdfile_start_of_record = s;

//cerr << "rdfile_start_of_record set to '" << s << "'\n";

  return 1;
}

static int
skip_to_rdfile_start_of_record (iwstring_data_source & input,
                                const IWString & rdfile_start_of_record)
{
  const_IWSubstring buffer;
  
  int records_read_here = 0;

  while (input.next_record(buffer))
  {
    records_read_here++;

    if (buffer.starts_with(rdfile_start_of_record))
    {
      input.push_record();
      return 1;;
    }
  }

  if (0 == records_read_here)
  {
    cerr << "read mol rdf eof\n";
    return 0;
  }

  cerr << "EOF reading RDFILE, cannot find start record '" << rdfile_start_of_record << "', tried " << records_read_here << "\n";
  return 0;
}

//#define DEBUG_READ_MOLECULE_RDF_DS

/*
  Lots of ways an RDFILE identifier can match
*/

static int
rdfile_record_matches (const IWString & buffer,
                       const IWString & rdfile_identifier)
{
#ifdef DEBUG_READ_MOLECULE_RDF_DSQ
  cerr << "Checking rdfile match '" << buffer << "' vs '" << rdfile_identifier << "'\n";
#endif

  if (buffer.starts_with(rdfile_identifier))
    return 1;

  if (buffer.starts_with("$DTYPE MOL:") && buffer.matches_at_position(11, rdfile_identifier))
    return 1;

  if (buffer.starts_with("$DTYPE ") && buffer.matches_at_position(7, rdfile_identifier))
      return 1;

  return 0;
}

/*
  An RDFILE is somewhat different, in that it has a bunch of records starting with $

  The connection table starts with '$DATUM $MFMT'
  Also have seen files where it starts with '$MFMT $MIREG' sometimes.
  We must tell read_molecule_mdl_ds to return on encountering 'M  END'
*/

int
Molecule::read_molecule_rdf_ds (iwstring_data_source & input)
{
  assert(ok());
  assert(input.good());

  if(input.eof())
    return 0;

#ifdef DEBUG_READ_MOLECULE_RDF_DSQ
  cerr << "Molecule::read_molecule_rdf_ds:input.record_buffered? " << input.record_buffered() << endl;
#endif

  if(rdfile_start_of_record.length())   // scan through file till we get to header
  {
    if (! skip_to_rdfile_start_of_record(input, rdfile_start_of_record))
      return 0;
  }

  IWString possible_name;

  if (! _read_molecule_rdf_ds(input, possible_name))
    return 0;

  if (0 == _molecule_name.length() && possible_name.length())
    set_name(possible_name);

  return 1;
}

int
Molecule::_read_molecule_rdf_ds (iwstring_data_source & input,
                                 IWString & possible_name)
{
  int nrdfid = rdfile_identifiers.number_elements();

// I've seen many kinds of RDF files. Some have the name in the name
// field, others have it as a data item of type LILLYNUMBER. Look for
// LILLYNUMBER things and save the next record.

  int records_read_here = 0;

  IWString buffer;    // must be an IWString, NOT a const_IWSubstring

  int have_read_structure = 0;

  while (input.next_record(buffer))
  {
    records_read_here++;

#ifdef DEBUG_READ_MOLECULE_RDF_DS
    cerr << "Just read '" << buffer << "'\n";
#endif

    if (buffer == "$DATUM $MFMT" || buffer.starts_with("$MFMT $MIREG"))
    {
      if(have_read_structure)
      {
        cerr << "Molecule::read_molecule_rdf_ds:multiple structures!!\n";
        return 0;
      }

#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Start reading mdl connection table at " << input.lines_read() << endl;
#endif

      if (! read_molecule_mdl_ds(input, 1))
        return 0;

      have_read_structure = 1;

#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Successfully read structure\n";
#endif

      if (0 == rdfile_start_of_record.length())
        return 1;
    }
    else if(nrdfid)
    {
      for (int i = 0; i < nrdfid; i++)
      {
        const IWString & rdfile_identifier = *(rdfile_identifiers[i]);

        if (rdfile_record_matches(buffer, rdfile_identifier))
        {
          const_IWSubstring tmp;
          if (! input.next_record(tmp))
          {
            cerr << "Molecule::read_molecule_rdf_ds: premature eof\n";
            return 0;
          }

          if (tmp.starts_with("$DATUM"))
            tmp.remove_leading_words(1);

          possible_name.append_with_spacer(rdfile_identifier);
          possible_name << ':' << tmp;
        }
      }
    }
    else if (buffer.contains("$DTYPE ") && buffer.contains(":LILLYNUMBER"))
    {
     (void) (input.next_record(possible_name));
      if (! possible_name.contains("$DATUM "))    // strip this some time.
      {
        cerr << "Molecule::read_molecule_rdf_ds: Yipes, no $datum '" << possible_name << "'\n";
        return 0;
      }
    }

#ifdef DEBUG_READ_MOLECULE_RDF_DS
    cerr << "Checking start '" << rdfile_start_of_record << "' and '" << buffer << "'\n";
#endif
    if(rdfile_start_of_record.length() && buffer.starts_with(rdfile_start_of_record) && records_read_here > 1)
    {
#ifdef DEBUG_READ_MOLECULE_RDF_DS
      cerr << "Found rdf start of record, structure? " << have_read_structure << endl;
#endif

      if (! have_read_structure)
        cerr << "Molecule::read_molecule_mdl_ds:no structure!\n";

      input.push_record();
      return 1;
    }
  }

  if (0 == records_read_here || have_read_structure)
  {
    cerr << "read mol rdf eof\n";
    return 0;
  }

  if(have_read_structure)
    return 1;

  if(possible_name.length())
    return 1;

  cerr << "EOF reading RDFILE\n";
  return 0;
}

int
Molecule::write_molecule_mdl (const char * fname, const char * comments) const
{
  assert(ok());
  assert (NULL != fname);

  ofstream output(fname);
  if (! output.good())
  {
    cerr << "Molecule::write_molecule_mdl: cannot open '" << fname << "'\n";
    return 0;
  }

  return write_molecule_mdl(output, comments);
}

/*
  Function to write an MDL molfile.
*/

// Feb 97, changed to write atom properties as M form by default

static int write_mdl_charges_as_m_chg = 1;

void
set_write_mdl_charges_as_m_chg (int i)
{
  write_mdl_charges_as_m_chg = i;

  return;
}

int
Molecule::write_molecule_mdl (ostream & os, const IWString & comments) const
{
  assert(ok());
  assert(os.good());

  if (_number_elements > 999 || write_v30_mdl_files)
    return write_molecule_mdl_v30(os, comments);

  os << _molecule_name << newline_string();

  if(isis_standard_records)
  {
    int dim;
    if(highest_coordinate_dimensionality() > 2)
      dim = 3;
    else 
      dim = 2;

    os << "  -ISIS-  0516971354" << dim << "D 1   1.00000     0.00000     1" << newline_string();
    os << newline_string();
  }
  else
  {
    if(comments.length())
      os << comments << newline_string();
    else
      os << "Blank" << newline_string();

    os << "Blank" << newline_string();
  }

  int rc = write_connection_table_mdl(os);

//if(isis_standard_records)
//  os << "M  END\n";

  if(::write_extra_text_info())
    write_extra_text_info(os);

  if(write_mdl_dollar_signs)
    os << "$$$$" << newline_string();

  if(flush_files_after_writing_each_molecule())
    os.flush();

  return rc;
}

int
Molecule::write_extra_text_info (ostream & os) const
{
  int ne = _text_info.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const IWString * info = _text_info[i];
    os <<(*info) << newline_string();
  }

  return os.good();
}

int
Molecule::write_extra_text_info (IWString & buffer) const
{
  int ne = _text_info.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const IWString * info = _text_info[i];

    buffer <<(*info) << newline_string();
  }

  return buffer.length();
}


/*
  Return the digit - 1,2,3 that is the stereo flag for a given atom
*/

int
Molecule::_mdl_atom_stereo_value (atom_number_t zatom) const
{
  Chiral_Centre * c = chiral_centre_at_atom(zatom);

  if (NULL == c)
    return 0;

  if (! c->chirality_known())
    return 3;

// If there is an explicit Hydrogen, that atom must be the highest numbered atom
// Who knows what you are supposed to do if there is an implicit and an explicit Hydrogen

  if(c->implicit_hydrogen_count())        // we assume no explicit hydrogen also
    return c->mdl_stereo_centre_value();

  if(c->lone_pair_count())                // who knows how we are supposed to handle a lone pair and an explicit hydrogen
    return c->mdl_stereo_centre_value();

  if (! _mdl_write_h_correct_chiral_centres)    // optional behaviour
    return c->mdl_stereo_centre_value();

// Look for an explicit Hydrogen atom. Remember, we pass the atoms SW N SE

  if (1 == _things[c->top_front()]->atomic_number())
    return c->mdl_stereo_centre_value(c->right_down(), c->top_back(), c->left_down());
  if (1 == _things[c->top_back()]->atomic_number())
    return c->mdl_stereo_centre_value(c->left_down(), c->top_front(), c->right_down());
  if (1 == _things[c->left_down()]->atomic_number())
    return c->mdl_stereo_centre_value(c->right_down(), c->top_front(), c->top_back());
  if (1 == _things[c->right_down()]->atomic_number())
    return c->mdl_stereo_centre_value(c->top_back(), c->top_front(), c->left_down());

  return c->mdl_stereo_centre_value();
}

/*
  For efficiency, we keep an array of 3 column digits
*/

static IWString * digits2 = NULL;
static IWString * digits3 = NULL;

static void
initialise_digits()
{
  digits2 = new IWString[100];
  digits3 = new IWString[1000];

  assert (NULL != digits3);

  for (int i = 0; i <= 999; i++)
  {
    IWString & d = digits3[i];

    if (i >= 100)
      d << i;
    else if (i >= 10)
      d << ' ' << i;
    else
      d << "  " << i;
  }

  for (int i = 0; i <= 99; i++)
  {
    IWString & d = digits2[i];
    if (i >= 10)
      d << i;
    else
      d << ' ' << i;
  }

  return;
}

/*
  I tried to have a class that would initialise the arrays of digits, but
  that doesn't work because there is no order in which various things
  get initialised. Therefore all we can do is have something that deletes
  the arrays if they have been formed
*/

class Digit_Deleter
{
  private:
  public:
    ~Digit_Deleter();
};

static void
delete_digits()
{
  if (NULL != digits2)
  {
    delete [] digits2;
    digits2 = NULL;
  }
  if (NULL != digits3)
  {
    delete [] digits3;
    digits3 = NULL;
  }

  return;
}

Digit_Deleter::~Digit_Deleter()
{
  delete_digits();
}

static Digit_Deleter digit_deleter;

void
delete_digits_objects_in_mdl_file ()
{
  delete_digits ();
}

/*
  One column of the SD file is the isotopic specification. The standard says
  only write as a number if it is in the range -3 +4. We don't really do it
  that way
*/

static void
append_isotope_information (IWString & output_buffer,
                            int iso,
                            int normal_isotope)
{
  int to_be_written;

  if(write_isotopes_as_numbers_rather_than_differences_from_normal)
    to_be_written = iso;
  else
    to_be_written = iso - normal_isotope;

  if (to_be_written > 99)
    output_buffer += digits2[0];
  else if (to_be_written >= 0)
    output_buffer += digits2[to_be_written];
  else if (to_be_written > -10)
    output_buffer << to_be_written;      // just write it
  else       // will be done in the M ISO records
    output_buffer += digits2[0];

  return;
}

/*
  Once reactions needed to be able to write part of an mdl file record,
  this becomes a stand-alone function
*/

int
Molecule::_write_mdl_atom_record_element_charge_and_chirality (atom_number_t i,
                        IWString & output_buffer) const
{
  const Atom * a = _things[i];

  const Element * e = a->element();

  if (! e->is_in_periodic_table() && 'R' == e->symbol()[0] && 2 == e->symbol().length() && isdigit(e->symbol()[1]))
    output_buffer = "R# ";
  else if (mdl_write_aromatic_atoms && NULL != _aromaticity && IS_AROMATIC_ATOM(_aromaticity[i]))
  {
    if (0 == e->aromatic_symbol().length())
      output_buffer = e->symbol();
    else
      output_buffer = e->aromatic_symbol();
    (void) output_buffer.extend(3, ' ');
  }
  else
  {
    output_buffer = e->symbol();
    (void) output_buffer.extend(3, ' ');
  }

  int iso = a->isotope();
  if (0 == iso)
    output_buffer += digits2[0];
  else
    append_isotope_information(output_buffer, iso, e->normal_isotope());
      
  output_buffer += digits3[convert_to_mdl_charge(a->formal_charge())];

  int s;
  if(include_chiral_info_in_mdl_outputs() && _chiral_centres.number_elements())
    s = _mdl_atom_stereo_value(i);
  else
    s = 0;

  if(s)
    output_buffer += digits3[s];
  else if(isis_standard_records)
    output_buffer += digits3[0];

  return 1;
}

/*
  Writing out the atoms and bonds record of an MDL file is complex
*/

int
Molecule::_mdl_write_atoms_and_bonds_record (ostream & os,
                                             int nfc,
                                             int iat,
                                             int isis_standard_records) const
{
  int write_stereo_info = write_mdl_chiral_flags();
  if(write_stereo_info)
    write_stereo_info = _chiral_centres.number_elements();

  if (NULL == digits2)
    initialise_digits();

  int nb = _bond_list.number_elements();

  os << digits3[_number_elements] << digits3[nb];
  if(isis_standard_records)
  {
    os << "  0";     // number of atoms lists
    os << "   ";     // obsolete
    if(_chiral_centres.number_elements() && write_stereo_info)
      os << "  1";
    else
      os << "  0";
    os << "  0";     // number of stext entries
    os << "   ";     // reaction components + 1
    os << "   ";     // number of reactants
    os << "   ";     // number of products
    os << "   ";     // number of intermediates

//  Work out the number of M lines

    int mmm = 0;

    if(nfc)
    {
      if (0 == nfc % 8)
        mmm = nfc / 8;
      else
        mmm = nfc / 8 + 1;
    }

    if(iat)
    {
      if (0 == iat % 8)
        mmm += iat / 8;
      else
        mmm += iat / 8 + 1;
    }

    mmm++;    // don't forget the 'M  END' line

    os << digits3[mmm];

    os << " V2000";  // Ctab version
  }

  os << newline_string();

  return os.good();
}

/*
  This function just writes out the connection table.
*/

/*
  Variable write_stereo_info is really a counter of the number
  of remaining stereo centres to write.
*/

int
Molecule::write_connection_table_mdl (ostream & os) const
{
  assert(ok());
  assert(os.good());

  if (NULL == digits3)
    initialise_digits();

  int nfc = number_formally_charged_atoms();
  int iat = number_isotopic_atoms();

  _mdl_write_atoms_and_bonds_record(os, nfc, iat, isis_standard_records);

  int nb = _bond_list.number_elements();
  int width = os.width();
  IWString output_buffer;
  output_buffer.resize(80);

  int number_r_groups_present = 0;

  for (int i = 0; i < _number_elements && os.good(); i++)
  {
    const Atom * a = _things[i];

    a->write_coordinates(os);

    _write_mdl_atom_record_element_charge_and_chirality(i, output_buffer);

    if(isis_standard_records)
      output_buffer += "  0  0  0           0  0  0";

    output_buffer += newline_string();

    os << output_buffer;

    if (! a->element()->is_in_periodic_table())    // is this an R# group
    {
      const IWString & s = a->element()->symbol();
      if ('R' == s[0] && 2 == s.length() && isdigit(s[1]))
        number_r_groups_present++;
    }
  }

  for (int i = 0; i < nb && os.good(); i++)
  {
    const Bond *b = bondi(i);
    output_buffer = digits3[b->a1() + 1];
    output_buffer += digits3[b->a2() + 1];

    if (mdl_write_aromatic_bonds && (b->is_aromatic() || b->is_permanent_aromatic()))
      output_buffer += digits3[4];
    else if(b->is_single_bond())
      output_buffer += digits3[1];
    else if(b->is_double_bond())
      output_buffer += digits3[2];
    else if(b->is_triple_bond())
      output_buffer += digits3[3];
    else if(b->is_aromatic())
      output_buffer += digits3[4];
    else
      output_buffer += digits3[1];   // defaults to single

    int directionality_written = 1;

    if(b->is_wedge_up())
      output_buffer += digits3[1];
    else if(b->is_wedge_down())
      output_buffer += digits3[6];
    else if(b->is_wedge_either())
      output_buffer += digits3[4];
    else if(b->is_cis_trans_either_double_bond())
      output_buffer += digits3[3];
    else
      directionality_written = 0;

    if(isis_standard_records)
    {
      if (! directionality_written)
        output_buffer << "  0";

      output_buffer += "     0  0";
    }

    output_buffer += newline_string();

    os << output_buffer;

    output_buffer.resize_keep_storage(0);
  }

  os.width(width);

  int need_m_end = 0;

  if (write_mdl_charges_as_m_chg && nfc)
  {
    _write_m_chg_records(os, nfc);
    need_m_end = 1;
  }

  if(iat)
  {
    _write_m_iso_records(os, iat);
    need_m_end = 1;
  }

  if(number_r_groups_present)
    _write_M_RGP_records(os);

  if (isis_standard_records || (need_m_end && write_mdl_m_end_record) || write_mdl_m_end_record > 1)
    os << "M  END" << newline_string();

  if (! os.good())
  {
    cerr << "Molecule::write_connection_table_mdl: cannot write\n";
    return 0;
  }

  return 1;
}

/*
  Write M CHG record(s)

  There can be a maximum of 8 charge entries per record.
*/

#define M_CHG_PER_RECORD 8

int
Molecule::_write_m_chg_records (ostream & os, int nc) const
{
  assert (nc > 0);

  int j = 0;    // index of atom being processed

  while(nc)
  {
    os << "M  CHG ";

    int items_this_record;
    if (nc > M_CHG_PER_RECORD)
      items_this_record = M_CHG_PER_RECORD;
    else
      items_this_record = nc;

    nc -= items_this_record;

    os << setw(2) << items_this_record;
    int items_written = 0;

    while (items_written < items_this_record)
    {
      formal_charge_t q = _things[j++]->formal_charge();
      if(q)
      {
        os << setw(4) << j << setw(4) << q;
        items_written++;
      }
    }

    os << newline_string();
  }

  return os.good();
}

int
Molecule::_write_m_iso_records (ostream & os, int n) const
{
  assert (n > 0);

  int j = 0;        // index of atom being processed

  while(n)
  {
    os << "M  ISO ";

    int items_this_record;
    if (n > M_CHG_PER_RECORD)
      items_this_record = M_CHG_PER_RECORD;
    else
      items_this_record = n;

    n -= items_this_record;

    os << setw(2) << items_this_record;

    int items_written = 0;
    while (items_written < items_this_record)
    {
      const Atom * a = _things[j++];

      if (0 == a->isotope())
        continue;

      int to_be_written;
      if(write_M_isotopes_as_numbers_rather_than_differences_from_normal)
        to_be_written = a->isotope();
      else
        to_be_written = a->isotope() - a->element()->normal_isotope();

      os << setw(4) << j << setw(4) << to_be_written;
      items_written++;
    }

    os << newline_string();
  }

  return os.good();
}

/*
  This function has only been used for mdl files.
  Discern cis-trans centres from examining the depictions.

  The bond is

     a1          a5
       \        /
        a3 -- a4
       /       \
     a2         a6

  In order to get consistent placement, we need to "grow" cis-trans perceptions once we get one
*/

int
Molecule::discern_cis_trans_bonds_from_depiction()
{
  if (_number_elements < 4)
    return 1;

 (void) compute_aromaticity_if_needed();   // path_scoring needs aromaticity

  int nb = _bond_list.number_elements();

  int rc = 0;
  for (int i = 0; i < nb; i++)
  {
    Bond * b = _bond_list[i];

    if (! b->is_double_bond())
      continue;

    atom_number_t a3 = b->a1();
    atom_number_t a4 = b->a2();

    int acon = _things[a3]->ncon();
    if (1 == acon || acon > 3)
      continue;

    acon = _things[a4]->ncon();
    if (1 == acon || acon > 3)
      continue;

    if (6 != _things[a3]->atomic_number())    // for compatability with Daylight, we only do Carbon atoms
      continue;

    if (6 != _things[a4]->atomic_number())
      continue;

    if (unspecified_double_bond_atoms.contains(a3))
      continue;

    if (is_ring_atom(a3) && is_ring_atom(a4))   // if a3 and a4 are both ring atoms, ignore
      continue;

    if(_discern_cis_trans_bond_from_depiction(b))
      rc++;
  }

  return rc;
}

/*
  We need special treatment for the case where there are multiple
  wedge bonds to the same atom.
*/

int
Molecule::_multiple_wedge_bonds_to_atom (atom_number_t zatom) const
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (b->is_wedge_definitive())
      rc++;
  }

  return rc > 1;
}

//#define DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS

int
Molecule::discern_chirality_from_wedge_bonds()
{
#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Discerning chirality from wedge bonds if present\n";
#endif

  _chiral_centres.resize_keep_storage(0);

  int nb = _bond_list.number_elements();

  int rc = 1;     // success until found otherwise

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = _bond_list[i];

    if (! b->is_wedge_definitive())
      continue;

    atom_number_t a1 = b->a1();

    int tmprc;
    if (_multiple_wedge_bonds_to_atom(a1))
      tmprc = _discern_chirality_from_multiple_wedge_bonds(a1);
    else if(b->is_wedge_up())
      tmprc = _discern_chirality_from_wedge_bond(b->a1(), b->a2(), 1);
    else if(b->is_wedge_down())
      tmprc = _discern_chirality_from_wedge_bond(b->a1(), b->a2(), -1);
    else if(b->is_wedge_any())
      tmprc = _create_unspecified_chirality_object(b->a1()); 
    else
      continue;

    if (0 == tmprc)
    {
      if (mdl_display_invalid_chiral_connectivity)
      {
        cerr << "Molecule::_discern_chirality_from_wedge_bonds: cannot determine chirality\n";
        cerr << "Atoms " << b->a1() << " and " << b->a2() << ", molecule '" << _molecule_name << "'\n";
      }
      rc = 0;
    }
  }

  return rc;
}

/*
  If there are 4 atoms connected, we need to determine the directionality from
  the 3 other atoms - those not at the edge of the wedge bond.

  If there are 3 atoms connected, the other connection must be an implicit
  Hydrogen, and the directionality is determined by looking at the 3 connected
  atoms, including the one with the wedge bond
*/

int
Molecule::_discern_chirality_from_wedge_bond (atom_number_t a1, 
                          atom_number_t a2,
                          int direction)
{
  const Atom * aa1 = _things[a1];

#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Bond from atom " << a1 << ' ' << aa1->atomic_symbol() << " to " << a2 << " " << _things[a2]->atomic_symbol() << " direction " << direction << endl;
  cerr << "ncon = " << aa1->ncon() << endl;
#endif

  if (4 == aa1->ncon())
    return _discern_chirality_from_wedge_bond_4(a1, a2, direction);

// We assume that the Hydrogen or lone lair is on the opposite side of the
// page from the atom at the end of the wedge bond, so we reverse the direction

  if (3 == aa1->ncon())
    return _discern_chirality_from_wedge_bond_4(a1, INVALID_ATOM_NUMBER, -direction);

// Connectivity is low, maybe we are just at the wrong end of the bond

  const Atom * aa2 = _things[a2];

  if (4 == aa2->ncon())
    return _discern_chirality_from_wedge_bond_4(a2, a1, -direction);

  if (3 == aa2->ncon())
    return _discern_chirality_from_wedge_bond_4(a2, INVALID_ATOM_NUMBER, direction);

  if (mdl_display_invalid_chiral_connectivity)
    cerr << "Molecule::_discern_chirality_from_wedge_bond: base atom " << aa1->atomic_symbol() << " has " << aa1->ncon() << " connections!\n";

  return return_code_depending_on_ignore_incorrect_chiral_input();
}

/*
  Atom ZATOM is 4 connected and has a wedge bond to atom A2
*/

int
Molecule::_discern_chirality_from_wedge_bond_4 (atom_number_t zatom,
                          atom_number_t a2,
                          int direction)
{
  const Atom * centre = _things[zatom];

// The other 3 atoms connected

  Coordinates c[3];
  atom_number_t a[3];

  int nfound = 0;
  for (int i = 0; i < centre->ncon(); i++)
  {
    atom_number_t j = centre->other(zatom, i);

    if (a2 == j)
      continue;

    a[nfound] = j;

    c[nfound] = *(_things[j]) - *centre;
    c[nfound].normalise();

    nfound++;
  }

  assert (3 == nfound);

// I observed that a positive dot product means anti-clockwise rotation from c[0]

  angle_t theta1 = c[0].angle_between_unit_vectors(c[1]);
  angle_t theta2 = c[0].angle_between_unit_vectors(c[2]);

  if (0.0 == theta1 && 0.0 == theta2)
  {
    cerr << "Molecule::_discern_chirality_from_wedge_bond_4: two zero angles encountered '" << _molecule_name << "'\n";
    return 1;     // ignore the problem
  }

  Coordinates x01(c[0]);
  x01.cross_product(c[1]);
  Coordinates x12(c[1]);
  x12.cross_product(c[2]);
  Coordinates x20(c[2]);
  x20.cross_product(c[0]);

  coord_t z01 = x01.z();
  coord_t z12 = x12.z();
  coord_t z20 = x20.z();

/*
  OK, there are bugs in this, but I don't have the time to chase them down. 
  If there are multiple wedge bonds to an atom, we will get multiple chiral
  centres. This example shows a case where the two chiral centre objects are
  incompatible - file was called t9b.mol


  -ISIS-  10300013162D

  5  4  0  0  0  0  0  0  0  0999 V2000
    1.4417   -3.2708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8292   -7.6583    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    5.8000  -13.8833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.1583   -1.6583    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   11.2000   -4.5250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  1  0  0  0
  1  2  1  0  0  0  0
  2  4  1  6  0  0  0
  2  5  1  0  0  0  0
M  END

  resolve this sometime if it ever matters.
*/

  int rotation;     // clockwise is positive - going 0 -> 1 -> 2

  if (z01 >= 0.0 && z12 >= 0.0)
  {
    rotation = 1;
  }
  else if (z01 < 0.0 && z12 < 0.0)
  {
    rotation = -1;
  }
  else if (z12 >= 0.0 && z20 >= 0.0)
  {
    rotation = 1;
  }
  else if (z12 < 0.0 && z20 < 0.0)
  {
    rotation = -1;
  }
  else if (z20 >= 0.0 && z01 >= 0.0)
  {
    rotation = 1;
  }
  else if (z20 < 0.0 && z01 < 0.0)
  {
    rotation = -1;
  }
  else
  {
    cerr << "Molecule::_discern_chirality_from_wedge_bond_4: unusual geometry z01 = " << z01 << " z12 = " << z12 << " z20 = " << z20 << endl;
    assert (NULL == "this should not happen");
  }

#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  for (int i = 0; i < 3; i++)
  {
    cerr << i << " is atom " << a[i] << ' ' << _things[a[i]]->atomic_symbol() << endl;
  }
  cerr << "z01 = " << z01 << " z12 = " << z12 << " z20 = " << z20 << " rotation " << rotation << endl;
#endif

  if (rotation < 0)
    return _create_chiral_centre(zatom, a[0], a[1], a[2], a2, direction);

  return _create_chiral_centre(zatom, a[0], a[2], a[1], a2, direction);
}

int
Molecule::_create_chiral_centre (atom_number_t zatom,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 atom_number_t a3,
                                 atom_number_t a4,
                                 int direction)
{
#ifdef DEBUG_DISCERN_CHIRALITY_FROM_WEDGE_BONDS
  cerr << "Molecule::_create_chiral_centre: " << zatom << " a1 = " << a1 << " a2 = " << a2 << " a3 = " << a3 << " a3 = " << a3 << " a4 = " << a4 << " direction " << direction << endl;
#endif

  Chiral_Centre * c = new Chiral_Centre(zatom);
  c->set_top_front(a1);
  c->set_chirality_known(1);

  if (INVALID_ATOM_NUMBER == a4)
  {
    int lp;

    if (1 == _things[zatom]->implicit_hydrogens())
      c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);
    else if (_things[zatom]->lone_pair_count(lp) && 1 == lp)
      c->set_top_back(CHIRAL_CONNECTION_IS_LONE_PAIR);
    else
    {
      delete c;
      if (mdl_display_invalid_chiral_connectivity)
        cerr << "Molecule::_create_chiral_centre: atom " << zatom << " strange chemistry, '" << _molecule_name << "'\n";
      return return_code_depending_on_ignore_incorrect_chiral_input();
    }
  }
  else
    c->set_top_back(a4);

  if (direction > 0)
  {
    c->set_left_down(a2);
    c->set_right_down(a3);
  }
  else
  {
    c->set_left_down(a3);
    c->set_right_down(a2);
  }

  _chiral_centres.add(c);

  return 1;
}

int
Molecule::_create_unspecified_chirality_object (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  Chiral_Centre * c = new Chiral_Centre(zatom);

  if (3 == a->ncon())
    c->set_top_back(CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN);

  for (int i = 0; i < a->ncon(); i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (0 == i)
      c->set_top_front(j);
    else if (1 == i)
      c->set_left_down(j);
    else if (2 == i)
      c->set_right_down(j);
    else 
      c->set_top_back(j);
  }

  c->set_chirality_known(0);

  _chiral_centres.add(c);

  return 1;
}

/*
  We have the case where there are multiple wedge bonds going to a given atom.
*/

int
Molecule::_discern_chirality_from_multiple_wedge_bonds (atom_number_t zatom)
{
  const Atom * a = _things[zatom];

  int acon = a->ncon();

  const Bond * wb1 = NULL;
  const Bond * wb2 = NULL;

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_wedge_definitive())
      continue;

    if (NULL == wb1)
      wb1 = b;
    else if (NULL == wb2)
      wb2 = b;
    else if (mdl_display_invalid_chiral_connectivity)
      cerr << "Molecule::_discern_chirality_from_multiple_wedge_bonds:warning, too many wedge bonds to atom " << zatom << ", ignored\n";
  }

  assert (NULL != wb2);

// Once upon a time I had a check here trying to figure out if the directionalities are
// consistent, but I'm not sure there are any inconsistencies possible with just two
// bonds being processed - see above

  if (wb1->is_wedge_up())
    return _discern_chirality_from_wedge_bond(wb1->a1(), wb1->a2(), 1);
  else
    return _discern_chirality_from_wedge_bond(wb1->a1(), wb1->a2(), -1);
}

int
Molecule::_assign_strange_atomic_symbol_to_atom (atom_number_t zatom,
                                                 const_IWSubstring s)
{
  const Element * e = get_element_from_symbol_no_case_conversion(s);

  if (NULL != e)
    ;
  else if (! auto_create_new_elements())
  {
    cerr << "Molecule::_assign_strange_atomic_symbol_to_atom:autocreate elements needed for R# and G group\n";
    return 0;
  }
  else
    e = create_element_with_symbol(s);

  if (NULL == e)
  {
    cerr << "Molecule::_assign_strange_atomic_symbol_to_atom:cannot create element '" << s << "'\n";
    return 0;
  }

// cerr << "Created element '" << e->symbol() << "'\n";

  _things[zatom]->set_element(e);

  return 1;
}

int
Molecule::_parse_M_RGP_record (const const_IWSubstring & buffer)
{
  Aprop atom_properties[MAX_PAIRS];

  int npairs;
  if (! fill_atom_property_array(buffer, npairs, atom_properties))
    return 0;

  for (int i = 0; i < npairs; i++)
  {
    const Aprop & api = atom_properties[i];

    atom_number_t a = api._atom_number - 1;
    int rnumber = api._property;

    IWString tmp;
    tmp << "R" << rnumber;

    if (! _assign_strange_atomic_symbol_to_atom(a, tmp))
      return 0;
  }

  _set_modified();

  return 1;
}

int
Molecule::_write_M_RGP_records (ostream & os) const
{
  for (int i = 0; i < _number_elements; i++)
  {
    const Element * e = _things[i]->element();

    if(e->is_in_periodic_table())
      continue;

    const IWString & s = e->symbol();

    if ('R' == s[0] && isdigit(s[1]))
      ;
    else
      continue;

    os << "M  RGP  1  " << setw(3) << (i + 1) << "   " << s[1] << newline_string();
  }

  return os.good();
}

int
Molecule::_set_elements_based_on_atom_aliases (const resizable_array_p<Atom_Alias> & a)
{
  set_atomic_symbols_can_have_arbitrary_length(1);
  set_auto_create_new_elements(1);

  for (int i = 0; i < a.number_elements(); i++)
  {
    const Atom_Alias * ai = a[i];

    atom_number_t j = ai->atom_number();

    assert (OK_ATOM_NUMBER (this, j));

    const IWString alias = ai->alias();

    int isotope_not_used;
    const Element * e = get_element_from_symbol(alias, isotope_not_used);

    if (NULL == e)
      e = create_element_with_symbol(alias);

    assert (NULL != e);

    _things[j]->set_element(e);
  }

  return 1;
}

int
Molecule::_process_mdl_g_record (const IWString & g,
                                 const const_IWSubstring & buffer)
{
  if (0 == buffer.length())
  {
    cerr << "Molecule::_process_mdl_g_record:empty G record\n";
    return 0;
  }

  assert (g.starts_with("G "));

  if (3 != g.nwords())
  {
    cerr << "Molecule::_process_mdl_g_record:Invalid G record\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;
  g.nextword(token, i);

  g.nextword(token, i);

  atom_number_t a;
  if (! token.numeric_value(a) || a < 1 || a > _number_elements)
  {
    cerr << "Molecule::_process_mdl_g_record:invalid atom number in G record\n";
    return 0;
  }

  return _assign_strange_atomic_symbol_to_atom(a - 1, buffer);
}

void
reset_mdl_file_scope_variables()
{
  isis_standard_records = 0;
  ignore_unrecognised_m_records = 0;
  report_unrecognised_records = 1;
  die_on_erroneous_m_input = 1;
  write_mdl_dollar_signs = 1;
  write_mdl_m_end_record = 1;
  write_v30_mdl_files = 0;
  _ignore_self_bonds = 0;
  _write_mdl_chiral_flags = 1;
  _include_chiral_info_in_mdl_outputs = 1;
  _mdl_read_h_correct_chiral_centres = 1;
  _mdl_write_h_correct_chiral_centres = 1;
  extract_isis_extregno = 0;
  fetch_all_sdf_identifiers = 0;
  take_first_tag_as_name = 0;
  prepend_sdfid = 1;
  discard_sdf_molecule_name = 0;
  _multi_record_tag_data_present = 1;
  mdl_write_aromatic_bonds = 0;
  mdl_write_aromatic_atoms = 0;
  accumulate_mdl_chirality_features = 0;
  display_non_organic_chirality_messages = 1;
  mdl_display_invalid_chiral_connectivity = 1;
  truncate_long_symbols = 0;
  discern_chirality_from_wedge_bonds = 0;
  write_isotopes_as_numbers_rather_than_differences_from_normal = 0;
  write_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;
  read_isotopes_as_numbers_rather_than_differences_from_normal = 0;
  read_M_isotopes_as_numbers_rather_than_differences_from_normal = 1;
  allow_deuterium = 0;
  allow_tritium = 0;
  if (NULL != input_bond_type_translation_table)
  {
    delete [] input_bond_type_translation_table;
    input_bond_type_translation_table = NULL;
  }

  name_in_m_tag.resize(0);
  unspecified_chiral_atoms_last_molecule.resize(0);
  unspecified_double_bond_atoms.resize(0);

  up_bonds_last_molecule.resize(0);
  down_bonds_last_molecule.resize(0);
  squiggle_bonds_last_molecule.resize(0);

  change_long_symbols_to.resize(0);

  mdl_g_records_hold_atom_symbols = 0;
  a_records_found = 0;
  g_records_found = 0;
  _set_elements_based_on_atom_aliases = 0;
  write_mdl_charges_as_m_chg = 1;

  aliases.resize(0);

  rdfile_identifiers.resize(0);
  rdfile_start_of_record.resize(0);

  return;
}
// arch-tag: 5295edde-b728-42a2-8d20-d284245b824e
