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
#include <stdlib.h>

#include "misc.h"
#include "iw_auto_array.h"

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

#define COMPILING_MDL_CC

#include "mdl.h"
#include "substructure.h"

#include "molecule_to_query.h"

#define COMPILING_MDL_CC

#include "mdl_molecule.h"
#include "mdl_atom_record.h"
#include "molecule_to_query.h"

static int convert_a_and_q_atoms_to_atom_lists = 0;

void
set_convert_a_and_q_atoms_to_atom_lists (int s)
{
  convert_a_and_q_atoms_to_atom_lists = s;
}

static int convert_not_atom_lists_to_organic_lists = 0;

void
set_convert_not_atom_lists_to_organic_lists (int s)
{
  convert_not_atom_lists_to_organic_lists = s;
}

MDL_Molecule::MDL_Molecule()
{

#ifdef USE_IWMALLOC
  cerr << "Checking MDL_Molecule\n";
  iwmalloc_check_all_malloced (stderr);
#endif

  return;
}

MDL_Molecule::~MDL_Molecule()
{
  return;
}

int
MDL_Molecule::_parse_atom_alias(iwstring_data_source & input,
                                const const_IWSubstring & buffer)
{
  assert (buffer.starts_with("A  "));

  const_IWSubstring tmp (buffer);

  tmp.remove_leading_chars (3);

  tmp.strip_leading_blanks();

  atom_number_t zatom;

  if (! tmp.numeric_value (zatom) || zatom < 1)
  {
    cerr << "MDL_Molecule::_parse_atom_alias: invalid atom number specification '" << buffer << "'\n";
    return 0;
  }

  zatom--;

  if (! input.next_record(tmp))
  {
    cerr << "MDL_Molecule::_parse_atom_alias:eof\n";
    return 0;
  }

  _mdl_atom[zatom]->set_alias (tmp);

  return 1;
}

int
MDL_Molecule::_parse_M_record (iwstring_data_source & input,
                               const const_IWSubstring & buffer,
                               ::resizable_array_p<ISIS_Link_Atom> & ltmp,
                               int & fatal)
{
  if (buffer.starts_with ("A  "))
  {
    if (! _parse_atom_alias (input, buffer))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with ("M  ALS"))
  {
    if (! _parse_atom_list (buffer))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with ("M  LIN"))
  {
    if (! _parse_link_record (buffer, ltmp))
    {
      fatal = 1;
      return 0;
    }

    return 1;;
  }

  if (buffer.starts_with ("M  CHG"))
  {
    Aprop atom_properties[MAX_PAIRS];

    int tokens;
    if (! fill_atom_property_array (buffer, tokens, atom_properties))
      return 0;

    if (0 == tokens)
      ;
    else if (! mdl_add_m_formal_charge (tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }
  
  if (buffer.starts_with ("M  ISO"))
  {
    Aprop atom_properties[MAX_PAIRS];

    int tokens;

    if (! fill_atom_property_array (buffer, tokens, atom_properties))
      return 0;

    if (0 == tokens)
      ;
    else if (! mdl_add_m_isotope (tokens, atom_properties))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with ("M  UNS"))
  {
    Aprop atom_properties[MAX_PAIRS];

    int tokens;
    if (! fill_atom_property_array (buffer, tokens, atom_properties))
      return 0;

    if (! _set_unsaturation_specifications (atom_properties, tokens))
    {
      fatal = 1;
      return 0;
    }
    return 1;
  }

  if (buffer.starts_with ("M  SUB"))
  {
    Aprop atom_properties[MAX_PAIRS];

    int tokens;
    if (! fill_atom_property_array (buffer, tokens, atom_properties))
      return 0;

    if (! _set_substitution_specifications (atom_properties, tokens))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }

  if (buffer.starts_with ("M  RBC"))
  {
    Aprop atom_properties[MAX_PAIRS];

    int tokens;
    if (! fill_atom_property_array (buffer, tokens, atom_properties))
      return 0;

    if (! _set_ring_bond_specifications (atom_properties, tokens))
    {
      fatal = 1;
      return 0;
    }

    return 1;
  }


  fatal = 0;
  return 0;    // not recognised
}

int
MDL_Molecule::_parse_atom_list (const IWString & buffer)
{
  assert (buffer.starts_with("M  ALS"));
  
  const_IWSubstring token;
  if (! buffer.word(2, token))
  {
    cerr << "MDL_Molecule::_parse_atom_list:invalid record '" << buffer << "'\n";
    return 0;
  }

  atom_number_t zatom;
  if (! token.numeric_value(zatom) || zatom < 1)
  {
    cerr << "MDL_Molecule::_parse_atom_list:invalid atom number '" << buffer << "'\n";
    return 0;
  }

  zatom--;

  if (! _mdl_atom.ok_index(zatom))
  {
    cerr << "MDL_Molecule::_parse_atom_list:invalid atom number in '" << buffer << "', natoms = " << _mdl_atom.number_elements() << endl;
    return 0;
  }

  return _mdl_atom[zatom]->build_atom_list(buffer);
}

int
MDL_Molecule::_parse_link_record (const IWString & buffer,
                                  ::resizable_array_p<ISIS_Link_Atom> & ltmp)
{
  assert (buffer.starts_with ("M  LIN"));

  Link_Atom * l = new Link_Atom;

  atom_number_t a;

  if (! l->initialise_from_mdl_record (buffer, Molecule::natoms(), a))
  {
    cerr << "MDL_Molecule::_parse_link_record:invalid link record '" << buffer << "'\n";
    delete l;
    return 0;
  }

  if (2 != Molecule::ncon(a))
  {
    cerr << "MDL_Molecule::_parse_link_record:link atoms must have 2 connections\n";
    return 0;
  }

// The bond types either side of the link atom must be the same

  const MDL_Bond_Data * b1 = mdl_bond_between_atoms (a, l->a1());
  const MDL_Bond_Data * b2 = mdl_bond_between_atoms (a, l->a2());

  if (b1->btype() != b2->btype())
  {
    cerr << "MDL_Molecule::_parse_link_record:inconsitent bond types " << b1->btype() << " vs " << b2->btype() << '\n';
    return 0;
  }

  if (b1->bond_topology() != b2->bond_topology())
  {
    cerr << "MDL_Molecule::_parse_link_record:inconsitent bond topologies " << b1->bond_topology() << " vs " << b2->bond_topology() << '\n';
    return 0;
  }

//cerr << "MDL_Molecule::_parse_link_record:set link atom bond type " << b1->btype() << ", bond topology " << b1->bond_topology() << endl;

  l->set_bond_type(b1->btype());
  l->set_bond_topology(b1->bond_topology());

  _link_atom.add(l);

  _mdl_atom[l->a1()]->increment_connections_lost();
  _mdl_atom[l->a2()]->increment_connections_lost();

  const IWString & s = atomi(a)->atomic_symbol();

  if (! l->set_symbol(s))
  {
    cerr << "MDL_Molecule::_parse_link_record:invalid element??? '" << s << "'\n";
    return 0;
  }

  l->set_mdl_atom_data (_mdl_atom[a]);

  remove_atom(a);    // the explicit link atom is no longer needed

  return 1;
}

/*
  Removing an atom is hard because of the arrays we must keep in sync
*/

int
MDL_Molecule::remove_atom(atom_number_t zatom)
{
  Set_of_Atoms check_ncon;

  atomic_number_t zremove = atomic_number(zatom);

  for (int i = Molecule::nedges() - 1; i >= 0; i--)
  {
    const Bond * b = Molecule::bondi(i);

    if (! b->involves(zatom))
      continue;

    check_ncon.add(b->other(zatom));

    Molecule::remove_bond(i);

    _mdl_bond.remove_item(i);
  }

  for (int i = 0; i < _link_atom.number_elements(); i++)
  {
    _link_atom[i]->atom_is_being_removed(zatom);
  }

  for (int i = 0; i < check_ncon.number_elements(); i++)
  {
    atom_number_t j = check_ncon[i];

    _mdl_atom[j]->connected_atom_is_being_removed(zremove);
  }

  _mdl_atom.remove_item(zatom);

  Molecule::remove_atom(zatom);

//cerr << "After removing atom " << zatom << " '" << smiles() << "'\n";

  return 1;
}

int
MDL_Molecule::remove_atoms (const int * to_remove)
{
  int rc = 0;

  for (int i = Molecule::natoms() - 1; i >= 0; i--)
  {
    if (to_remove[i])
    {
      MDL_Molecule::remove_atom(i);
      rc++;
    }
  }

  return 1;
}

/*
  If we are reading an MDL type file, we can process it here. Otherwise
  go to the default processing
*/

int
MDL_Molecule::read_molecule_ds(iwstring_data_source & input,
                               int input_type)
{
//cerr << "MDL_Molecule::read_molecule_ds:type " << input_type << endl;

  if (SDF == input_type)
    ;
  else if (MDL == input_type)
    ;
  else
    return Molecule::read_molecule_ds(input, input_type);

  return MDL_Molecule::read_molecule_mdl_ds(input);
}

static const Element * element_a = NULL;
static const Element * element_q = NULL;
static const Element * element_l = NULL;

void
do_create_special_elements_for_mdl_stuff()
{
  int isave = auto_create_new_elements();

  if (0 == isave)
    set_auto_create_new_elements(1);

  element_a = get_element_from_symbol_no_case_conversion ("A");
  if (NULL == element_a)
    element_a = create_element_with_symbol("A");

  element_q = get_element_from_symbol_no_case_conversion ("Q");
  if (NULL == element_q)
    element_q = create_element_with_symbol("Q");

  element_l = get_element_from_symbol_no_case_conversion ("L");
  if (NULL == element_l)
    element_l = create_element_with_symbol("L");

  set_auto_create_new_elements(isave);

  return;
}

/*
  Very dismaying that there isn't code sharing between mdl.cc and
  rxnfile.cc, but just too hard to figure out how to do it in a robust
  and maintainable form!
*/

int
MDL_Molecule::read_molecule_mdl_ds (iwstring_data_source & input,
                                    int return_on_m_end)
{
#ifdef USE_IWMALLOC
  cerr << "MDL_Molecule at " << this << endl;

  iwmalloc_check_all_malloced (stderr);
#endif

  if (NULL == element_a)
    do_create_special_elements_for_mdl_stuff();

  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    return 0;
  }

  set_name(buffer);
//cerr << "Name is '" << buffer << "'\n";

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::read_molecule_ds:eof reading header\n";
    return 0;
  }

//cerr << "2nd record '" << buffer << "'\n";

  if (! input.next_record(_third_line_of_input_sdf_file))
  {
    cerr << "MDL_Molecule::read_molecule_ds:eof reading header\n";
    return 0;
  }

//cerr << "Third line " << _third_line_of_input_sdf_file << "'\n";

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::read_molecule_ds:eof reading header\n";
    return 0;
  }

//cerr << "NA NB '" << buffer << "'\n";

  int na, nb;

  if (2 != int3d (buffer, na, nb))
  {
    cerr << "MDL_Molecule::read_molecule_ds: error from int3d '" << buffer << "'\n";
    cerr << "Line " << input.lines_read() << endl;
    return 0;
  }

  if (0 == na && 0 == nb && buffer.contains("V3000"))
    return _read_v3000 (input);

  assert (na >= 0 && (nb >= 0));

  if (! allocate_arrays(na, nb))
    return 0;

  MDL_Atom_Record mdlar;

  for (int i = 0; i < na; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "MDL_Molecule::read_molecule_ds:premature eof\n";
      return 0;
    }

    if (! mdlar.build(buffer))
    {
      cerr << "MDL_Molecule::read_molecule_ds:invalid record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    const_IWSubstring s = mdlar.atomic_symbol();

    Atom * a = mdlar.create_atom();

    if (NULL == a)
    {
      cerr << "MDL_Molecule::read_molecule_ds:cannot create atom, line " << input.lines_read() << endl;
      return 0;
    }

    Molecule::add(a);

    _mdl_atom[i]->extract_info_from_mdl_file_record(mdlar);

    if (0 == mdlar.astere())     // the most common case, no chirality
      ;
    else if (! _mdl_atom_is_chiral_centre (Molecule::natoms() - 1, mdlar.astere()))
    {
      cerr << "MDL_Molecule::read_molecule_ds:invalid chirality on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  int wedge_bonds_present = 0;

  MDL_Bond_Record mdlbr;

  for (int i = 0; i < nb; i++)
  {
    if (! input.next_record(buffer))
    {
      cerr << "MDL_Molecule::read_molecule_ds:premature eof\n";
      return 0;
    }

    if (! mdlbr.build(buffer, na))
    {
      cerr << "MDL_Molecule::read_molecule_ds:invalid bond record, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    bond_type_t bt_for_molecule;
    if (! mdlbr.bond_type_for_molecule(bt_for_molecule))
    {
      cerr << "MDL_Molecule::read_molecule_ds:invalid bond type\n";
      return 0;
    }

    if (! _mdl_bond[i]->extract_info_from_mdl_file_record(mdlbr))
      return 0;

    Molecule::add_bond(mdlbr.a1(), mdlbr.a2(), bt_for_molecule, 1);

    if (mdlbr.directionality())
    {
      _mdl_set_bond_directionality(mdlbr.a1(), mdlbr.a2(), mdlbr.directionality());
      wedge_bonds_present++;
    }
  }

  if (nb > 0)
    check_bonding();

  if (ignore_all_chiral_information_on_input())
    Molecule::remove_all_chiral_centres();
  else if (wedge_bonds_present)
    discern_chirality_from_wedge_bonds();

  if (0 == chiral_centres())    // none to worry about
    ;
  else if (_complete_chiral_centres_from_mdl_files())    // good
    ;
  else               // OOPS, bad chirality info
  {
    cerr << "MDL_Molecule::read_molecule_mdl_ds: erroneous chiral input '" << Molecule::name() << "'\n";
    Molecule::remove_all_chiral_centres();
    if (! ignore_incorrect_chiral_input())
      return 0;
  }

  int got_mend = 0;
  int got_dollars = 0;

  ::resizable_array_p<ISIS_Link_Atom> ltmp;

  while (input.next_record(buffer))
  {
    if ("$$$$" == buffer)
    {
      got_dollars = 1;
      break;
    }

    if ("M  END" == buffer)
    {
      got_mend = 1;
      if (return_on_m_end)
        break;
      else
        continue;
    }

    int fatal;
    if (Molecule::_common_parse_M_record(buffer, fatal))
      continue;
    else if (fatal)
      return 0;

    if (_parse_M_record(input, buffer, ltmp, fatal))
      continue;
    else if (fatal)
    {
      cerr << "MDL_Molecule::read_molecule_ds:fatal error, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }

    if (read_extra_text_info())
      Molecule::add_extra_text_info(buffer);
  }

  na = natoms();    // link atoms remove the atom, so NA may have changed. Recompute

  for (int i = 0; i < na; i++)
  {
    MDL_Atom_Data * mdla = _mdl_atom[i];

    if (Molecule::is_aromatic(i))
      mdla->set_aromatic(1);
    else
      mdla->set_aromatic(0);
  }

  if (convert_a_and_q_atoms_to_atom_lists)
  {
    for (int i = 0; i < na; i++)
    {
      MDL_Atom_Data * mdla = _mdl_atom[i];

      const IWString & s = atomic_symbol(i);

      mdla->convert_a_or_q_atoms_to_atom_list(s);
    }
  }

  if (convert_not_atom_lists_to_organic_lists)
  {
    for (int i = 0; i < na; i++)
    {
      _mdl_atom[i]->convert_not_atom_lists_to_organic_lists();
    }
  }

#ifdef ECHO_ATOM_LIST_INFO
  cerr << "After reading\n";
  for (int i = 0; i < na; i++)
  {
    cerr << "Atom " << i << " atom list contains " << _mdl_atom[i]->atom_list().number_elements() << " elements\n";
  }
#endif

// For now, I'm ignoring the trailing records
//reset_mdl_molecule_file_scope_variables();
  return na;
}

const ISIS_Atom_List *
MDL_Molecule::atom_list_for_atom (atom_number_t a) const
{
  assert (ok_atom_number (a));

  const ISIS_Atom_List & rc = _mdl_atom[a]->atom_list();

  if (! rc.active())
    return NULL;

  return &rc;
}

int
MDL_Molecule::_common_set_from_aprop (const Aprop * atom_properties,
                                      int ntokens,
                                      int * dest)
{
  for (int i = 0; i < ntokens; i++)
  {
    const Aprop & a = atom_properties[i];

    dest[a._atom_number - 1] = a._property;
  }

  return 1;
}

int
MDL_Molecule::_set_unsaturation_specifications (const Aprop * atom_properties,
                                                          int tokens)
{
  for (int i = 0; i < tokens; i++)
  {
    const Aprop & a = atom_properties[i];

    atom_number_t n = a._atom_number - 1;

    _mdl_atom[n]->set_unsaturation(a._property);
  }

  return 1;
}

int
MDL_Molecule::_set_substitution_specifications (const Aprop * atom_properties,
                                                          int tokens)
{
  for (int i = 0; i < tokens; i++)
  {
    const Aprop & a = atom_properties[i];

    atom_number_t n = a._atom_number - 1;

    _mdl_atom[n]->set_substitution(a._property);
  }

  return 1;
}

int
MDL_Molecule::_set_ring_bond_specifications (const Aprop * atom_properties,
                                                          int tokens)
{
  for (int i = 0; i < tokens; i++)
  {
    const Aprop & a = atom_properties[i];

    atom_number_t n = a._atom_number - 1;

    _mdl_atom[n]->set_ring_bond(a._property);
  }

  return 1;
}







/*
  When we remove explicit hydrogens, we can increment the min_hcount value
  of the atoms to which the H's were attached
*/

int
MDL_Molecule::remove_explicit_hydrogens(atomic_number_t z)
{
  int matoms = natoms();

  if (! arrays_allocated())
    build(*this);

  int number_hydrogens = Molecule::natoms(z);

  if (0 == number_hydrogens)
    return 0;

  int rc = 0;

  for (int i = matoms - 1; i >= 0; i--)
  {
    const MDL_Atom_Data * mdlad = _mdl_atom[i];

    const ISIS_Atom_List & l = mdlad->atom_list();    // don't delete atom lists that start with H

//  cerr << "Atom list on atom " << i << " contains " << l.number_elements() << " items\n";
    if (l.number_elements() > 0)
      continue;

    const Atom * ai = Molecule::atomi(i);

    if (z != ai->atomic_number())
      continue;

    if (1 != ai->ncon())   // I could actually handle the other cases, just messy and not needed
      continue;

    atom_number_t o = ai->other(i, 0);

    if (1 == z)
    {
      _mdl_atom[o]->increment_min_hcount();
      _mdl_atom[o]->increment_explicit_hydrogen_atoms_removed();
    }
    else
      _mdl_atom[o]->increment_connections_lost();

    MDL_Molecule::remove_atom(i);
    rc++;

    if (rc == number_hydrogens)
      break;
  }

  return rc;
}

int
MDL_Molecule::not_atom_lists_present() const
{
  int rc = 0;

  int matoms = Molecule::natoms();

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * madi = MDL_File_Data::mdl_atom_data(i);

    const ISIS_Atom_List & iali = madi->atom_list();

    if (! iali.active())
      continue;

    if (! iali.normal_list())
      rc++;
  }

  return rc;
}

MDL_Molecule::MDL_Molecule(const Molecule & m) : Molecule(m)
{
  MDL_File_Data::build(*this);

  return;
}

MDL_Molecule::MDL_Molecule (const MDL_Molecule & m) : Molecule (m),
                                                MDL_File_Data (m)
{
  return;
}

/*
  Someone has a molecule with R groups. The atoms to which these
  are attached are the only attachment points possible.
*/

int
MDL_Molecule::change_R_groups_to_substitutions (Element_Matcher & rgroup)
{
  int matoms = natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = atomi (i);

    const Element * e = ai->element();

    _mdl_atom[i]->set_substitution(ai->ncon());

    if (! rgroup.matches(e))
      continue;

    if (1 != ai->ncon())   // possible limitation, hmmmm
      continue;

    atom_number_t j = ai->other(i, 0);

//  cerr << "Processing R group, bonded to atom " << j << endl;

    _mdl_atom[j]->set_min_ncon(ncon(j));
    _mdl_atom[j]->set_substitution(0);

    _substitution_points.add(j);

    _substitution_points.adjust_for_loss_of_atom(i);

    MDL_Molecule::remove_atom(i);
    i--;
    matoms--;
    rc++;
  }

  if (0 == rc)
    cerr << "None of the rgroup element matches matched '" << name() << "'\n";

  return rc;
}

int
MDL_Molecule::only_allow_substitutions_at_isotopic_atoms(const Molecule_to_Query_Specifications & mqs)
{
  int matoms = natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = atomi(i);

    int iso = ai->isotope();

    MDL_Atom_Data * mad = _mdl_atom[i];

#ifdef DEBUG_ONLY_ALLOW_SUBSTITUTIONS_AT_ISOTOPIC_ATOMS
    cerr << "Atom " << i << " is isotope " << iso << ", currently " << mad->substitution() << endl;
#endif

    if (0 == iso)
    {
      mad->set_substitution(ai->ncon());
      continue;
    }

    if (hcount(i))
      ;
    else if (16 == atomic_number(i))   // can have any valence it wants
      ;
    else
      cerr << "MDL_Molecule::only_allow_substitutions_at_isotopic_atoms:no open valence in '" << name() << "'\n";

    if (must_have_substituent_at_every_isotopic_atom())
    {
      if (0 != mad->substitution())   // already set
        ;
      else if (isotope_count_means_extra_connections())
        mad->set_substitution(ai->ncon() + iso);
      else
        mad->set_min_ncon(ai->ncon() + 1);
    }
    else
    {
      if (isotope_count_means_extra_connections())
        mad->set_min_ncon(ai->ncon() + iso);
      else
        mad->set_min_ncon(ai->ncon());

      mad->set_substitution(0);
    }

    if (mqs.environment_near_substitution_points_specified())  // do I need this test, or should it always be done
      _substitution_points.add (i);

    set_isotope(i, 0);

    rc++;
  }

#ifdef DEBUG_ONLY_ALLOW_SUBSTITUTIONS_AT_ISOTOPIC_ATOMS
  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * ma = _mdl_atom[i];
    cerr << " atom " << i << " substitution " << ma->substitution() << endl;
  }
#endif

  return rc;
}

int
MDL_Molecule::only_allow_substitutions_at_non_isotopic_atoms()
{
  int matoms = natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = atomi(i);

    int iso = ai->isotope();

    MDL_Atom_Data * mad = _mdl_atom[i];

#ifdef DEBUG_ONLY_ALLOW_SUBSTITUTIONS_AT_ISOTOPIC_ATOMS
    cerr << "Atom " << i << " is isotope " << iso << ", currently " << mad->substitution() << endl;
#endif

    if (0 != iso)   // labelled, set max ncon
    {
      mad->set_substitution(ai->ncon());
      set_isotope(i, 0);
      continue;
    }

//  unlabelled, 

    if (hcount(i))
      ;
    else if (16 == atomic_number(i))   // can have any valence it wants
      ;
    else
      cerr << "MDL_Molecule::only_allow_substitutions_at_non isotopic_atoms:no open valence in '" << name() << "'\n";

    rc++;
  }

#ifdef DEBUG_ONLY_ALLOW_SUBSTITUTIONS_AT_ISOTOPIC_ATOMS
  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * ma = _mdl_atom[i];
    cerr << " atom " << i << " substitution " << ma->substitution() << endl;
  }
#endif

  return rc;
}


int
MDL_Molecule::determine_attachment_points_by_query(Molecule_to_Query_Specifications & mqs)
{
  Substructure_Results sresults;
  int nhits = mqs.substitutions_only_at().substructure_search (*this, sresults);

//cerr << "MDL_Molecule::determine_attachment_points_by_query:nhits = " << nhits << endl;

  if (0 == nhits)
  {
    cerr << "MDL_Molecule::determine_attachment_points_by_query:no hits to substitutions_only_at query\n";
    return 0;
  }

  int matoms = natoms();

  for (int i = 0; i < matoms; i++)
  {
    _mdl_atom[i]->set_substitution(ncon(i));
  }

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding (i);
    
    for (int j = 0; j < e->number_elements(); j++)
    {
      atom_number_t k = e->item(j);

      _mdl_atom[k]->set_substitution(0);

      if (0 == _mdl_atom[k]->min_ncon())
        _mdl_atom[k]->set_min_ncon(ncon(k));
    }

//  should we also set the unsaturation flag?????

    if (mqs.environment_near_substitution_points_specified())
      _substitution_points.add_non_duplicated_elements (*e);
  }

  return 1;
}

int
MDL_Molecule::swap_atoms (atom_number_t a1,
                          atom_number_t a2)
{
  Molecule::swap_atoms(a1, a2);

  return MDL_File_Data::swap_atoms(a1, a2);
}

const MDL_Bond_Data *
MDL_Molecule::mdl_bond_between_atoms (atom_number_t a1,
                                      atom_number_t a2) const
{
  int j = Molecule::which_bond(a1, a2);

  if (j < 0)
  {
    cerr << "MDL_Molecule::mdl_bond_between_atoms:no bond between " << a1 << " and " << a2 << endl;
    return NULL;
  }

  return _mdl_bond[j];
}

int
MDL_Molecule::add (const Element * e)
{
  if (! Molecule::add(e))
    return 0;

  if (! arrays_allocated())
    return 1;

  MDL_Atom_Data * t = new MDL_Atom_Data();
  _mdl_atom.add(t);

  return 1;
}

/*
  The atomic symbol is something other than a known atomic symbol
*/

int
MDL_Molecule::add_atom_based_on_symbol (const IWString & s)
{
  set_auto_create_new_elements(1);
  if (! s.starts_with('[') && s.length() <= 2)
  {
    const Element * e = get_element_from_symbol_no_case_conversion(s);
    return MDL_Molecule::add(e);
  }
  else   // hope that somethingin the MDL_Atom_Data specifies this
  {
    const Element * e = get_element_from_symbol_no_case_conversion("*");
    return MDL_Molecule::add(e);
  }
}

int
MDL_Molecule::add_bond (atom_number_t a1, atom_number_t a2, 
                        bond_type_t bond_for_molecule,
                        bond_type_t query_bond,
                        int bond_topology)
{
  if (! Molecule::add_bond (a1, a2, bond_for_molecule))
    return 0;

  if (! arrays_allocated())
    return 1;

  MDL_Bond_Data * b = new MDL_Bond_Data();

  b->set_btype(query_bond);

  if (bond_topology >= 0)
    b->set_bond_topology(bond_topology);

  _mdl_bond.add(b);

  return 1;
}

int
MDL_Molecule::remove_bond_between_atoms (atom_number_t a1, atom_number_t a2)
{
  int i = Molecule::which_bond(a1, a2);

  if (i < 0)
  {
    cerr << "MDL_Molecule::remove_bond_between_atoms:atoms " << a1 << " and " << a2 << " not bonded\n";
    return 0;
  }

  if (! Molecule::remove_bond_between_atoms(a1, a2))   // how could that fail?
    return 0;

  _mdl_bond.remove_item(i);

// We don't have any link atoms, this function was built for help enumerating link atoms

  return 1;
}

int
MDL_Molecule::set_mdl_atom_data (int i, const MDL_Atom_Data * a)
{
  if (! arrays_allocated())
    build(*this);

  assert (_mdl_atom.ok_index(i));

  *(_mdl_atom[i]) = *a;

  return 1;
}

int
MDL_Molecule::_read_v3000 (iwstring_data_source & input)
{
  IWString buffer;

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:no data\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 BEGIN CTAB"))
  {
    cerr << "MDL_Molecule::_read_v3000:not V30 BEGIN CTAB '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:no counts\n";
    return 0;
  }

// M  V30 COUNTS 4 3 0 0 0

  if (! buffer.starts_with("M  V30 COUNTS ") || buffer.nwords() < 5)
  {
    cerr << "MDL_Molecule::_read_v3000:not counts '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring token = buffer.word(3);
  int na;
  if (! token.numeric_value(na) || na < 1)
  {
    cerr << "MDL_Molecule::_read_v3000:invalid atom count '" << buffer << "'\n";
    return 0;
  }

  token = buffer.word(4);
  int nb;
  if (! token.numeric_value(nb) || nb < 0)
  {
    cerr << "MDL_Molecule::_read_v3000:invalid bond count '" << buffer << "'\n";
    return 0;
  }

  _fill_empty_molecule_with_null_atoms (na);

  if (! allocate_arrays(na, nb))
    return 0;

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:atom data missing\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 BEGIN ATOM"))
  {
    cerr << "MDL_Molecule::_read_v3000:not BEGIN ATOM '" << buffer << "'\n";
    return 0;
  }

  int isave = atomic_symbols_can_have_arbitrary_length();
  set_atomic_symbols_can_have_arbitrary_length(1);
  int rc = _read_v3000(input, na, nb);
  set_atomic_symbols_can_have_arbitrary_length(isave);

  return rc;
}

int
MDL_Molecule::_read_v3000 (iwstring_data_source & input,
                           int na, 
                           int nb)
{
  IWString buffer;

  for (int i = 0; i < na; i++)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "MDL_Molecule::_read_v3000:premature eof during atoms, expected " << na << endl;
      return 0;
    }

    if (! _parse_v30_atom_record (buffer, 0))
    {
      cerr << "MDL_Molecule::_read_v3000:invalid atom record '" << buffer << "'\n";
      return 0;
    }

    const_IWSubstring token = buffer.word(2);

    int ndx;
    if (! token.numeric_value(ndx) || ndx < 1 || ndx > na)  // how could that happen?
      return 0;

    ndx--;

    token = buffer.word(3);
    if (! _convert_symbol_to_element(ndx, token))
    {
      cerr << "MDL_Molecule::_read_v3000:invalid atomic symbol '" << token << "'\n";
      cerr << buffer << endl;
      return 0;
    }

    if (! _look_for_atom_query_directives(ndx, buffer))
    {
      cerr << "MDL_Molecule::_read_v3000:invalid query specifiers '" << buffer << "'\n";
      return 0;
    }
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:premature eof after atoms\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 END ATOM"))
  {
    cerr << "MDL_Molecule::_read_v3000:should be END ATOM '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:premature eof at bonds\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 BEGIN BOND"))
  {
    cerr << "MDL_Molecule::_read_v3000:should be BEGIN BOND '" << buffer << "'\n";
    return 0;
  }

  int * aromatic_atom = new_int(na); iw_auto_array<int> free_aromatic_atom(aromatic_atom);
  int * aromatic_bond = new_int(na); iw_auto_array<int> free_aromatic_bond(aromatic_bond);

  for (int i = 0; i < nb; i++)
  {
    if (! read_next_v30_record(input, buffer))
    {
      cerr << "MDL_Molecule::_read_v3000:premature eof during bonds, expected " << nb << endl;
      return 0;
    }

    if (! _parse_v30_bond_record (buffer, aromatic_atom, aromatic_bond[i], 1))  // last arg means reading query file
    {
      cerr << "MDL_Molecule::_read_v3000:cannot parse bond record\n";
      cerr << buffer << endl;
      return 0;
    }

    _mdl_bond[i]->set_btype(bondi(i)->btype());

    if (! _look_for_bond_query_directives (i, buffer))
    {
      cerr << "MDL_Molecule::_read_v3000:invalid bond directives\n";
      cerr << buffer << endl;
      return 0;
    }
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:premature eof after bonds\n";
    return 0;
  }

  if (! buffer.starts_with("M  V30 END BOND"))
  {
    cerr << "MDL_Molecule::_read_v3000:should be END BOND '" << buffer << "'\n";
    return 0;
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:premature eof near end\n";
    return 0;
  }

  while (! buffer.starts_with ("M  V30 END CTAB"))
  {
    if (! input.next_record(buffer))
    {
      cerr << "MDL_Molecule::_read_v3000:no END CTAB found\n";
      return 0;
    }
  }

  if (! input.next_record(buffer))
  {
    cerr << "MDL_Molecule::_read_v3000:end record missing\n";
    return 0;
  }

  if (! buffer.starts_with("M  END"))
  {
    cerr << "MDL_Molecule::_read_v3000:not M  END '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
MDL_Molecule::_parse_v3000_atom_record (const const_IWSubstring & buffer)
{
  return 1;
}

int
MDL_Molecule::_convert_symbol_to_element (int ndx,
                                          const IWString & s) 
{
  if (s.length() <= 2)
  {
    int notused;
    const Element * e = get_element_from_symbol(s, notused);
    if (NULL == e)
    {
      cerr << "MDL_Molecule::_convert_symbol_to_element:cannot get element for '" << s << "'\n";
      return 0;
    }

    Molecule::set_element(ndx, e);

    return 1;
  }

  return _mdl_atom[ndx]->initialise_atom_list_from_symbol(s);
}

//M  V30 3 N 2.075 -0.975 0 0 SUBST=1

int
MDL_Molecule::_look_for_atom_query_directives (int ndx,
                                               const IWString & buffer)
{
  int i = 0;
  const_IWSubstring token;
  int col = 0;

  MDL_Atom_Data * mad = _mdl_atom[ndx];

  while (buffer.nextword(token, i))
  {
    col++;
    if (col < 7)
      continue;

    const_IWSubstring directive;
    int v;

    if (token.index('=') <= 0)
      continue;

    if (! token.split_into_directive_and_value(directive, '=', v))
      continue;

    if ("SUBST" == directive)
    {
      mad->set_substitution(v);
      continue;
    }
    else if ("UNSAT" == directive)
    {
      mad->set_unsaturation(v);
      continue;
    }
    else if ("RBCNT" == directive)
    {
      mad->set_ring_bond(v);
      continue;
    }
    else
    {
      cerr << "MDL_Molecule::_look_for_query_directives:unrecognised directive '" << token << "'\n";
      return 0;
    }
  }

  return 1;
}

int
MDL_Molecule::_look_for_bond_query_directives (int ndx,
                                               const IWString & buffer)
{
  int i = 0;
  const_IWSubstring token;

  MDL_Bond_Data * b = _mdl_bond[ndx];

  while (buffer.nextword(token, i))
  {
    if (token.index('=') <= 0)
      continue;

    const_IWSubstring directive;
    int v;
    if (! token.split_into_directive_and_value(directive, '=', v))
      continue;

    if ("TOPO" == directive)
    {
      b->set_bond_topology(v);
      continue;

    }
  }

  return 1;
}

/*
  Atom lists cause problems because there is no atomic number
*/

int
MDL_Molecule::compute_aromaticity_handle_atom_lists ()
{
  int matoms = natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const MDL_Atom_Data * mdlat = mdl_atom_data(i);

    const ISIS_Atom_List & alist = mdlat->atom_list();

    if (! alist.active())
      continue;

    if (atomic_number(i) > 0)
      continue;

    set_element(i, alist.elementi(0));
    rc++;
  }

  if (rc)
    compute_aromaticity();
  else
    compute_aromaticity_if_needed();
  return 1;
}

void
reset_mdl_molecule_file_scope_variables()
{
  convert_a_and_q_atoms_to_atom_lists=0;
  convert_not_atom_lists_to_organic_lists=0;
  element_a = NULL;
  element_q = NULL;
  element_l = NULL;
}