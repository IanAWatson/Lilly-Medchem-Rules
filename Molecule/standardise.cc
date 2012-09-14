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
/*
  Chemical standardisation routines.
*/

#include <stdlib.h>

#include "cmdline.h"
#include "iw_auto_array.h"
#include "misc.h"

#include "molecule.h"
#include "misc2.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "path.h"
#include "toggle_kekule_form.h"

class Molecule_Data_for_Standardisation
{
  private:
    int _natoms;
    int _nrings;
    int * _ncon;
    int * _atom_is_aromatic;
    int * _ring_is_aromatic;
    atomic_number_t _atomic_number;

  public:
    Molecule_Data_for_Standardisation();
    ~Molecule_Data_for_Standardisation();

    int initialise (Molecule & m);
};

/*
  In a multi-threaded environment, we cannot update accumulators
*/

static int update_accumulators = 1;

void
set_update_chemical_standardisation_accumulators (int s)
{
  update_accumulators = s;
}

Chemical_Transformation::Chemical_Transformation()
{
  _active = 0;
  _groups_changed = 0;
  _molecules_changed = 0;

  return;
}

Chemical_Transformation::~Chemical_Transformation()
{
  _active = 0;
  _groups_changed = -1;

  return;
}

/*
  A molecule has been processed, and N groups were changed in it.
*/

void
Chemical_Transformation::extra(int n)
{
  if (! update_accumulators)
    ;
  else if (n)
  {
    _groups_changed += n;
    _molecules_changed++;
  }

  return;
}

int
Chemical_Transformation::report(ostream & os) const
{
  os << "changed " << _groups_changed << " groups";
  if (_groups_changed)
    os << " in " << _molecules_changed << " molecules";

  os << endl;

  return os.good();
}

void
Chemical_Standardisation::_default_values()
{
  _ok = 1;
  _active = 0;

  _append_string_depending_on_what_changed = 0;

  _remove_hydrogens_attached_to_chiral_centres = 1;

  return;
}

Chemical_Standardisation::Chemical_Standardisation()
{
  _default_values();

  return;
}

Chemical_Standardisation::~Chemical_Standardisation()
{
  _ok = 0;

  return;
}

int
Chemical_Standardisation::ok() const
{
  return _ok;
}

int
Chemical_Standardisation::debug_print (ostream & os) const
{
  os << "Chemical Standardisation object\n";

  return os.good();
}

int
display_standard_chemical_standardisation_options (ostream & os, char zoption)
{
  os << "  -" << zoption << " <qualifier> chemical standardisations, enter \"-" << zoption << " help\" for usage\n";

  return os.good();
}

void
Chemical_Standardisation::activate_all()
{
  _transform_amines.activate();
  _transform_nplus_ominus.activate();    // no need to do nitro, as they are a subset
  _transform_n_charge_sep.activate();
  _remove_hydrogens.activate();
  _protonate_carboxyllic_acids.activate();
  _protonate_no.activate();
  _protonate_sulfonic_acids.activate();
  _protonate_sulfinic_acids.activate();
  _transform_splus_cminus.activate();
  _transform_ominus.activate();
  _transform_nminus.activate();
  _transform_covalent_metals.activate();
  _transform_guanidine.activate();
  _transform_guanidine_ring.activate();
  _transform_tetrazole.activate();
  _transform_azid.activate();
  _transform_misdrawn_urea.activate();
  _transform_imidazole.activate();
  _transform_pyrazole.activate();
  _transform_lactim_lactam.activate();
  _transform_lactim_lactam_ring.activate();
  _transform_triazole.activate();

  _active = 1;

  return;
}

#define CS_NITRO "nitro"
#define CS_NpOm  "n+o-"
#define CS_NpNm  "n+n-"
#define CS_SpCm  "s+c-"
#define CS_ALLpm "all+-"
#define CS_XH    "xh"
#define CS_NpH3  "n+h3"
#define CS_AMINE "amine"
#define CS_Om    "o-"
#define CS_Nm    "n-"
#define CS_ALL   "all"
#define CS_NRMCH "nrmch"
#define CS_COVM  "covm"
#define CS_ISOLC "isolc"
#define CS_GUAND "guan"
#define CS_GUANDR "Rguan"
#define CS_SPOM  "s+o-"
#define CS_ACID  "acid"
#define CS_EHLST "ehlast"
#define CS_FMRK  "fmrk"
#define CS_AZID  "azid"
#define CS_MSDUR "msdur"
#define CS_FCRN  "fcor"
#define CS_RNPNM "Rn+n-"
#define CS_FWIH  "fwih"
#define CS_IMIDAZOLE  "imidazole"
#define CS_PYRAZOLE  "pyrazole"
#define CS_TRIAZOLE  "triazole"
#define CS_TETRAZOLE  "tetrazole"
#define CS_LACTIM_LACTAM "ltlt"
#define CS_LACTIM_LACTAM_RING "ltltr"
#define CS_REVERSE_NITRO "rvnitro"
#define CS_REVERSE_NV5 "rvnv5"

int
display_all_chemical_standardisation_options (ostream & os, char zoption)
{
  os << "  -" << zoption << ' ' << CS_NITRO << "       transform nitro groups to N(=O)=O\n";
  os << "  -" << zoption << ' ' << CS_NpOm << "        transform charge separated [N+]-[O-] (includes nitro)\n";
  os << "  -" << zoption << ' ' << CS_NpNm << "        transform charge separated [N+]-[N-] to N=N\n";
  os << "  -" << zoption << ' ' << CS_SpCm << "        transform [S+]-[C-] to S=C\n";
  os << "  -" << zoption << ' ' << CS_ALLpm << "       transform all [X+]-[Y-] to X=Y\n";
  os << "  -" << zoption << ' ' << CS_XH << "          remove hydrogens\n";
  os << "  -" << zoption << ' ' << CS_AMINE << "       change all amines\n";
  os << "  -" << zoption << ' ' << CS_Om << "          protonate all O- groups\n";
  os << "  -" << zoption << ' ' << CS_Nm << "          protonate all N- groups\n";
  os << "  -" << zoption << ' ' << CS_ALL << "         ALL the above standardistions\n";
  os << "  -" << zoption << ' ' << CS_NRMCH << "       remove all hydrogens except those to chiral centres\n";
  os << "  -" << zoption << ' ' << CS_COVM << "        break covalent bonds between Oxygen and Na,K\n";
  os << "  -" << zoption << ' ' << CS_ISOLC << "       assign formal charges to isolated Na, K, .. and Halogens\n";
  os << "  -" << zoption << ' ' << CS_GUAND << "        convert guanidines to -N-C(=N)-N form\n";
  os << "  -" << zoption << ' ' << CS_GUANDR << "        convert Ring type guanidines to -N-C(=N)-N form\n";
  os << "  -" << zoption << ' ' << CS_AZID << "        convert charge separated azids [N-]=[N+]=N to N#N=N\n";
  os << "  -" << zoption << ' ' << CS_MSDUR << "       convert misdrawn ureas, O-C(=N)-N to O=C(-N)-N\n";
  os << "  -" << zoption << ' ' << CS_FCRN << "        for converting back from corina mangled structures\n";
  os << "  -" << zoption << ' ' << CS_EHLST << "      move all explicit Hydrogen atoms to last in the connection table\n";
  os << "  -" << zoption << ' ' << CS_FMRK << "        reverse transformations applied to .mrk files\n";
  os << "  -" << zoption << ' ' << CS_RNPNM << "       transform 5 valent N=N to charge separated form\n";
  os << "  -" << zoption << ' ' << CS_FWIH << "        fix obviously wrong implicit hydrogen settings\n";
  os << "  -" << zoption << ' ' << CS_IMIDAZOLE << "   convert imidazoles to have nH near cD3\n";
  os << "  -" << zoption << ' ' << CS_PYRAZOLE << "    convert pyrazoles to have nH near electron withdrawing\n";
  os << "  -" << zoption << ' ' << CS_TRIAZOLE << "    convert triazoles to have nH near electron withdrawing\n";
  os << "  -" << zoption << ' ' << CS_TETRAZOLE << "    convert tetrazoles to have nH near attachment\n";
  os << "  -" << zoption << ' ' << CS_LACTIM_LACTAM << "        convert lactim to lactam form (non ring)\n";
  os << "  -" << zoption << ' ' << CS_LACTIM_LACTAM_RING << "       convert lactim to lactam form (ring)\n";
  os << "  -" << zoption << ' ' << CS_REVERSE_NITRO << "     convert O=N=O nitro groups to charge separated\n";
  os << "  -" << zoption << ' ' << CS_REVERSE_NV5 << "       convert all 5 valent N atoms to charge separated\n";
  os << "  -" << zoption << ' ' << "APP=<xxx>   append 'xxx' to the name of changed molecules\n";
  os << "  -" << zoption << ' ' << "APP=EACH    append the reason for each change\n";

  return os.good();
}

#include "cmdline.h"

int
Chemical_Standardisation::construct_from_command_line (Command_Line & cl,
                                            int verbose,
                                            char flag)
{
  assert (ok());

  _verbose = verbose;

  if (verbose)
    set_update_chemical_standardisation_accumulators(1);

  int i = 0;
  IWString tmp;
  int rc = 0;
  while (cl.value (flag, tmp, i++))
  {
    if (CS_NITRO == tmp)
    {
      _transform_nitro.activate();
      if (verbose)
        cerr << "CS: nitro groups will be transformed to N(=O)=O\n";

      _active = 1;
      rc++;
    }
    else if (CS_NpOm == tmp)
    {
      _transform_nplus_ominus.activate();
      if (verbose)
        cerr << "CS: charge separated [N+]-[O-] will be transformed to N=O\n";

      _active = 1;
      rc++;
    }
    else if (CS_NpNm == tmp)
    {
      _transform_n_charge_sep.activate();
      if (verbose)
        cerr << "CS: charge separated [N+]-[N-] will be transformed to N=N\n";

      _active = 1;
      rc++;
    }
    else if (CS_ALLpm == tmp)
    {
      _transform_plus_minus.activate();
      if (verbose)
        cerr << "CS: all charge separated [X+]-[Y-] will be transformed to X=Y\n";

      _active = 1;
      rc++;
    }
    else if (CS_XH == tmp)
    {
      _remove_hydrogens.activate();
      if (verbose)
        cerr << "Hydrogens will be removed\n";

      _active = 1;
      rc++;
    }
    else if (CS_SpCm == tmp)
    {
      _transform_splus_cminus.activate();
      if (verbose)
        cerr << "[S+]-[C-] will be transformed to S=C\n";

      _active = 1;
      rc++;
    }
    else if (CS_Om == tmp)
    {
      _transform_ominus.activate();
      if (verbose)
        cerr << "All free [O-] groups will be protonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_Nm == tmp)
    {
      _transform_nminus.activate();
      if (verbose)
        cerr << "All [N-] groups will be protonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_AMINE == tmp)
    {
      _transform_amines.activate();
      if (verbose)
        cerr << "All charged amines will be deprotonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_COVM == tmp)
    {
      _transform_covalent_metals.activate();
      if (verbose)
        cerr << "Will break bonds to covalently bonded Na and K\n";

      _active = 1;
      rc++;
    }
    else if (CS_ISOLC == tmp)
    {
      _transform_single_atom_ions.activate();

      if (verbose)
        cerr << "Isolated metals and halogens will be assigned charges\n";

      _active = 1;

      rc++;
    }
    else if (CS_GUAND == tmp)
    {
      _transform_guanidine.activate();

      if (verbose)
        cerr << "Guanidines will be transformed\n";

      _active = 1;

      rc++;
    }
    else if (CS_GUANDR == tmp)
    {
      _transform_guanidine_ring.activate();

      if (verbose)
        cerr << "Ring guanidines will be transformed\n";

      _active = 1;

      rc++;
    }
    else if (CS_ALL == tmp)
    {
      activate_all();
      if (verbose)
        cerr << "All chemical standardisation transformations enabled\n";

      rc++;
    }
    else if (CS_NRMCH == tmp)
    {
      _remove_hydrogens.activate();
      _remove_hydrogens_attached_to_chiral_centres = 0;
      _active = 1;
      rc++;
    }
    else if (CS_ACID == tmp)
    {
      _protonate_carboxyllic_acids.activate();
      _protonate_sulfinic_acids.activate();
      _protonate_sulfonic_acids.activate();
      _protonate_sulfur_acids.activate();
      _protonate_phosphorous_acids.activate();
      _active = 1;
      rc++;

      if (verbose)
        cerr << "All acids will be protonated\n";
    }
    else if (CS_EHLST == tmp)
    {
      _explicit_hydrogens_last.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FMRK == tmp)
    {
      _from_mrk_standardisations.activate();
      _active = 1;
      rc++;
    }
    else if (CS_RNPNM == tmp)
    {
      _transform_back_to_nplus_nminus.activate();
      _transform_to_charge_separated_azid.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FWIH == tmp)
    {
      _transform_obvious_implicit_hydrogen_errors.activate();
      _active = 1;
      rc++;
    }
    else if (CS_AZID == tmp)
    {
      _transform_azid.activate();
      _active = 1;
      rc++;
    }
    else if (CS_MSDUR == tmp)
    {
      _transform_misdrawn_urea.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FCRN == tmp)
    {
      _transform_nitro.activate();
      _transform_azid.activate();
      _transform_nplus_ominus.activate();
      _active = 1;
      rc++;
    }
    else if (CS_IMIDAZOLE == tmp)
    {
      _transform_imidazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_TETRAZOLE == tmp)
    {
      _transform_tetrazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_PYRAZOLE == tmp)
    {
      _transform_pyrazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_TRIAZOLE == tmp)
    {
      _transform_triazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_LACTIM_LACTAM == tmp)
    {
      _transform_lactim_lactam.activate();
      _active = 1;
      rc++;
    }
    else if (CS_LACTIM_LACTAM_RING == tmp)
    {
      _transform_lactim_lactam_ring.activate();
      _active = 1;
      rc++;
    }
    else if (CS_REVERSE_NITRO == tmp)
    {
      _transform_nitro_reverse.activate();
      _active = 1;
      rc++;
    }
    else if (CS_REVERSE_NV5 == tmp)
    {
      _transform_nv5_to_charge_separated.activate();
      _active = 1;
      rc++;
    }
    else if (tmp.starts_with ("APP="))
    {
      tmp.remove_leading_chars (4);

      if ("EACH" == tmp)
      {
        _append_string_depending_on_what_changed = 1;

        if (verbose)
          cerr << "Will append standardisation details to molecule names\n";
      }
      else
      {
        _append_to_changed_molecules = tmp;

        if (_verbose)
          cerr << "Will append '" << _append_to_changed_molecules << "' to changed molecules\n";
      }
    }
    else if ("help" == tmp)
    {
      display_all_chemical_standardisation_options (cerr, flag);
      exit (0);    // note the very different behaviour in this case!!
    }
    else
    {
      cerr << "Chemical_Standardisation::construct...: unrecognised directive '" << tmp << "'\n";
      _ok = 0;
      _active = 0;
      return 0;
    }
  }

  return rc;
}

int
Chemical_Standardisation::_do_remove_hydrogens (Molecule & m)
{
  assert (ok());

  if (0 == m.natoms (1))     // quick check
    return 0;

  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Atom * a = m.atomi (i);

    if (1 != a->atomic_number())
      continue;

//  cerr << "H atom " << i << " has ncon " << a->ncon() << endl;
    if (a->ncon() > 1)     // don't touch those bridging hydrogens
      continue;

    if (a->formal_charge() && a->ncon())     // formally charged Hydrogen attached to something
    {
      if (_verbose)
        cerr << "Chemical_Standardisation::_do_remove_hydrogens: formally charged, covalently bonded Hydrogen not removed\n";
      continue;
    }

    if (a->is_isotope())
      continue;

    if (0 == _remove_hydrogens_attached_to_chiral_centres && a->ncon())
    {
      atom_number_t j = a->other (i, 0);
      if (NULL != m.chiral_centre_at_atom (j))
        continue;
    }

    m.remove_atom (i);
    rc++;
  }

  if (rc)
  {
    _remove_hydrogens.extra(rc);

    if (_verbose) 
      cerr << "Removed " << rc << " hydrogen atoms\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:remove H";
  }

  return rc;
}

/*
  Transform [N+][O-]=O to N(=O)=O
*/

int
Chemical_Standardisation::_do_transform_nitro (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;      // the number of nitro groups we change
  int natoms = m.natoms();

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])    // central atom must be N
      continue;
    if (3 != ncon[nitrogen])            // must have three connections
      continue;

    const Atom * an = atoms[nitrogen];

    if (1 != an->formal_charge())
      continue;
    if (4 != an->nbonds())          // all cases have 4 connections
      continue;

//  Now we must examine the atoms connected.

    int oxygens_attached = 0;
    int singly_bonded_oxygen = 0;
    int singly_bonded_oxygen_with_neg_charge = 0;
    int doubly_bonded_oxygen = 0;
    atom_number_t negative_oxygen = INVALID_ATOM_NUMBER;
    int negative_oxygen_j = -1;

    for (int j = 0; j < ncon[nitrogen]; j++)
    {
      const Bond * b = an->item (j);

      atom_number_t k = b->other (nitrogen);
      if (8 != z[k] || 1 != ncon[k])
        continue;

      oxygens_attached++;

      const Atom * ak = atoms[k];

      if (b->is_double_bond())
        doubly_bonded_oxygen++;

      else if (-1 == ak->formal_charge() && b->is_single_bond())
      {
        singly_bonded_oxygen_with_neg_charge++;
        negative_oxygen = k;
        negative_oxygen_j = j;
      }
      else if (0 == ak->formal_charge() && b->is_single_bond())
      {
        singly_bonded_oxygen++;
        negative_oxygen = k;
        negative_oxygen_j = j;
      }
    }

    if (2 != oxygens_attached || 1 != doubly_bonded_oxygen)
      continue; 

    if (singly_bonded_oxygen_with_neg_charge && singly_bonded_oxygen)   // cannot have both
      continue;

    if (singly_bonded_oxygen_with_neg_charge > 1 ||     // can have only one of either kind
        singly_bonded_oxygen > 1)
      continue;

    if (0 == singly_bonded_oxygen_with_neg_charge &&     // must have either kind
        0 == singly_bonded_oxygen)
      continue;

//  We have a nitro!

    m.set_formal_charge (nitrogen, 0);
    m.set_formal_charge (negative_oxygen, 0);
    m.set_bond_type_between_atoms  (nitrogen, negative_oxygen, DOUBLE_BOND);

//  adjust the global counters - not the oxygen because maybe it was neutral..

    current_molecule_data.change_nplus(-1);
    current_molecule_data.change_npos(-1);
    rc++;
  }

  if (rc)
  {
    _transform_nitro.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " nitro groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Nitro";
  }

  return rc;
}

/*
  Transform [N+][O-] to N(=O)
*/

int
Chemical_Standardisation::_do_transform_nplus_ominus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])    // central atom must be N
      continue;

    const Atom * an = atoms[nitrogen];

    if (1 != an->formal_charge())
      continue;

    for (int j = 0; j < ncon[nitrogen]; j++)
    {
      const Bond * b = an->item (j);
      if (! b->is_single_bond())
        continue;

      atom_number_t oxygen = b->other (nitrogen);

      if (8 != z[oxygen] || 1 != ncon[oxygen])
        continue;

      if (-1 != atoms[oxygen]->formal_charge())
        continue;

//    We have an [N+]-[O-] pair

      m.set_formal_charge (nitrogen, 0);
      m.set_formal_charge (oxygen, 0);
      m.set_bond_type_between_atoms (nitrogen, oxygen, DOUBLE_BOND);
      rc++;

      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_ominus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);


      break;     // this N is no longer positively charged! - Comgenex example of [O-][N+]([O-])=O for a nitro group
    }
  }

  if (rc)
  {
    _transform_nplus_ominus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N+]-[O-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N+]-[O-]";
  }

  return rc;
}

/*
  After all other transformations are done, there may be other [O-]
  atoms to transform

  Jul 97, also do S-
  Sept 98. Make sure these are singly connected - it was changing
  things like C[N+]1=CC(=NC)[O-]=N1 
  which destroyed the aromaticity
*/

int
Chemical_Standardisation::_do_transform_ominus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++)
  {
    if (8 == z[i])
      ;
    else if (16 == z[i])
      ;
    else
      continue;

    if (-1 != atoms[i]->formal_charge())
      continue;

    if (2 == ncon[i])
      continue;

    m.set_formal_charge (i, 0);
    current_molecule_data.change_nneg(-1);
    if (8 == z[i])
      current_molecule_data.change_ominus(-1);
    else
      current_molecule_data.change_sminus(-1);

    rc++;
  }

  if (rc)
  {
    _transform_ominus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [O-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[O-]";
  }

  return rc;
}

static int
bonded_to_positively_charged_nitrogen (const Molecule & m,
                 atom_number_t nitrogen)
{
  const Atom * a = m.atomi (nitrogen);

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (nitrogen, i);
    
    const Atom * aj = m.atomi (j);

    if (1 != aj->formal_charge())
      continue;

    if (7 != aj->atomic_number())
      continue;

    return 1;
  }

  return 0;
}

/*
  Change N- that aren't next to anything positive
*/

int
Chemical_Standardisation::_do_transform_nminus (Molecule & m,
                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])
      continue;

    const Atom * n = atoms[nitrogen];

    if (-1 != n->formal_charge())
      continue;

    if (bonded_to_positively_charged_nitrogen (m, nitrogen))
      continue;

    m.set_formal_charge (nitrogen, 0);
    rc++;
  }

  if (rc)
  {
    _transform_nminus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N-]";
  }

  return rc;
}

/*
  Transform [N+]-[N-] to N=N
  Make sure that the N+ has 4 bonds, and the N- two

  Ran into a big problem with this in the case where the N+ and N-
  are in an aromatic ring, so we don't do the transformation in
  those circumstances.

  Oct 2001. Don't want to do things like

  
*/

static int
ok_to_change_charge_separated_pair (const Molecule & m,
                        atom_number_t a1,
                        atom_number_t a2)
{
  const Atom * aa1 = m.atomi (a1);
  if (7 != aa1->atomic_number())     // not a Nitrogen, OK
    return 1;

  const Atom * aa2 = m.atomi (a2);
  if (7 != aa2->atomic_number())     // not a Nitrogen, OK
    return 1;

// we have a charge separated pair of Nitrogens

  int a1con = m.ncon (a1);
  int a2con = m.ncon (a2);

  if (2 == a1con && 1 == a2con)     // An azide in the form =[N+]=[N-]
    return 1;

  return 0;     // all other cases, don't change
}

int
Chemical_Standardisation::_do_transform_n_charge_sep (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change

// Ran into problems with LLY ...... where the transformation
// destroyed aromaticity. I then fixed the aromaticity rules, but
// to be careful, let's recompute the aromaticity if likely

  int need_to_recompute_aromaticity = 0;

  int natoms = m.natoms();
  for (int n1 = 0; n1 < natoms; n1++)
  {
    if (7 != z[n1])    // central atom must be N
      continue;

    Atom * a1 = const_cast<Atom *>(atoms[n1]);

    if (1 != a1->formal_charge())
      continue;
    if (4 != (a1->nbonds() + a1->implicit_hydrogens()))
      continue;

//  beware cases like C12=CC=CC=C1[N-][N+H2][N-]2 p10, don't change them

    int negative_nitrogens_attached = 0;   // multiple N- attached to one N+
    int negative_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[n1]; j++)
    {
      const Bond * b = a1->item (j);
      if (! b->is_single_bond())
        continue;

      atom_number_t n2 = b->other (n1);

      if (7 != z[n2])
        continue;

      Atom * an2 = const_cast<Atom *>(atoms[n2]);

      if (-1 != an2->formal_charge())
        continue;

      if (2 != (an2->nbonds() + an2->implicit_hydrogens()))
        continue;

//    We have an [N+]-[N-] pair. Are they in an aromatic ring
     
//    if (m.in_same_aromatic_ring (n1, n2))    August 2002, decided to drive out as much charge separation as possible
//      continue;

      negative_nitrogens_attached++;
      if (1 == negative_nitrogens_attached)
        negative_nitrogen = n2;
    }

    if (1 == negative_nitrogens_attached)
    {
      m.set_formal_charge(n1, 0);
      m.set_formal_charge(negative_nitrogen, 0);
      m.set_bond_type_between_atoms(n1, negative_nitrogen, DOUBLE_BOND);
      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);
      rc++;
    }
  }

  if (need_to_recompute_aromaticity)
    m.compute_aromaticity();

  if (rc)
  {
    _transform_n_charge_sep.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N+]-[N-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N+]-[N-]";
  }

  return rc;
}

static int
two_negatively_charged_connections (atom_number_t zatom,
                                    const Atom & a,
                                    const Atom * const * atoms,
                                    atom_number_t & n1,
                                    atom_number_t & n2)
{
  n1 = INVALID_ATOM_NUMBER;
  n2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < a.ncon(); i++)
  {
    const Bond * b = a.item (i);
    if (! b->is_single_bond())
      continue;

    atom_number_t n = b->other (zatom);

    if (-1 != atoms[n]->formal_charge())
      continue;

    if (INVALID_ATOM_NUMBER == n1)
      n1 = n;
    else
      n2 = n;
  }

  if (INVALID_ATOM_NUMBER == n1 || INVALID_ATOM_NUMBER == n2)
    return 0;

  return 1;
}

void
Chemical_Standardisation::_do_transform_plus_minus_pair (Molecule & m,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 IWStandard_Current_Molecule & current_molecule_data)
{
  m.set_formal_charge (a1, 0);
  m.set_formal_charge (a2, 0);
  m.set_bond_type_between_atoms (a1, a2, DOUBLE_BOND);

  current_molecule_data.change_npos(-1);
  current_molecule_data.change_nneg(-1);

  return;
}

/*
  Transform all [X+]-[Y-] to X=Y
*/

int
Chemical_Standardisation::_do_transform_plus_minus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++)
  {
    const Atom * ai = atoms[i];

    if (1 != ai->formal_charge())
      continue;

    atom_number_t n1, n2;
    (void) two_negatively_charged_connections (i, *ai, atoms, n1, n2);

    if (INVALID_ATOM_NUMBER == n1)    // no negatively charged neighbours
      continue;

    if (INVALID_ATOM_NUMBER == n2)    // only N1 is a negatively charged neighbour
    {
      if (ok_to_change_charge_separated_pair (m, i, n1))
        continue;

      _do_transform_plus_minus_pair (m, i, n1, current_molecule_data);
      rc++;
      continue;
    }

//  We have two negatively charged neighbours. Which one should we use?

    int ok1 = ok_to_change_charge_separated_pair (m, i, n1);
    int ok2 = ok_to_change_charge_separated_pair (m, i, n2);
//  cerr << ok1 << " and " << ok2 << endl;
    if (ok1 && ! ok2)
      _do_transform_plus_minus_pair (m, i, n1, current_molecule_data);
    else if (ok2 && ! ok1)
      _do_transform_plus_minus_pair (m, i, n2, current_molecule_data);
    else if (atoms[n1]->nbonds() < atoms[n2]->nbonds())
      _do_transform_plus_minus_pair (m, i, n1, current_molecule_data);
    else if (atoms[n1]->nbonds() > atoms[n2]->nbonds())
      _do_transform_plus_minus_pair (m, i, n2, current_molecule_data);
    else
      _do_transform_plus_minus_pair (m, i, n2, current_molecule_data);

    rc++;
  }

  if (rc)
  {
    _transform_plus_minus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [X+]-[Y-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[X+]-[Y-]";
  }

  return rc;
}

/*
  Change C(=O)[O-] to C(=O)O
*/

int
Chemical_Standardisation::_do_protonate_carboxyllic_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a carbon.

    const Bond * b = ai->item (0);

    if (! b->is_single_bond())
      continue;

    atom_number_t c = b->other (i);

    if (6 != z[c])
      continue;

    if (3 != ncon[c])
      continue;

    const Atom * ac = atoms[c];

    int found_doubly_bonded_oxygen = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = ac->item (j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t o = b2->other (c);

      if (8 == z[o] && 1 == ncon[o])
      {
        found_doubly_bonded_oxygen = 1;
        break;
      }
    }

    if (! found_doubly_bonded_oxygen)
      continue;

    m.set_formal_charge (i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_carboxyllic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " carboxyllic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[OH]-C=O";
  }

  return rc;
}

/*
  Change [O-,S-]-P=[O,S] to [O,S]-P=[O,S]
*/

int
Chemical_Standardisation::_do_protonate_phosphorous_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * ai = atoms[i];

//  Find a negatively charged Oxygen or Sulphur, and a doubly bonded O or S - let's hope we don't see too many P~S bonds, ugly!

    atom_number_t negatively_charged_OS = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_OS = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item (j);

      atom_number_t k = b->other (i);

      const Atom * ak = atoms[k];

      if (8 == ak->atomic_number())
        ;
      else if (16 == ak->atomic_number())
        ;
      else
        continue;

      if (b->is_single_bond())
      {
        if (-1 == ak->formal_charge())
          negatively_charged_OS = k;
      }
      if (b->is_double_bond())
        doubly_bonded_OS = k;
      else if (b->is_single_bond() && -1 == ak->formal_charge())
        negatively_charged_OS = k;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_OS || INVALID_ATOM_NUMBER == negatively_charged_OS)
      continue;

    m.set_formal_charge (negatively_charged_OS, 0);
    rc++;

    current_molecule_data.change_nneg(-1);

    if (8 == atoms[negatively_charged_OS]->atomic_number())
      current_molecule_data.change_ominus(-1);
    else
      current_molecule_data.change_sminus(-1);

    current_molecule_data.change_phosphorus(-1);
    if (0 == current_molecule_data.phosphorus())
      break;
  }

  if (rc)
  {
    _protonate_phosphorous_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " phosphorus acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[OH,SH]-P=[O,S]";
  }

  return rc;
}

/*
  Change [S-]-[C,P]=[O,S] to S-[C,P]=[O,S]
*/

int
Chemical_Standardisation::_do_protonate_sulfur_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (16 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

    const Bond * b = ai->item (0);

    if (! b->is_single_bond())
      continue;

    atom_number_t cps = b->other (i);

    if (6 == z[cps])
      ;
    else if (15 == z[cps])
      ;
    else if (16 == z[cps])
      ;
    else
      continue;

    if (ncon[cps] < 3)
      continue;

    const Atom * acps = atoms[cps];

    int found_doubly_bonded_OS = 0;
    for (int j = 0; j < ncon[cps]; j++)
    {
      const Bond * b2 = acps->item (j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t os = b2->other (cps);

      if (1 != ncon[os])
        continue;

      if (8 == z[os] || 16 == z[os])
      {
        found_doubly_bonded_OS = 1;
        break;
      }
    }

    if (! found_doubly_bonded_OS)
      continue;

    m.set_formal_charge (i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_sminus(-1);
  }
  if (rc)
  {
    _protonate_sulfur_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " Sulphur acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[SH]-[C,P,S]=[O,S]";
  }

  return rc;
}

int
Chemical_Standardisation::_do_protonate_sulfonic_acids (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a sulphur.

    const Bond * b1 = ai->item (0);

    if (! b1->is_single_bond())
      continue;

    atom_number_t s = b1->other (i);

    if (16 != z[s])
      continue;

    if (4 != ncon[s])
      continue;

    const Atom * as = atoms[s];

    int doubly_bonded_oxygen_count = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = as->item (j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t o = b2->other (s);

      if (8 == z[o] && 1 == ncon[o])
        doubly_bonded_oxygen_count++;
    }

    if (2 != doubly_bonded_oxygen_count)
      continue;

    m.set_formal_charge (i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_sulfonic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " sulfonic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Sulfonic";
  }

  return rc;
}

/*
  A sulfinic acid is S(=O)O
*/

int
Chemical_Standardisation::_do_protonate_sulfinic_acids (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a sulphur.

    const Bond * b = ai->item (0);
    if (! b->is_single_bond())
      continue;

    atom_number_t s = b->other (i);

    if (16 != z[s])
      continue;

    if (4 != ncon[s])
      continue;

    const Atom * as = atoms[s];

    int doubly_bonded_oxygen_count = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = as->item (j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t o = b2->other (s);

      if (8 == z[o] && 1 == ncon[o])
        doubly_bonded_oxygen_count++;
    }

    if (1 != doubly_bonded_oxygen_count)
      continue;

    m.set_formal_charge (i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_sulfinic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " sulfinic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Sulfinic";
  }

  return rc;
}

int
Chemical_Standardisation::_do_protonate_no (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a Nitrogen.

    const Bond * b = ai->item (0);
    if (! b->is_single_bond())
      continue;

    atom_number_t n = b->other (i);;
    if (7 != z[n])
      continue;

    m.set_formal_charge (i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_no.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " NO- groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:N-[O-]";
  }

  return rc;
}

/*
  Transform [C-]-[S+] to C=S
  if the S+ has 5 bonds
  October 2000. Relax the restrictions on the environment of the S. Now allow either
    4 connections and 5 bonds
  or
    3 connections and 3 bonds

  Feb 2001. Extend so that we also transform [O-]-[S+]
*/

int
Chemical_Standardisation::_do_transform_splus_cminus (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * ai = atoms[i];

    if (1 != ai->formal_charge())
      continue;

    if (4 == ncon[i] && 5 == ai->nbonds())     // yep, do this one
      ;
    else if (3 == ncon[i] && 3 == ai->nbonds())     // yep, do this one
      ;
    else
      continue;

//  At this stage, we have an appropriate S+. Is it bonded to a C- or O-

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item (j);
      if (! b->is_single_bond())
        continue;

      atom_number_t c = b->other (i);
      if (6 == z[c])
        ;
      else if (8 == z[c])
        ;
      else
        continue;

      if (-1 != atoms[c]->formal_charge())
        continue;

      m.set_formal_charge (i, 0);
      m.set_formal_charge (c, 0);
      m.set_bond_type_between_atoms (i, c, DOUBLE_BOND);
      rc++;

      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);

      break;
    }
  }

  if (rc)
  {
    _transform_splus_cminus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [S+]-[C-]\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[S+]-[C-]";
  }

  return rc;
}

/*
  This next one is greatly complicated by the possibility of explicit hydrogens
*/

int
Chemical_Standardisation::_do_transform_amines (Molecule & m,
                                        Set_of_Atoms & atoms_to_be_removed,
                                        IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    Atom * a = const_cast<Atom *>(atoms[i]);

    if (1 != a->formal_charge())
      continue;

//  It must have a hydrogen.

    int hc = m.hcount (i);     // check both implicit and explicit

    if (0 == hc)     // must have at least 1 Hydrogen attachment
      continue;

    int implicit_hydrogens;
    if (hc)
      implicit_hydrogens = a->implicit_hydrogens();
    else
      implicit_hydrogens = 0;

//  If all hydrogens are explicit, we must find one to remove. Check also for
//  adjacent N- which do not get processed

    atom_number_t explicit_hydrogen_to_remove = INVALID_ATOM_NUMBER;
    int found_negative_nitrogen = 0;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other (i, j);
      if (1 == z[k])
        explicit_hydrogen_to_remove = k;
      else if (7 == z[k] && -1 == atoms[k]->formal_charge())
      {
        found_negative_nitrogen = 1;
        break;
      }
    }

    if (found_negative_nitrogen)
      continue;

    if (implicit_hydrogens)
      ;
    else if (INVALID_ATOM_NUMBER == explicit_hydrogen_to_remove)   // should not happen
      continue;
    else
      atoms_to_be_removed.add_if_not_already_present(explicit_hydrogen_to_remove);

    m.set_formal_charge (i, 0);

    rc++;

    current_molecule_data.change_nplus(-1);
    current_molecule_data.change_npos(-1);
  }

  if (rc)
  {
    _transform_amines.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " protonated amines\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:amines";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_covalent_metals (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  int matoms = m.natoms();

  int rc = 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int i = 0; i < matoms; i++)
  {
    if (1 != ncon[i])       // we only process single atoms
      continue;

    if (11 == z[i] || 19 == z[i])     // Na and K only
    {
      atom_number_t o = atoms[i]->other (i, 0);

      if (8 == z[o])
        ;
      else if (16 == z[o])
        ;
      else
        continue;

      if (2 != ncon[o])
        continue;

      m.remove_bond_between_atoms (i, o);
      m.set_formal_charge (i, 1);
      m.set_formal_charge (o, -1);
      ncon[i] = 0;
      ncon[o] = 1;
      current_molecule_data.change_nneg(1);
      if (8 == z[o])
        current_molecule_data.change_ominus(1);
      else
        current_molecule_data.change_sminus(1);

      rc++;

      current_molecule_data.change_singly_connected_metal(-1);

      if (0 == current_molecule_data.singly_connected_metal())
        break;
    }
  }

  if (rc)
  {
    _transform_covalent_metals.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " covalently bonded metals\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:cv_metal";
  }

  return rc;
}

/*
  We convert -N=C(-NH2)-NH2  to -N-C(=NH)-NH2
*/

int
Chemical_Standardisation::_do_transform_guanidine (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();

  int matoms = m.natoms();

//cerr << "Checking " << matoms << " atoms for guanidines\n";

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (2 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (3 != ai->nbonds())
      continue;

    if (ring_membership[i])
      continue;

    atom_number_t doubly_bonded_carbon = INVALID_ATOM_NUMBER;

//  cerr << "Atom " << i << " is a candidate Nitrogen atom\n";

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = ai->item (j);

      if (! b->is_double_bond())
        continue;

      if (b->part_of_cis_trans_grouping())
        continue;

      atom_number_t k = b->other (i);

      if (3 == ncon[k] && 6 == z[k] && 0 == ring_membership[k])
      {
        doubly_bonded_carbon = k;
        break;
      }
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_carbon)
      continue;

//  We have -N=[CD3]. Is it a guanidine?

    int nh2_found = 0;
    atom_number_t second_nh2 = INVALID_ATOM_NUMBER;

    const Atom * c = atoms[doubly_bonded_carbon];
    if (4 != c->nbonds())
      continue;

    for (int l = 0; l < 3; l++)
    {
      atom_number_t n = c->other (doubly_bonded_carbon, l);

      if (i == n)
        continue;

      if (7 == z[n] && 1 == ncon[n])
      {
        second_nh2 = n;
        nh2_found++;
      }
    }

    if (2 != nh2_found)
      continue;

    m.set_bond_type_between_atoms (i, doubly_bonded_carbon, SINGLE_BOND);
    m.set_bond_type_between_atoms (doubly_bonded_carbon, second_nh2, DOUBLE_BOND);
    m.set_implicit_hydrogens_known (i, 0);
    m.set_implicit_hydrogens_known (second_nh2, 0);

//  Lose any formal charge already applied

    if (1 == ai->formal_charge())
      m.set_formal_charge (i, 0);

    rc++;
  }

  if (rc)
  {
    _transform_guanidine.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " guanidine groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:guanidine";
  }

  return rc;
}

static int
compute_guanidine_bond_acceptance_desirability (Molecule & m,
                                                atom_number_t n,
                                                const int * ring_membership,
                                                const int * ncon)
{
  if (0 == m.hcount(n))
    return 0;
  if (1 == ncon[n])
    return 10;
  if (0 == ring_membership[n])
    return 5 + m.attached_heteroatom_count(n);

   return m.attached_heteroatom_count(n);
}

/*
  We have something like


  N     N
   \   /
    \ /
     C
     ||
     N

  where the N=C bond is in a ring. We want to put the double bond outside the ring
*/

int
Chemical_Standardisation::_do_transform_ring_guanidine (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  if (0 == m.nrings())
    return 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  m.compute_aromaticity();    // must compute it. Molecule may have aromaticity definition computed with Pearlman rules

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (2 != ncon[i])
      continue;

    if (0 == ring_membership[i])   // we only consider these in a ring
      continue;

//  cerr << "Atom " << i << " aromaticity " << m.is_aromatic(i) << endl;

    if (atom_is_aromatic[i])
      continue;

    const Atom * ai = atoms[i];

    if (3 != ai->nbonds())
      continue;

    atom_number_t doubly_bonded_carbon = INVALID_ATOM_NUMBER;

//  cerr << "Atom " << i << " is a candidate Nitrogen atom\n";

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = ai->item (j);

      if (! b->is_double_bond())
        continue;

      if (b->part_of_cis_trans_grouping())
        continue;

      atom_number_t k = b->other (i);

      if (3 == ncon[k] && 6 == z[k] && ring_membership[k])
      {
        doubly_bonded_carbon = k;
        break;
      }
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_carbon)
      continue;

//  We have -N=[CD3]. Need to find the other two Nitrogens.
//  Give preference to singly connected exocyclic

    const Atom * c = atoms[doubly_bonded_carbon];
    if (4 != c->nbonds())
      continue;

    atom_number_t n1 = INVALID_ATOM_NUMBER;
    int n1_desirability = 0;
    atom_number_t n2 = INVALID_ATOM_NUMBER;
    int n2_desirability = 0;

    for (int l = 0; l < 3; l++)
    {
      atom_number_t n = c->other (doubly_bonded_carbon, l);

      if (7 != z[n])   // all attachments must be N atoms
        break;

      if (i == n)    // the first Nitrogen
        continue;

      if (INVALID_ATOM_NUMBER == n1)
      {
        n1 = n;
        n1_desirability = compute_guanidine_bond_acceptance_desirability (m, n, ring_membership, ncon);
      }
      else
      {
        n2 = n;
        n2_desirability = compute_guanidine_bond_acceptance_desirability (m, n, ring_membership, ncon);
      }
    }

    if (INVALID_ATOM_NUMBER == n1 || INVALID_ATOM_NUMBER == n2)
      continue;

    int n3_desirability = compute_guanidine_bond_acceptance_desirability(m, i, ring_membership, ncon);

    if (n3_desirability >= n1_desirability && n3_desirability >= n2_desirability)
      continue;

    atom_number_t n;
    if (n2_desirability > n1_desirability && n2_desirability > n3_desirability)
      n = n2;
    else if (n1_desirability > n2_desirability && n1_desirability > n3_desirability)
      n = n1;
    else
      continue;

//  cerr << "Ring guanidine " << i << " (" << doubly_bonded_carbon << ")(" << n1 << ")" << n2 << endl;

    m.set_bond_type_between_atoms (i, doubly_bonded_carbon, SINGLE_BOND);
    m.set_bond_type_between_atoms (doubly_bonded_carbon, n, DOUBLE_BOND);
    m.set_implicit_hydrogens_known (i, 0);
    m.set_implicit_hydrogens_known (n, 0);

//  Lose any formal charge already applied

    if (1 == ai->formal_charge())
      m.set_formal_charge (i, 0);

    rc++;
  }

  if (rc)
  {
    _transform_guanidine_ring.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " ring guanidine groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Rguanidine";
  }

  return rc;
}

/*
  For Research Records, we often want counterions to be properly charged

  We put + charges on
   Na, K
  Two + on
   Ca
  Negative charges on
   F, Cl, Br, I
*/

int
Chemical_Standardisation::_do_transform_single_atom_ions (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  int matoms = m.natoms();

  int rc = 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 != ncon[i])       // we only process isolated atoms
      continue;

    if (atoms[i]->formal_charge())    // already set, we don't change it
      continue;

    if (atoms[i]->implicit_hydrogens_known() && const_cast<Atom *>(atoms[i])->implicit_hydrogens() > 0)    // don't change it
      continue;

    if (atoms[i]->element()->is_halogen())
    {
      m.set_formal_charge (i, -1);

      rc++;
    }
    else if (! atoms[i]->element()->is_metal())     // we only process metals here
      ;
    else if (3 == z[i] || 11 == z[i] || 19 == z[i])     // Li, Na, K
    {
      m.set_formal_charge (i, 1);

      rc++;
    }
    else if (12 == z[i] || 20 == z[i] || 30 == z[i] || 56 == z[i])   // Mg, Ca, Zn, Ba
    {
      m.set_formal_charge (i, 2);

      rc++;
    }
  }

  if (rc)
  {
    _transform_single_atom_ions.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " isolated metals/halogens\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:isolated";
  }

  return rc;
}

int
Chemical_Standardisation::_do_explicit_hydrogens_last (Molecule & m)
{
  int matoms = m.natoms();

// Since other standardisations may have removed Hydrogen atoms, we count
// the number still present

  atom_number_t last_non_hydrogen_atom = INVALID_ATOM_NUMBER;

  Set_of_Atoms eh;
  for (int i = 0; i < matoms; i++)
  {
    if (1 == m.atomic_number (i))
      eh.add (i);
    else 
      last_non_hydrogen_atom = i;
  }

  if (0 == eh.number_elements())    // none present, that's easy
    return 0;

  if (INVALID_ATOM_NUMBER == last_non_hydrogen_atom)    // Huh, a molecule with just Hydrogen atoms
    return 0;

// Any explicit Hydrogen that is already past the last non-hydrogen atom doesn't need to be removed

  for (int i = eh.number_elements() - 1; i >= 0; i--)
  {
    atom_number_t h = eh[i];

    if (h > last_non_hydrogen_atom)   // doesn't need to be moved around
      eh.remove_item (i);
  }

  int neh = eh.number_elements();

  if (0 == neh)    // order must have already been OK
    return 1;

  for (int i = neh - 1; i >= 0; i--)
  {
    atom_number_t h = eh[i];

    m.move_atom_to_end_of_atom_list (h);
  }

  _explicit_hydrogens_last.extra(neh);

  if (neh)
  {
    if (_verbose)
      cerr << "Moved " << neh << " explicit Hydrogens\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:ehlast";
  }

  return neh;
}


int
Chemical_Standardisation::_process (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;

  if (_transform_nitro_reverse.active())
    rc += _do_transform_reverse_nitro (m, current_molecule_data);

  if (_transform_azid_reverse.active())
    rc += _do_transform_reverse_azid (m, current_molecule_data);

  if (_transform_nv5_to_charge_separated.active())
    rc += _do_nv5_to_charge_separated(m, current_molecule_data);

  if (_transform_misdrawn_urea.active())
    rc += _do_transform_misdrawn_urea (m, current_molecule_data);

  if (_transform_nitro.active() && current_molecule_data.nplus() && current_molecule_data.ominus())
    rc += _do_transform_nitro (m, current_molecule_data);

  if (_transform_covalent_metals.active() && current_molecule_data.singly_connected_metal())
    rc += _do_transform_covalent_metals (m, current_molecule_data);

  if (_transform_single_atom_ions.active())
    rc += _do_transform_single_atom_ions (m, current_molecule_data);

  if (_transform_nplus_ominus.active() && current_molecule_data.nplus() && current_molecule_data.ominus())
    rc += _do_transform_nplus_ominus (m, current_molecule_data);

  if (_transform_plus_minus.active() && current_molecule_data.nneg() && current_molecule_data.npos())
    rc += _do_transform_plus_minus (m, current_molecule_data);

  if (_transform_azid.active() && current_molecule_data.nneg() && current_molecule_data.npos() && current_molecule_data.nitrogens() >= 2)
    rc += _do_transform_azid (m, current_molecule_data);

  if (_transform_n_charge_sep.active() && current_molecule_data.nneg() && current_molecule_data.nplus())
    rc += _do_transform_n_charge_sep (m, current_molecule_data);

  if (_protonate_no.active() && current_molecule_data.ominus())
    rc += _do_protonate_no (m, current_molecule_data);

  if (_protonate_carboxyllic_acids.active() && (current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_protonate_carboxyllic_acids (m, current_molecule_data);

// The order of the sulphur acid types is important.
// Make sure sulfonic are done before sulfinic

  if (_protonate_sulfonic_acids.active() && current_molecule_data.ominus() && current_molecule_data.sulphur())
    rc += _do_protonate_sulfonic_acids (m, current_molecule_data);

  if (_protonate_sulfinic_acids.active() && current_molecule_data.ominus() && current_molecule_data.sulphur())
    rc += _do_protonate_sulfonic_acids (m, current_molecule_data);

  if (_protonate_sulfur_acids.active() && current_molecule_data.sminus())
    rc += _do_protonate_sulfur_acids (m, current_molecule_data);

  if (_protonate_phosphorous_acids.active() && current_molecule_data.phosphorus() && (current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_protonate_phosphorous_acids (m, current_molecule_data);

  if (_transform_splus_cminus.active() && current_molecule_data.nneg() && current_molecule_data.npos() && current_molecule_data.sulphur() && current_molecule_data.cminus())
    rc += _do_transform_splus_cminus (m, current_molecule_data);

  if (_transform_guanidine.active() && current_molecule_data.possible_guanidine())
    rc += _do_transform_guanidine (m, current_molecule_data);

  if (_transform_guanidine_ring.active() && current_molecule_data.possible_guanidine())
    rc += _do_transform_ring_guanidine (m, current_molecule_data);

  if (_transform_tetrazole.active())
    rc += _do_tetrazole (m, current_molecule_data);

  if (_transform_imidazole.active())
    rc += _do_imidazole (m, current_molecule_data);

  if (_transform_pyrazole.active())
    rc += _do_pyrazole (m, current_molecule_data);

  if (_transform_triazole.active())
    rc += _do_triazole (m, current_molecule_data);

  if (_from_mrk_standardisations.active())
    rc += _do_from_mrk_standardisations (m, current_molecule_data);

// Explicit Hydrogens may get removed from Amines

  Set_of_Atoms atoms_to_be_removed;

  if (_transform_amines.active() && current_molecule_data.nplus())
    rc += _do_transform_amines (m, atoms_to_be_removed, current_molecule_data);

// Make sure we do O- after everything else

  if (_transform_ominus.active() && (current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_transform_ominus (m, current_molecule_data);

  if (_transform_lactim_lactam.active() && current_molecule_data.possible_lactim())
    rc += _do_transform_lactim(m, 0, current_molecule_data);

  if (_transform_lactim_lactam_ring.active() && current_molecule_data.possible_lactim())
    rc += _do_transform_lactim(m, 1, current_molecule_data);

  if (_transform_back_to_nplus_nminus.active() && current_molecule_data.nitrogens() > 1)
    rc += _do_transform_back_to_nplus_nminus (m, current_molecule_data);

  if (_transform_to_charge_separated_azid.active() && current_molecule_data.nitrogens() > 2)
    rc += _do_transform_azid_to_charge_separated (m, current_molecule_data);

  if (_transform_obvious_implicit_hydrogen_errors.active())
    rc += _do_transform_implicit_hydrogen_known_errors (m, current_molecule_data);

// Nminus must also be done after most other things

  if (_transform_nminus.active() && current_molecule_data.nneg())
    rc += _do_transform_nminus (m, current_molecule_data);

  assert (current_molecule_data.npos() >= 0);
  assert (current_molecule_data.nneg() >= 0);
  assert (current_molecule_data.nplus() >= 0);
  assert (current_molecule_data.ominus() >= 0);
  assert (current_molecule_data.sminus() >= 0);

  if (atoms_to_be_removed.number_elements())
    m.remove_atoms (atoms_to_be_removed);

// Since this changes the order of the atoms we don't use any of the arrays

  if (current_molecule_data.explicit_hydrogen_count() && _explicit_hydrogens_last.active())
    rc += _do_explicit_hydrogens_last (m);

  return rc;
}

int
Chemical_Standardisation::process (Molecule & m) 
{
  if (0 == _active)
    return 1;

  int asave = global_aromaticity_type();

  if (Daylight != asave)
    set_global_aromaticity_type(Daylight);

  int rc = _process(m);

  set_global_aromaticity_type(asave);

  return rc;
}

int
Chemical_Standardisation::_process(Molecule & m)
{
  _molecules_processed++;

// Removing atoms can mess up anything which comes after, so make 
// sure we do that before anything else

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules.resize_keep_storage (0);

  int rc = 0;

  if (_remove_hydrogens.active())
    rc += _do_remove_hydrogens (m);

  if (0 == m.natoms())
    return 0;

  IWStandard_Current_Molecule current_molecule_data;

  if (! current_molecule_data.initialise(m))
  {
    cerr << "Chemical_Standardisation::process:could not initialise '" << m.name() << "'\n";
    return 0;
  }

  if (! _processing_needed (current_molecule_data))
    return 0;

  rc += _process (m, current_molecule_data);

  if (rc)
  {
    _molecules_changed++;
    m.recompute_implicit_hydrogens();

    if (_append_to_changed_molecules.length())
      m.append_to_name (_append_to_changed_molecules);
  }

  return rc;
}

int
Chemical_Standardisation::report(ostream & os) const
{
  os << "Report on chemical standardisation object\n";

  if (0 == _molecules_processed)
    return os.good();

  if (_transform_amines.active())
  {
    os << "  Amines ";
    _transform_amines.report(os);
  }
  if (_transform_nitro.active())
  {
    os << "  Nitro groups ";
    _transform_nitro.report(os);
  }

  if (_transform_nplus_ominus.active())
  {
    os << "  N+O- ";
    _transform_nplus_ominus.report(os);
  }

  if (_transform_nv5_to_charge_separated.active())
  {
    os << "  Nv5 ";
    _transform_nv5_to_charge_separated.report(os);
  }

  if (_transform_plus_minus.active())
  {
    os << "  Plus Minus ";
    _transform_plus_minus.report(os);
  }

  if (_transform_covalent_metals.active())
  {
    os << "  Covalent metals ";
    _transform_covalent_metals.report(os);
  }

  if (_transform_n_charge_sep.active())
  {
    os << "  Nitrogen charge separated ";
    _transform_n_charge_sep.report(os);
  }

  if (_remove_hydrogens.active())
  {
    os << "  Remove hydrogens ";
    _remove_hydrogens.report(os);
  }

  if (_protonate_carboxyllic_acids.active())
  {
    os << "  Protonate carboxyllic acids ";
    _protonate_carboxyllic_acids.report(os);
  }

  if (_protonate_no.active())
  {
    os << "  Protonate NO ";
    _protonate_no.report(os);
  }

  if (_protonate_sulfonic_acids.active())
  {
    os << "  Protonate sulfonic acids ";
    _protonate_sulfonic_acids.report(os);
  }

  if (_protonate_sulfinic_acids.active())
  {
    os << "  Protonate sulfinic acids ";
    _protonate_sulfinic_acids.report(os);
  }

  if (_transform_splus_cminus.active())
  {
    os << "  [S+]-[C-] ";
    _transform_splus_cminus.report(os);
  }

  if (_transform_ominus.active())
  {
    os << "  Protonate O- ";
    _transform_ominus.report(os);
  }

  if (_transform_nminus.active())
  {
    os << "  Protonate N- ";
    _transform_nminus.report(os);
  }

  if (_protonate_sulfinic_acids.active())
  {
    os << "  Protonate [S-]-[C,S,P]=[C,S] ";
    _protonate_sulfinic_acids.report(os);
  }

  if (_protonate_phosphorous_acids.active())
  {
    os << "  Protonate [O,S]-P=[O,S] ";
    _protonate_phosphorous_acids.report(os);
  }

  if (_explicit_hydrogens_last.active())
  {
    os << "  Explicit Hydrogens last ";
    _explicit_hydrogens_last.report(os);
  }

  if (_transform_tetrazole.active())
  {
    os << "  Tetrazole ";
    _transform_tetrazole.report(os);
  }

  if (_transform_triazole.active())
  {
    os << "  Triazole ";
    _transform_triazole.report(os);
  }

  if (_transform_imidazole.active())
  {
    os << "  Imidazole ";
    _transform_imidazole.report(os);
  }

  if (_transform_pyrazole.active())
  {
    os << "  Pyrazole ";
    _transform_pyrazole.report(os);
  }

  if (_transform_guanidine.active())
  {
    os << "  Guanidine ";
    _transform_guanidine.report(os);
  }

  if (_transform_guanidine_ring.active())
  {
    os << "  Guanidine Ring ";
    _transform_guanidine_ring.report(os);
  }

  if (_transform_lactim_lactam.active())
  {
    os << "  Lactim_lactam ";
    _transform_lactim_lactam.report(os);
  }

  if (_transform_lactim_lactam_ring.active())
  {
    os << "  Lactim_lactam ring ";
    _transform_lactim_lactam_ring.report(os);
  }

  if (_transform_azid.active())
  {
    os << "  Azid ";
    _transform_azid.report(os);
  }

  if (_transform_misdrawn_urea.active())
  {
    os << "  MSDUR ";
    _transform_misdrawn_urea.report(os);
  }

  if (_transform_back_to_nplus_nminus.active())
  {
    os << "=N=N- to =[N+]-[N-] ";
    _transform_back_to_nplus_nminus.report(os);
  }

  return os.good();
}

int
Chemical_Standardisation::deactivate (const const_IWSubstring & d)
{
  if (CS_NITRO  == d)
    _transform_nitro.deactivate();
  else if (CS_NpOm == d)
    _transform_nplus_ominus.deactivate();
  else if (CS_NpNm == d)
      _transform_nplus_ominus.deactivate();
  else if (CS_SpCm == d)
    _transform_splus_cminus.deactivate();
  else if (CS_ALLpm == d)
    _transform_plus_minus.deactivate();
  else if (CS_XH == d)
    _remove_hydrogens.deactivate();
  else if (CS_AMINE == d)
    _transform_amines.deactivate();
  else if (CS_Om == d)
    _transform_ominus.deactivate();
  else if (CS_Nm == d)
    _transform_nminus.deactivate();
  else
  {
    cerr << "Chemical_Standardisation::deactivate: unrecognised directive '" << d << "'\n";
    return 0;
  }

  return 1;
}

int
Chemical_Standardisation::activate_all_except_hydrogen_removal()
{
  activate_all();

  _remove_hydrogens.deactivate();

  return 1;
}

/*
  Change N#N=N to [N-]=[N+]=N
*/

int
Chemical_Standardisation::_do_transform_reverse_azid (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    const Bond * b = ai->item (0);

    if (! b->is_triple_bond())
      continue;

    atom_number_t n2 = b->other (i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * an2 = atoms[n2];

    if (5 != an2->nbonds())
      continue;

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = an2->item (j);

      if (! b->is_double_bond())
        continue;

      atom_number_t n3 = b->other (n2);

      if (7 != z[n3])
        continue;

      m.set_bond_type_between_atoms (i, n2, DOUBLE_BOND);
      m.set_formal_charge (i, -1);
      m.set_formal_charge (n2, 1);

      rc++;
    }
  }

  if (rc)
  {
    _transform_azid_reverse.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " azid groups to charge separated\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:azid+-";
  }

  return rc;
}

/*
  Change [N-]=[N+]=N to N#N=N
*/

int
Chemical_Standardisation::_do_transform_azid (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

    const Bond * b = ai->item (0);

    if (! b->is_double_bond())
      continue;

    atom_number_t n2 = b->other (i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * an2 = atoms[n2];

    if (1 != an2->formal_charge())
      continue;

    if (4 != an2->nbonds())
      continue;

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = an2->item (j);

      if (! b->is_double_bond())    // how could that happen?
        continue;

      atom_number_t n3 = b->other (n2);
      if (i == n3)
        continue;

      if (7 == z[n3] || 6 == z[n3])   // not sure what else could be there...
        ;
      else
        continue;

      if (-1 == atoms[n3]->formal_charge())    // pathological case of [N+](=[N-])=[N-]
        continue;

      m.set_bond_type_between_atoms (i, n2, TRIPLE_BOND);
      m.set_formal_charge (i, 0);
      m.set_formal_charge (n2, 0);
      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);

      rc++;

      break;
    }
  }

  if (rc)
  {
    _transform_azid.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " charge separated azid groups to neutral\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:azid";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_reverse_nitro (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * o1 = atoms[i];

    const Bond * b1 = o1->item (0);

    if (! b1->is_double_bond())
      continue;

    atom_number_t n = b1->other (i);

    if (7 != z[n])
      continue;
     
    if (3 != ncon[n])
      continue;

    const Atom * na = atoms[n];

    if (5 != na->nbonds())
      continue;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = na->item (j);

      if (! b->is_double_bond())
        continue;

      atom_number_t o2 = b->other (n);

      if (i == o2)
        continue;

      if (8 != z[o2])
        continue;

      m.set_bond_type_between_atoms (i, n, SINGLE_BOND);
      m.set_formal_charge (i, -1);
      m.set_formal_charge (n, +1);

      rc++;
    }
  }

  if (rc)
  {
    _transform_nitro_reverse.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " nitro groups to charge separated\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:NO2+-";
  }

  return rc;
}

int
Chemical_Standardisation::_do_nv5_to_charge_separated(Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

    const Atom * a = atoms[i];

    if (5 != a->nbonds())
      continue;

    if (0 != a->formal_charge())
      continue;

    atom_number_t doubly_bonded_singly_connected = INVALID_ATOM_NUMBER;
    atom_number_t triply_connected_n = INVALID_ATOM_NUMBER;
    bond_type_t bond_to_be_placed = SINGLE_BOND;
    atom_number_t double_bonded_2_connected_n = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (1 == ncon[k] && 8 == z[k] && b->is_double_bond())
        doubly_bonded_singly_connected = k;
      else if (1 == ncon[k] && 7 == z[k] && b->is_triple_bond())
      {
        triply_connected_n = k;
        bond_to_be_placed = DOUBLE_BOND;
      }
      else if (2 == ncon[k] && 7 == z[k] && b->is_double_bond())
        double_bonded_2_connected_n = k;
    }

    if (INVALID_ATOM_NUMBER != doubly_bonded_singly_connected)
      ;
    else if (INVALID_ATOM_NUMBER != triply_connected_n)
      doubly_bonded_singly_connected = triply_connected_n;
    else if (INVALID_ATOM_NUMBER != double_bonded_2_connected_n)
      doubly_bonded_singly_connected = double_bonded_2_connected_n;
    else
      continue;

    m.set_formal_charge(i, 1);
    m.set_formal_charge(doubly_bonded_singly_connected, -1);
    m.set_bond_type_between_atoms(i, doubly_bonded_singly_connected, bond_to_be_placed);
    rc++;
  }

  if (rc)
  {
    _transform_nv5_to_charge_separated.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " 5 valent Nitrogens to charge separated form\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RNV5";
  }

  return rc;
}

int
Chemical_Standardisation::_do_from_mrk_standardisations (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();

  if (0 == current_molecule_data.nneg() && 0 == current_molecule_data.npos())
    return 0;

  int rc = 0;

  if (current_molecule_data.phosphorus() && current_molecule_data.nneg() > 1 && current_molecule_data.npos())   //  [O-]-[P+]-[O-]
  {
    for (int i = 0; i < matoms; i++)
    {
      if (15 != z[i])
        continue;

      if (4 != ncon[i])
        continue;

      const Atom * ai = atoms[i];

      if (ai->formal_charge() < 1)
        continue;

      atom_number_t n1;
      atom_number_t n2;

      if (! two_negatively_charged_connections (i, *ai, atoms, n1, n2))
        continue;

      m.set_formal_charge (i, 0);
      m.set_formal_charge (n1, 0);
      m.set_formal_charge (n2, 0);
      m.set_bond_type_between_atoms (i, n1, DOUBLE_BOND);
      m.set_bond_type_between_atoms (i, n2, DOUBLE_BOND);
      current_molecule_data.change_nneg(-2);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_phosphorus(-1);
      rc++;
    } 
  }

  if ((current_molecule_data.sulphur() || current_molecule_data.phosphorus()) && current_molecule_data.nneg() > 1 && current_molecule_data.npos())    // try to fix [*-]-[S,++]-[*-]
  {
    for (int i = 0; i < matoms; i++)
    {
      if (16 != z[i])
        continue;

      if (4 != ncon[i])
        continue;

      const Atom * ai = atoms[i];

      if (2 != ai->formal_charge())
        continue;

      atom_number_t n1;
      atom_number_t n2;

      if (! two_negatively_charged_connections (i, *ai, atoms, n1, n2))
        continue;

      m.set_formal_charge (i, 0);
      m.set_formal_charge (n1, 0);
      m.set_formal_charge (n2, 0);
      m.set_bond_type_between_atoms (i, n1, DOUBLE_BOND);
      m.set_bond_type_between_atoms (i, n2, DOUBLE_BOND);
      current_molecule_data.change_nneg(-2);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_sulphur(-1);
      rc++;
    }
  }

  if (rc)
  {
    _from_mrk_standardisations.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " moities that had been changed in ->MRK form\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:MRK";
  }

  return rc;
}

/*
  Make sure that the Nitrogen with the Hydrogen is adjacent to the carbon
*/

int
Chemical_Standardisation::_do_tetrazole (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();

  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    if (4 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    const Ring * r = m.ringi (i);

    if (5 != r->number_elements())
      continue;

    if (r->is_fused())
      continue;

    rc += _do_tetrazole (m, *r, current_molecule_data);
  }

  if (rc)
  {
    _transform_tetrazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " tetrazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:tetrazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_triazole (Molecule & m,
                                        const Ring & r,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();

//#define DEBUG_DO_TRIAZOLE
#ifdef DEBUG_DO_TRIAZOLE
  cerr << "Checking possible triazole " << r << endl;
#endif

  int n1_index_in_ring = -1;
  int n2_index_in_ring = -1;
  int n3_index_in_ring = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (7 != z[j])
      continue;

    if (2 != ncon[j])    // If there is a 3 connected Nitrogen, we cannot change the bonding in the ring
      return 0;

    if (n1_index_in_ring < 0)   // there must be exactly three nitrogens in the ring
      n1_index_in_ring = i;
    else if (n2_index_in_ring < 0)
      n2_index_in_ring = i;
    else if (n3_index_in_ring < 0)
      n3_index_in_ring = i;
    else
      return 0;
  }

#ifdef DEBUG_DO_TRIAZOLE
  cerr << " Indices N1 " << n1_index_in_ring << " N2 " << n2_index_in_ring << " and N3 " << n3_index_in_ring << endl;
#endif

  if (n3_index_in_ring < 0)
    return 0;

// We must decide if we have a 1,2,4 or 1,2,3 triazole

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 2 == n3_index_in_ring)
    return _do_123_triazole(m, r, 0, 1, 2, 3, 4, current_molecule_data);

  if (0 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 3, 4, 0, 1, 2, current_molecule_data);

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 4, 0, 1, 2, 3, current_molecule_data);

  if (1 == n1_index_in_ring && 2 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_123_triazole(m, r, 1, 2, 3, 4, 0, current_molecule_data);

  if (2 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 2, 3, 4, 0, 1, current_molecule_data);

// Several more 134 types

  if (0 == n1_index_in_ring && 2 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_134_triazole(m, r, 0, 1, 2, 3, 4, current_molecule_data);

  if (1 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 1, 2, 3, 4, 0, current_molecule_data);

  if (1 == n1_index_in_ring && 2 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 4, 3, 2, 1, 0, current_molecule_data);

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_134_triazole(m, r, 3, 2, 1, 0, 4, current_molecule_data);

  if (0 == n1_index_in_ring && 2 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 2, 3, 4, 0, 1, current_molecule_data);

  cerr << "Unrecognised triazole form n1 = " << n1_index_in_ring << " n2 = " << n2_index_in_ring << " n3 = " << n3_index_in_ring << endl;
  return 0;
}

static int
acyl_group_attached (atom_number_t zatom,
                     const Atom * const * atoms)
{
  const Atom * a = atoms[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(zatom);

    return 8 == atoms[j]->atomic_number();
  }

  return 0;
}

static int
is_cf3 (atom_number_t zatom,
        const Atom * const * atoms)
{
  const Atom * a = atoms[zatom];

  assert (4 == a->ncon());

  int nfluorine = 0;

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (9 != atoms[j]->atomic_number())
      continue;

    nfluorine++;

    if (3 == nfluorine)
      return 1;
  }

  return 0;
}

/*
  Sept 2007. Google search for "electron withdrawing groups" yielded.

  http://www.mhhe.com/physsci/chemistry/carey/student/olc/graphics/carey04oc/ref/ch12substituenteffects.html


*/

static double
compute_electron_donating_power (Molecule & m,
                                 atom_number_t astart,
                                 atom_number_t avoid,
                                 const Atom * const * atoms)
{
#ifdef DEBUG_COMPUTE_ELECTRON_DONATING_POWER
  cerr << "Computing electron donating tendency for atom " << astart << " avoid " << avoid << endl;
#endif

  const Atom * a = atoms[astart];

  atomic_number_t z = a->atomic_number();

  int acon = a->ncon();

  if (1 == acon)
  {
    if (6 == z)
      return 1.3;

    if (7 == z)
      return 3.3;

    if (8 == z)
      return 3.5;

    if (9 == z)
      return -1.4;
    else if (17 == z)
      return -1.3;
    else if (35 == z)
      return -1.2;
    else if (53 == z)
      return -1.1;

    if (16 == z)
      return 3.25;

    return 0.0;   // hmm, what is this?
  }

// We must differentiate NR2 and NO2

  if (7 == z && 3 == acon)    // NR2 or NO2
  {
    if (3 == a->nbonds())   // NR2
      return 3.4;
    else if (5 == a->nbonds())    // presumably NO2
      return -3.6;
  }

  if (7 == z && 4 == acon)
    return -3.5;

  if (16 == z)
  {
    if (4 == acon)   // most likely SO2
      return -3.3;
    else
      return 2.0;   // maybe a thioether?
  }

  if (6 == z)
  {
    if (2 == acon && 4 == a->nbonds())   // cyano or acetylene. Actually reference lists cyano but not acetylene...
      return -3.2;

    if (4 == acon)
    {
      if (is_cf3(astart, atoms))
        return -3.0;
      else
        return 1.3;
    }
  }

// Looking at what is left, it seems that if there is an =O group attached
// we are withdrawing.
// IAW makes up a bunch of other heuristics

  if (acon == a->nbonds())   // saturated
    return static_cast<double>(m.attached_heteroatom_count(astart)) * -0.20;   // IAW
  else if (acyl_group_attached(astart, atoms))
    return -2.0;

// the unsaturated case, which may include aromatic. Cannot compute aromaticity

  int ahc = m.attached_heteroatom_count(astart);

  if (m.multiple_bond_to_heteroatom(astart))   // most likely a Nitrogen
    return 0.57;

  if (ahc > 0)
    return -0.21 * static_cast<double>(ahc);

  return 0.11 * static_cast<double>(acon);

  return 0.0;
}

/*
  The caller must present things with this numbering

       N2
      /  \
     /    \
    N1    N3
    |      |
    |      |
    C5----C4

  We try to put the Hydrogen on N1
*/

static int
place_123_triazole_bonds (Molecule & m,
                          atom_number_t a1,
                          atom_number_t a2,
                          atom_number_t a3,
                          atom_number_t a4,
                          atom_number_t a5,
                          const Atom * const * atoms)
{
  if (2 == atoms[a1]->nbonds())   // already set up properly
    return 0;

  m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
  m.set_bond_type_between_atoms(a2, a3, DOUBLE_BOND);
  m.set_bond_type_between_atoms(a3, a4, SINGLE_BOND);
  m.set_bond_type_between_atoms(a4, a5, DOUBLE_BOND);
  m.set_bond_type_between_atoms(a5, a1, SINGLE_BOND);

  return 1;
}

static atom_number_t
identify_extra_ring_atom (const Atom * a, 
                          atom_number_t zatom,
                          atom_number_t avoid1,
                          atom_number_t avoid2)
{
  if (a->ncon() < 3)
    return INVALID_ATOM_NUMBER;

  for (int i = 0; i < 3; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == avoid1)
      continue;

    if (j == avoid2)
      continue;

    return j;
  }

  return INVALID_ATOM_NUMBER;   // not sure how this could happen
}

/*
  The caller must present things with this numbering

       N2
      /  \
     /    \
    N1    N3
    |      |
    |      |
    C5----C4

  We try to put the Hydrogen on N1 or N3.

  BUT, we run into problems with uniqueness. Therefore, if
  we cannot resolve N1 and N3, we put the Hydrogen on N2
*/

int
Chemical_Standardisation::_do_123_triazole (Molecule & m,
                                        const Ring & r,
                                        int n1_index_in_ring, 
                                        int n2_index_in_ring,
                                        int n3_index_in_ring,
                                        int c4_index_in_ring,
                                        int c5_index_in_ring,
                                        IWStandard_Current_Molecule & current_molecule_data) const
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  atom_number_t n1 = r[n1_index_in_ring];
  atom_number_t n2 = r[n2_index_in_ring];

//#define DEBUG_123_TRIAZOLE
#ifdef DEBUG_123_TRIAZOLE
  cerr << "Checking possible 123 triazole " << r << endl;
  cerr << "nbonds " << atoms[n2]->nbonds() << " ncon " << ncon[n2] << endl;
#endif

  if (3 == ncon[n2])   // rare, cannot change anything
    return 0;

  atom_number_t n3 = r[n3_index_in_ring];

  if (3 == ncon[n1] || 3 == ncon[n3])   // cannot change anything
    return 0;

  atom_number_t c4 = r[c4_index_in_ring];
  atom_number_t c5 = r[c5_index_in_ring];

  if (6 != z[c4])
    return 0;

  if (6 != z[c5])
    return 0;

// do we have unsubstituted 123-triazole

  if (2 == ncon[c4] && 2 == ncon[c5])   // 123-triazole
    return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms);

// Now we need to decide whether to put the H on N1 or N3.
// We know that at least one of C4 and C5 are substituted, 

  atom_number_t cc4 = identify_extra_ring_atom(atoms[c4], c4, n3, c5);
  atom_number_t cc5 = identify_extra_ring_atom(atoms[c5], c5, n1, c4);

  if (INVALID_ATOM_NUMBER == cc4)   // substituted at c5
  {
    double electron_donating_power = compute_electron_donating_power(m, cc5, c5, atoms);
    if (electron_donating_power > 0.0)   // put H on N1
      return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms);
    else    // put H on N3
      return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms);
  }

  if (INVALID_ATOM_NUMBER == cc5)   // sustituted at c4
  {
    double electron_donating_power = compute_electron_donating_power(m, cc4, c4, atoms);
    if (electron_donating_power > 0)   // put H on N3
      return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms);
    else
      return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms);
  }

// Substituted at both positions

  double ed4 = compute_electron_donating_power (m, cc4, c4, atoms);
  double ed5 = compute_electron_donating_power (m, cc5, c5, atoms);

#ifdef DEBUG_123_TRIAZOLE
  cerr << "Atoms C4 " << c4 << " and C5 " << c5 << endl;
  cerr << "ed4 " << ed4 << " ed5 " << ed5 << endl;
#endif

// Wierd compiler but with gcc-4.0.2 on Linux. The comparisons failed!!!!
// C1(=CC=C(F)C=C1)C1=C(N=NN1)C1=CC=NC=C1 
// Restructure code so equality is checked first. Bizzare stuff

  if (fabs(ed4 - ed5) < 1.0e-09)
    return place_123_triazole_bonds(m, n2, n3, c4, c5, n1, atoms);
  else if (ed4 < ed5)   // put Hydrogen on N1
    return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms);
  else
    return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms);
}

/*
  The caller must present things with this numbering

       N1
      /  \
     /    \
    C5    C2
    |      |
    |      |
    N4----N3

  We put double bonds between C2=N3 and N4=C5
*/

int
Chemical_Standardisation::_do_134_triazole (Molecule & m,
                                        const Ring & r,
                                        int n1_index_in_ring, 
                                        int c2_index_in_ring,
                                        int n3_index_in_ring,
                                        int n4_index_in_ring,
                                        int c5_index_in_ring,
                                        IWStandard_Current_Molecule & current_molecule_data) const
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

//#define DEBUG_134_TRIAZOLE
#ifdef DEBUG_134_TRIAZOLE
  cerr << "Checking 134 triazole " << n1_index_in_ring << ',' << c2_index_in_ring << ',' << n3_index_in_ring << ',' <<  n4_index_in_ring << ',' << c5_index_in_ring << endl;
#endif

  atom_number_t n1 = r[n1_index_in_ring];

#ifdef DEBUG_134_TRIAZOLE
  cerr << "What is n1? bonds " << atoms[n1]->nbonds() << " ncon " << ncon[n1] << endl;
#endif

  if (atoms[n1]->nbonds() == ncon[n1])   // already single bonds at N1
    return 0;

  atom_number_t c2 = r[c2_index_in_ring];
  atom_number_t n3 = r[n3_index_in_ring];
  atom_number_t n4 = r[n4_index_in_ring];
  atom_number_t c5 = r[c5_index_in_ring];

#ifdef DEBUG_134_TRIAZOLE
  cerr << "Carbons " << z[c2] << " and " << z[c5] << endl;
#endif

  if (6 != z[c2])
    return 0;

  if (6 != z[c5])
    return 0;

// Beware N1=C(N)NNC1=S MFCD08273670

  if (2 == atoms[n3]->nbonds() && 2 == atoms[n4]->nbonds())
    return 0;

  m.set_bond_type_between_atoms(n1, c5, SINGLE_BOND);
  m.set_bond_type_between_atoms(n1, c2, SINGLE_BOND);
  m.set_bond_type_between_atoms(n3, n4, SINGLE_BOND);
  m.set_bond_type_between_atoms(c2, n3, DOUBLE_BOND);
  m.set_bond_type_between_atoms(n4, c5, DOUBLE_BOND);

  return 1;
}

int
Chemical_Standardisation::_do_triazole (Molecule & m,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  int nr = m.nrings();

  int rc = 0;

#ifdef DEBUG_DO_TRIAZOLE
  cerr << "_possible_triazole " << _possible_triazole << ", nr = " << nr << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    if (3 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    const Ring * r = m.ringi (i);

    if (5 != r->number_elements())
      continue;

    if (r->is_fused())
      continue;

    rc += _do_triazole (m, *r, current_molecule_data);
  }

  if (rc)
  {
    _transform_triazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " triazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:triazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();

  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    if (2 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    const Ring * r = m.ringi (i);

    if (5 != r->number_elements())
      continue;

    rc += _do_imidazole (m, *r, current_molecule_data);
  }

  return rc;
}

static int
expand_shell (const Molecule & m,
              int * visited)
{
#define VISITED_HERE 9

  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (1 != visited[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == visited[k])
      {
        visited[k] = VISITED_HERE;
        rc++;
      }
    }
  }

  if (0 == rc)
    return 0;

  rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (VISITED_HERE == visited[i])
    {
      visited[i] = 1;
      rc += 5 * m.atomic_number(i) + m.ncon(i);   // just some arbitrary thing
    }
  }

  return rc;
}

int 
Chemical_Standardisation::_swap_imidazole (Molecule & m,
                                           atom_number_t n1,
                                           atom_number_t c,
                                           atom_number_t n2) const
{
  m.set_bond_type_between_atoms(n1, c, DOUBLE_BOND);
  m.set_bond_type_between_atoms(c, n2, SINGLE_BOND);
  m.set_implicit_hydrogens_known(n1, 0);
  m.set_implicit_hydrogens_known(n2, 0);

  return 1;
}

int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         const Ring & r,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const Atom * const * atoms = current_molecule_data.atoms();

  int ndxh0 = -1;         // index of nitrogen atom with zero hydrogens
  int ndxh1 = -1;         // index of nitrogen atom with one  hydrogens
  int c31 = -1;         // indices of carbon atoms
  int c32 = -1;
  int c33 = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t a = r[i];

    atomic_number_t za = z[a];

    if (7 == za)
    {
      if (m.hcount(a))
      {
        if (ndxh1 >= 0)
          return 0;

        ndxh1 = i;
      }
      else
      {
        if (ndxh0 >= 0)
          return 0;

        ndxh0 = i;
      }
    }
    else if (6 == za)
    {
      if (c31 < 0)
        c31 = i;
      else if (c32 < 0)
        c32 = i;
      else if (c33 < 0)
        c33 = i;
      else              // Not an imidazole
        return 0;
    }
    else
      return 0;
  }

  if (ndxh0 < 0 || ndxh1 < 0)
    return 0;

  if (c31 < 0 || c32 < 0 || c33 < 0)
    return 0;

// Make sure this is an imidazole - nitrogens separated by one atom

  atom_number_t nh1, c1, nh0, c2, c3;
  if ((ndxh0 + 2) % 5 == ndxh1)
  {
    nh0 = r[ndxh0];
    c1  = r[(ndxh0 + 1) % 5];
    nh1 = r[(ndxh0 + 2) % 5];
    c3  = r[(ndxh0 + 3) % 5];
    c2  = r[(ndxh0 + 4) % 5];
  }
  else if ((ndxh1 + 2) % 5 == ndxh0)
  {
    nh1 = r[ndxh1];
    c1  = r[(ndxh1 + 1) % 5];
    nh0 = r[(ndxh1 + 2) % 5];
    c2  = r[(ndxh1 + 3) % 5];
    c3  = r[(ndxh1 + 4) % 5];
  }
  else         // nitrogens not separated by two bonds
    return 0;

#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Atoms " << nh1 << " " << c1 << " " << nh0 << " " << c2 << " " << c3 << endl;
#endif

// probably not necessary to check the bonds since the ring is aromatic

  const Bond * b = atoms[c2]->bond_to_atom(c3);
  if (! b->is_double_bond())
    return 0;

  b = atoms[c1]->bond_to_atom(nh0);
  if (! b->is_double_bond())
    return 0;

// At this stage we have an imidazole. See if we can resolve it by
// number of connections at c2 and c3

/*             c1
             /    \
            /      \
         nh1        nh0
           |        |
          c3 ------ c2
*/

  assert (m.are_bonded(c2, nh0));
  assert (m.are_bonded(c3, nh1));
  assert (m.are_bonded(c2, c3));

#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Ncon " << atoms[c2]->ncon() << " and " << atoms[c3]->ncon() << endl;
#endif

  if (atoms[c2]->ncon() < atoms[c3]->ncon())   // already as we want it
    return 0;

  if (atoms[c2]->ncon() > atoms[c3]->ncon())
    return _swap_imidazole(m, nh1, c1, nh0);

  if (! r.is_fused())   // cannot expand around the join points
    return 0;

// Not differentiated by the connectivity, start expanding

  int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); iw_auto_array<int> free_tmp(tmp);

  int * tmp1 = tmp;
  int * tmp2 = tmp + matoms;

// Mark all atoms in the ring as visited

  r.set_vector(tmp1, 1);
  r.set_vector(tmp2, 1);

// Mark the ring atoms as off limits for each expansion

  tmp1[c2] = -1;
  tmp2[c3] = -1;

#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Beginning out of ring imidazole detection\n";
#endif

  while (1)
  {
    int e1 = expand_shell(m, tmp1);
    int e2 = expand_shell(m, tmp2);

#ifdef DEBUG_DO_IMIDAZOLE
    cerr << " e1 " << e1 << " and " << e2 << endl;
#endif

    if (e1 > e2)    // correct as is
      return 0;
    else if (e1 < e2)
      return _swap_imidazole(m, nh1, c1, nh0);

    if (0 == e1)    // no more expansion, cannot be resolved
      return 0;
  }

  return 1;
}

/*
  Make sure that the Nitrogen with the Hydrogen is adjacent to the carbon
*/

/*int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         const atomic_number_t * z,
                                         const int * ncon,
                                         Atom ** atoms,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  int nr = m.nrings();

  int rc = 0;

//cerr << "_possible_imidazole " << _possible_imidazole << ", nr = " << nr << endl;

  for (int i = current_molecule_data.possible_imidazole(); i < nr; i++)
  {
    const Ring * r = m.ringi (i);
    if (5 != r->number_elements())
      continue;

    if (r->is_fused())
      continue;

    rc += _do_imidazole (m, *r, z, ncon, atoms);
  }

  if (rc)
  {
    _transform_imidazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " imidazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:imidazole";
  }

  return rc;
}*/

int
Chemical_Standardisation::_do_pyrazole (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  int nr = m.nrings();

  int rc = 0;

//#define DEBUG_DO_PYRAZOLE
#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Processing pyrazoles, nrings " << nr << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    if (2 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    const Ring * r = m.ringi (i);

    if (5 != r->number_elements())
      continue;

//  cerr << "Ring " << (*r) << endl;
//  if (r->is_fused())
//    continue;

    rc += _do_pyrazole (m, *r, current_molecule_data);
  }

  if (rc)
  {
    _transform_pyrazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " pyrazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:pyrazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_tetrazole (Molecule & m,
                                         const Ring & r,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int c = -1;     // index of the Carbon atom
  
  atom_number_t nitrogen_with_one_hydrogen = INVALID_ATOM_NUMBER;
  int nitrogens_with_no_implicit_hydrogens = 0;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    Atom * aj = const_cast<Atom *> (atoms[j]);    // the implicit_hydrogens method is possibly non-const

    if (aj->formal_charge() > 0)   // we may have a negatively charged N in the ring
      continue;

    if (7 == z[j] && 2 == ncon[j])
    {
      if (-1 == aj->formal_charge())
      {
        if (INVALID_ATOM_NUMBER != nitrogen_with_one_hydrogen)
          return 0;

        nitrogen_with_one_hydrogen = j;
        continue;
      }

      int hcount = aj->implicit_hydrogens();
      if (0 == hcount)
      {
        nitrogens_with_no_implicit_hydrogens++;
        continue;
      }

      if (INVALID_ATOM_NUMBER != nitrogen_with_one_hydrogen)
        return 0;

      nitrogen_with_one_hydrogen = j;
    }
    else if (6 == z[j] && 3 == ncon[j] && 4 == aj->nbonds())
    {
      if (c >= 0)     // can be only one Carbon atom in a tetrazole
        return 0;

      c = i;
    }
    else
      return 0;
  }

  if (c < 0 || nitrogen_with_one_hydrogen < 0 || 3 != nitrogens_with_no_implicit_hydrogens)
    return 0;

// Identify the four nitrogens

  atom_number_t n1 = r.next_after_wrap (c);

  if (nitrogen_with_one_hydrogen == n1)
    return 0;

  atom_number_t n2 = r.next_after_wrap (c);
  atom_number_t n3 = r.next_after_wrap (c);
  atom_number_t n4 = r.next_after_wrap (c);

  if (nitrogen_with_one_hydrogen == n4)
    return 0;

// Not in the right form. 

  if (2 == atoms[n3]->nbonds())
    m.set_bond_type_between_atoms (n1, n2, SINGLE_BOND);
  else
    m.set_bond_type_between_atoms (n3, n4, SINGLE_BOND);

  m.set_bond_type_between_atoms (n2, n3, DOUBLE_BOND);

  if (-1 == m.formal_charge (nitrogen_with_one_hydrogen))
    m.set_formal_charge (nitrogen_with_one_hydrogen, 0);

  m.set_implicit_hydrogens_known(nitrogen_with_one_hydrogen, 0);

  return 1;
}

/*
  Make sure the NH is adjacent to the carbon
*/

static int
switch_pyrazole (Molecule & m,
                  atom_number_t c1,
                  atom_number_t n1,
                  atom_number_t n2,
                  atom_number_t c2,
                  atom_number_t c_opposite)
{
  assert (0 == m.hcount(n1));

  m.set_bond_type_between_atoms(n1, c1, SINGLE_BOND);
  m.set_bond_type_between_atoms(c2, c_opposite, SINGLE_BOND);

  m.set_implicit_hydrogens_known(n2, 0);

  m.set_bond_type_between_atoms(c1, c_opposite, DOUBLE_BOND);
  m.set_bond_type_between_atoms(n2, c2, DOUBLE_BOND);

  return 1;
}

/*
  We have a pyrazole and need to know if is part of a fused
  aromatic system
*/

static int
determine_fused_aromatic_pyrazole (const Molecule & m,
                                   const Ring & r,
                                   atom_number_t c1,
                                   atom_number_t c_opposite,
                                   atom_number_t c2,
                                   IWStandard_Current_Molecule & current_molecule_data)
{
  if (! r.is_fused())
    return 0;

  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  if (3 != ncon[c_opposite])   // must be fused in some strange way
    return 0;

  for (int i = 0; i < 3; i++)
  {
    atom_number_t j = atoms[c_opposite]->other(c_opposite, i);

    if (c1 == j || c2 == j)
      continue;

    return atom_is_aromatic[j];
  }

  return 0;   // should never come here
}

/*
  Basically make sure there is a double bond between atoms A1 and A2
*/

static int
switch_fused_pyrazole (Molecule & m,
                       atom_number_t a1,
                       atom_number_t a2)
{
  Toggle_Kekule_Form tkf;

  tkf.set_allow_pyrrole_to_change(1);

  tkf.set_display_error_messages(0);

  int changed;

  if (! tkf.process(m, a1, a2, DOUBLE_BOND, changed))
    return 0;

  return changed;
}

int
Chemical_Standardisation::_do_pyrazole (Molecule & m,
                                        const Ring & r,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Checking pyrazole " << r << endl;
#endif

  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int n1_index_in_ring = -1;
  int n2_index_in_ring = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (2 != ncon[j])
      continue;

    Atom * aj = const_cast<Atom *> (atoms[j]);    // the implicit_hydrogens method is possibly non-const

    if (7 != aj->atomic_number())
      continue;

    if (n1_index_in_ring < 0)   // there must be exactly two nitrogens in the ring
      n1_index_in_ring = i;
    else if (n2_index_in_ring < 0)
      n2_index_in_ring = i;
    else
      return 0;
  }

  if (n2_index_in_ring < 0)
    return 0;

  if (n1_index_in_ring > n2_index_in_ring)
  {
    int tmp = n1_index_in_ring;
    n1_index_in_ring = n2_index_in_ring;
    n2_index_in_ring = tmp;
  }
//  std::swap(n1_index_in_ring, n2_index_in_ring);

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Got two nitrogens, n1 " << n1_index_in_ring << " and n2 " << n2_index_in_ring << endl;
#endif

// Thet two nitrogens must be adjacent. Identify the adjoining carbons

  int c1_index = -1;;
  int c2_index = -1;;
  int c_opposite_index = -1;

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring)
  {
    c1_index = 4;
    c2_index = 2;
    c_opposite_index = 3;
  }
  else if (1 == n1_index_in_ring && 2 == n2_index_in_ring)
  {
    c1_index = 0;
    c2_index = 3;
    c_opposite_index = 4;
  }
  else if (2 == n1_index_in_ring && 3 == n2_index_in_ring)
  {
    c1_index = 1;
    c2_index = 4;
    c_opposite_index = 0;
  }
  else if (3 == n1_index_in_ring && 4 == n2_index_in_ring)
  {
    c1_index = 2;
    c2_index = 0;
    c_opposite_index = 1;
  }
  else if (0 == n1_index_in_ring && 4 == n2_index_in_ring)
  {
    c1_index = 1;
    c2_index = 3;
    c_opposite_index = 2;
  }
  else
    return 0;

// Now we need to differentiate the two nitrogens...

  atom_number_t c1 = r[c1_index];
  atom_number_t n1 = r[n1_index_in_ring];
  atom_number_t n2 = r[n2_index_in_ring];
  atom_number_t c2 = r[c2_index];

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole atoms c1 " << c2 << " n1 " << n1 << " n2 " << n2 << " c2 " << c2 << endl;
#endif

  if (6 != atoms[c1]->atomic_number())
    return 0;

  if (6 != atoms[c2]->atomic_number())
    return 0;

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "First two carbons OK\n";
#endif

  atom_number_t c_opposite = r[c_opposite_index];

  if (6 != atoms[c_opposite]->atomic_number())
    return 0;

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole atoms c1 " << c2 << " n1 " << n1 << " n2 " << n2 << " c2 " << c2 << " opposite " << c_opposite << endl;
#endif

  int h1 = m.hcount(n1);
  int h2 = m.hcount(n2);

  if (0 == h1 && 0 == h2)
    return 0;

  if (h1 && h2)
    return 0;

  atom_number_t n_atom_now_with_hydrogen;

  if (h1)
  {
    if (3 != atoms[n2]->nbonds())
      return 0;

    n_atom_now_with_hydrogen = n1;
  }
  else
  {
    if (3 != atoms[n1]->nbonds())
      return 0;

    n_atom_now_with_hydrogen = n2;
  }

#ifdef DEBUG_DO_PYRAZOLE
  cerr << n_atom_now_with_hydrogen << " has the hydrogen, connections " << ncon[c1] << " and " << ncon[c2] << endl;
#endif

  if (2 == ncon[c1] && 2 == ncon[c2])    // cannot do anything, symmetric
    return 0;

// Many pyrazoles are fused, but to an aliphatic ring, so 
// they do not need to be treated by Toggle_Kekule_Form

  int fused_aromatic = determine_fused_aromatic_pyrazole(m, r, c1, c_opposite, c2, current_molecule_data);

// See if we can resolve things by connectivity

  if (ncon[c1] < ncon[c2])
  {
    if (n2 == n_atom_now_with_hydrogen)
      return 0;

    if (fused_aromatic)
      return switch_fused_pyrazole(m, c2, n2);
    else
      return switch_pyrazole(m, c2, n2, n1, c1, c_opposite);
  }
  else if (ncon[c1] > ncon[c2])
  {
    if (n1 == n_atom_now_with_hydrogen)
      return 0;

    if (fused_aromatic)
      return switch_fused_pyrazole(m, c1, n1);
    else
      return switch_pyrazole(m, c1, n1, n2, c2, c_opposite);
  }

// Both adjacent carbons have 3 connections. Resolve by shells

  int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); iw_auto_array<int> free_tmp(tmp);
  int * tmp1 = tmp;
  int * tmp2 = tmp1 + matoms;

  r.set_vector(tmp1, 1);
  r.set_vector(tmp2, 1);
  tmp1[n2] = -1;
  tmp2[n1] = -1;

  while (1)
  {
    int e1 = expand_shell(m, tmp1);
    int e2 = expand_shell(m, tmp2);

    if (e1 == e2)
    {
      if (0 == e1)   // done, cannot be resolved
        return 0;

      continue;
    }

    if (e1 > e2)
      return 0;

    if (fused_aromatic)
    {
      if (n_atom_now_with_hydrogen == n1)
        return 0;

      return switch_fused_pyrazole(m, c2, n2);
    }
    else if (n_atom_now_with_hydrogen == n1)
      return switch_pyrazole(m, c2, n2, n1, c1, c_opposite);
    else
      return switch_pyrazole(m, c1, n1, n2, c2, c_opposite);
  }

  return 0;   // cannot figure out what to do
}


/*
  Take N=N to [N+]-[N-]

  We also do CN(C)(C)=N
*/

int
Chemical_Standardisation::_do_transform_back_to_nplus_nminus (Molecule & m,
                                          IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * ai = atoms[i];

    if (0 != ai->formal_charge())
      continue;

    if (5 != ai->nbonds())
      continue;

    atom_number_t doubly_bonded_nitrogen = INVALID_ATOM_NUMBER;
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item (j);
      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other (i);
      if (7 != z[k])
        continue;

      if (0 != atoms[k]->formal_charge())
        continue;

      doubly_bonded_nitrogen = k;
      break;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_nitrogen)
      continue;

    m.set_formal_charge (i, 1);
    m.set_formal_charge (doubly_bonded_nitrogen, -1);
    m.set_bond_type_between_atoms (i, doubly_bonded_nitrogen, SINGLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_back_to_nplus_nminus.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " =N=N- to =[N+]-[N-]-\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RN+N-";
  }

  return rc;
}

/*
  Convert N#N=N- to [N-]=[N+]=N-
*/

int 
Chemical_Standardisation::_do_transform_azid_to_charge_separated (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (0 != a->formal_charge())
      continue;

    const Bond * b = a->item(0);

    if (! b->is_triple_bond())
      continue;

    atom_number_t n2 = b->other(i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * a2 = atoms[n2];

    if (5 != a2->nbonds())
      continue;

    atom_number_t n3;

    if (i == a2->other(n2, 0))
      n3 = a2->other(n2, 1);
    else
      n3 = a2->other(n2, 0);

    if (7 != z[n3])
      continue;

    if (0 != atoms[n3]->formal_charge())
      continue;

    m.set_formal_charge(i, -1);
    m.set_formal_charge(n2, 1);
    m.set_bond_type_between_atoms(i, n2, DOUBLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_to_charge_separated_azid.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " N#N=N- to N=[N+]=[N-]-\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RAZID";
  }

  return rc;
}

/*
  Change -N=C(-[OH])-N to -N-C(=O)-N
*/

int
Chemical_Standardisation::_do_transform_misdrawn_urea (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (1 != ai->nbonds())
      continue;

    atom_number_t c = ai->other (i, 0);

    if (3 != ncon[c])
      continue;

    if (atom_is_aromatic[c])
      continue;

    const Atom * ac = atoms[c];

    if (4 != ac->nbonds())
      continue;

    atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[c]; j++)
    {
      const Bond * b = ac->item (j);

      atom_number_t k = b->other (c);

      if (7 != z[k])
        continue;

// don't change N1C(=O)C=CC2=C1C=CC=C2.C(=O)(O)/C=C/C(=O)O PTNT74937358

      if (b->is_double_bond())
      {
        if (! b->part_of_cis_trans_grouping())
          doubly_bonded_nitrogen = k;
      }
      else
        singly_bonded_nitrogen = k;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_nitrogen || INVALID_ATOM_NUMBER == singly_bonded_nitrogen)
      continue;

    m.set_bond_type_between_atoms (i, c, DOUBLE_BOND);
    m.set_bond_type_between_atoms (c, doubly_bonded_nitrogen, SINGLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_misdrawn_urea.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " misdrawn ureas\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:MSDUR";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_implicit_hydrogen_known_errors (Molecule & m,
                                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    Atom * ai = const_cast<Atom *>(atoms[i]);

    if (! ai->implicit_hydrogens_known())
      continue;

    if (ai->valence_ok())
      continue;

//  we have an atom with implicit hydrogens known, and an invalid valence.

    int ih = ai->implicit_hydrogens();

    if (ncon[i] >= ai->element()->normal_valence() + ai->formal_charge())
      continue;

    int hdiff = ncon[i] - ai->element()->normal_valence() - ai->formal_charge();

//  cerr << " ih = " << ih << " hdiff " << hdiff << endl;

    ai->set_implicit_hydrogens (ih - hdiff, 1);
    rc++;
  }

  return rc;
}

static int
single_bond_to_oxygen (const Molecule & m,
                       atom_number_t n)
{
  const Atom * a = m.atomi(n);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    int o = b->other(n);

    if (8 == m.atomic_number(o))
      return 1;
  }

  return 0;
}

//#define DEBUG_DO_TRANSFORM_LACTIM

/*
  We have a molecule in the form
 ............
  that is, no double bond between the N and the carbon attached to the oxygen
*/

int
Chemical_Standardisation::_toggle_kekule_forms_to_lactim_form (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) const
{
//cerr << "Into _toggle_kekule_forms_to_lactim_form, smiles '" << m.smiles() << "'\n";
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Bond * b = atoms[i]->item(0);

    if (! b->is_single_bond())
      continue;

    atom_number_t c = b->other(i);

    if (6 != z[c])    // must be carbon
      continue;

    if (0 == ring_membership[c])
      continue;

    if (! current_molecule_data.atom_is_aromatic()[c])
      continue;

    const Atom * ac = atoms[c];

    if (3 != ac->ncon())
      continue;

    atom_number_t single_bond_to_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = ac->item(j);

      if (! b->is_single_bond())
        continue;

      atom_number_t n = b->other(c);

      if (i == n)
        continue;

      if (7 != z[n])
        continue;

      if (3 == ncon[n])   // we need to put a double bond there
        continue;

      single_bond_to_nitrogen = n;
      break;
    }

    if (INVALID_ATOM_NUMBER == single_bond_to_nitrogen)
      continue;

    Toggle_Kekule_Form tkf;
    tkf.set_display_error_messages(0);

    int changed;
    tkf.process(m, single_bond_to_nitrogen, c, DOUBLE_BOND, changed);

//  cerr << "After toggling Kekule form '" << m.smiles() << endl;

    if (changed)
      current_molecule_data.initialise(m);   // need to reset everything
  }

  return 1;
}

int
Chemical_Standardisation::_do_transform_lactim (Molecule & m,
                                    int ring_version,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
#ifdef DEBUG_DO_TRANSFORM_LACTIM
  cerr << "Into _do_transform_lactim, rv " << ring_version << endl;
#endif

  if (ring_version && ! _toggle_kekule_forms_to_lactim_form (m, current_molecule_data))
    return 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * o = atoms[i];

    if (1 != o->nbonds())
      continue;

    atom_number_t c = o->other(i, 0);

    if (6 != z[c])
    {
      current_molecule_data.change_possible_lactim(-1);
      continue;
    }

#ifdef DEBUG_DO_TRANSFORM_LACTIM
    cerr << "Atom " << i << " ring_version " << ring_version << " ring_membership " << ring_membership[c] << " possible " << current_molecule_data.possible_lactim() << endl;
#endif

    if ((ring_version && 0 == ring_membership[c]) ||
        (0 == ring_version && ring_membership[c]))
      continue;    // do not decrement _possible_lactim, may be ring/non-ring form waiting to be processed

    const Atom * ac = atoms[c];

    int ccon = ac->ncon();

#ifdef DEBUG_DO_TRANSFORM_LACTIM
    cerr << "Carbon atom has " << ccon << " connections\n";
#endif

    if (1 == ccon || ccon + 1 != ac->nbonds())
    {
      current_molecule_data.change_possible_lactim(-1);
      continue;
    }

    atom_number_t n = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ccon; j++)
    {
      const Bond * b = ac->item(j);
      if (! b->is_double_bond())
        continue;

      if (b->part_of_cis_trans_grouping())    // C(=NNC(N)=N)(C=CC1=CC=C(O1)N(=O)=O)C(=O)NO MDDR238982 - shonw without cis trans bonds
        continue;

      atom_number_t k = b->other(c);

      if (7 != z[k])
      {
        current_molecule_data.change_possible_lactim(-1);
        break;
      }

      if (INVALID_ATOM_NUMBER != n)   // cannot process N-C(=O)-N, possibly ambiguous. Fix sometime...
        continue;

      n = k;
    }

    if (INVALID_ATOM_NUMBER == n)
      continue;

    if (single_bond_to_oxygen(m, n))
    {
      current_molecule_data.change_possible_lactim(-1);
      continue;
    }

    m.set_bond_type_between_atoms(c, n, SINGLE_BOND);
    m.set_bond_type_between_atoms(i, c, DOUBLE_BOND);
    rc++;
    current_molecule_data.change_possible_lactim(-1);
    if (0 == current_molecule_data.possible_lactim())
      break;
  }

  if (0 == rc)
    ;
  else if (ring_version)
  {
    _transform_lactim_lactam_ring.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " lactim->lactam ring\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:LTLTR";
  }
  else
  {
    _transform_lactim_lactam.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " lactim->lactam\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:LTLT";
  }

  return rc;
}

int
Chemical_Standardisation::activate_from_corina_transformations()
{
  _transform_nitro.activate();
  _transform_azid.activate();
  _transform_nplus_ominus.activate();
  _active = 1;

  return 1;
}

/*
  Corina puts a negative charge on -P(=O)(=O)-
*/

#ifdef PMINUS

int
Chemical_Standardisation::_do_p_minus (Molecule & m,
                                       IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;
  }
}
#endif

IWStandard_Current_Molecule::IWStandard_Current_Molecule()
{
  _atomic_number = NULL;
  _ncon = NULL;
  _ring_membership = NULL;
  _atom_is_aromatic = NULL;
  _ring_is_aromatic = NULL;
  _atom = NULL;
  _ring_nitrogen_count = NULL;

  _npos = 0;
  _nneg = 0;
  _ominus = 0;
  _sminus = 0;
  _nplus = 0;
  _cminus = 0;
  _phosphorus = 0;
  _sulphur = 0;
  _isolated_metal = 0;
  _isolated_halogen = 0;
  _singly_connected_metal = 0;
  _possible_guanidine = 0;
  _phosphorus = 0;
  _explicit_hydrogen_count = 0;
  _possible_valence_errors = 0;

  _nitrogens = 0;
  _oxygens = 0;

  _possible_lactim = 0;

  return;
}

IWStandard_Current_Molecule::~IWStandard_Current_Molecule ()
{
  if (NULL != _atomic_number)
    delete [] _atomic_number;

  if (NULL != _ncon)
    delete [] _ncon;

  if (NULL != _ring_membership)
    delete [] _ring_membership;

  if (NULL != _atom_is_aromatic)
    delete [] _atom_is_aromatic;

  if (NULL != _ring_is_aromatic)
    delete [] _ring_is_aromatic;

  if (NULL != _atom)
    delete [] _atom;

  if (NULL != _ring_nitrogen_count)
    delete [] _ring_nitrogen_count;

  return;
}

int
IWStandard_Current_Molecule::initialise (Molecule & m)
{
  _matoms = m.natoms();

  if (0 == _matoms)
    return 0;

  _atom = new const Atom * [_matoms];

  m.atoms(_atom);

  _atomic_number = new atomic_number_t[_matoms];
  _ncon = new int[_matoms];

  _nrings = m.nrings();

  _atom_is_aromatic = new_int(_matoms);

  if (_nrings)
  {
    _ring_membership = new int[_matoms];
    m.ring_membership(_ring_membership);
    _ring_is_aromatic = new_int(_nrings);
    _ring_nitrogen_count = new_int(_nrings);
    m.compute_aromaticity_if_needed();

    for (int i = 0; i < _nrings; i++)
    {
      const Ring * ri = m.ringi(i);
      if (ri->is_aromatic())
      {
        _ring_is_aromatic[i] = 1;
        ri->set_vector(_atom_is_aromatic, 1);
      }
    }
  }
  else
    _ring_membership = new_int(_matoms);

  for (int i = 0; i < _matoms; i++)
  {
    Atom * ai = const_cast<Atom *>(_atom[i]);

    if (! ai->valence_ok())
      _possible_valence_errors++;

    formal_charge_t fc = ai->formal_charge();

    _ncon[i] = ai->ncon();

    atomic_number_t z = ai->atomic_number();

    _atomic_number[i] = z;

    if (6 == z)
    {
      if (fc < 0)
        _cminus++;
    }
    else if (7 == z)
    {
      _nitrogens++;
      if (fc > 0)
        _nplus++;
    }
    else if (8 == z)
    {
      _oxygens++;
      if (ai->formal_charge() < 0)
        _ominus++;
      else if (1 == _ncon[i] && 2 == ai->nbonds() && 1 == _ring_membership[ai->other(i, 0)])
        _possible_lactim++;
    }
    else if (16 == z)
    {
      _sulphur++;
      if (fc < 0)
        _sminus++;
    }
    else if (15 == z)
      _phosphorus++;
    else if (ai->element()->is_halogen() && 0 == _ncon[i])
      _isolated_halogen++;
    else if (ai->element()->is_metal())
    {
      if (0 == _ncon[i])
        _isolated_metal++;
      else if (1 == _ncon[i])
        _singly_connected_metal++;
    }
    else if (1 == z)
      _explicit_hydrogen_count++;

    if (0 == fc)
      continue;

    if (fc < 0)
      _nneg++;
    else if (fc > 0)
      _npos++;
  }

  if (_nitrogens > _possible_lactim)
    _possible_lactim = _nitrogens;

  _possible_guanidine = _nitrogens / 3;

  if (_nitrogens >= 2 && _nrings)
  {
    for (int i = 0; i < _nrings; i++)
    {
      const Ring * ri = m.ringi(i);

      if (5 != ri->number_elements())
        continue;

      int nitrogens = 0;
      for (int j = 0; j < 5; j++)
      {
        atom_number_t k = ri->item(j);
        if (7 == _atomic_number[k])
          nitrogens++;
      }

      if (nitrogens < 2)
        continue;

      _ring_nitrogen_count[i] = nitrogens;
    }
  }

  return 1;
}

int
Chemical_Standardisation::_processing_needed (const IWStandard_Current_Molecule & current_molecule_data) const
{
  if (current_molecule_data.nneg() || current_molecule_data.npos() || 
     (current_molecule_data.singly_connected_metal() && _transform_covalent_metals.active()) ||
     (current_molecule_data.isolated_metal() && _transform_single_atom_ions.active()) ||
     (current_molecule_data.isolated_halogen() && _transform_single_atom_ions.active()) ||
      current_molecule_data.possible_guanidine() ||
      current_molecule_data.explicit_hydrogen_count() ||
      current_molecule_data.aromatic_rings_with_multiple_nitrogens() > 1 ||
      (current_molecule_data.nitrogens ()> 1 && _transform_back_to_nplus_nminus.active()) ||
      (current_molecule_data.nitrogens() > 1 && current_molecule_data.oxygens() > 0 && _transform_misdrawn_urea.active()) ||
      current_molecule_data.possible_lactim() ||
      current_molecule_data.possible_valence_errors())

    return 1;

  return 0;
}

int
IWStandard_Current_Molecule::aromatic_rings_with_multiple_nitrogens () const
{
  if (0 == _nrings)
    return 0;

  int rc = 0;

  for (int i = 0; i < _nrings; i++)
  {
    if (! _ring_is_aromatic[i])
      continue;

    if (_ring_nitrogen_count[i] > 1)
      rc++;
  }

  return rc;
}
