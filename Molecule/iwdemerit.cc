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
#include <memory>
using namespace std;

#if defined (IWSGI) || (__GNUC_MINOR__ == 95)
#include <strstream>
#else
#include <sstream>
#endif

#include <assert.h>
#include <math.h>
#include <errno.h>

#include "iwconfig.h"

#include "cmdline.h"

#include "rwsubstructure.h"
#include "target.h"
#include "smiles.h"
#include "qry_wstats.h"
#include "aromatic.h"
#include "output.h"
#include "istream_and_type.h"
#include "charge_assigner.h"
#include "misc2.h"
#include "rmele.h"

#include "substructure_demerits.h"
#include "demerit.h"


//#define USE_IWMALLOC

#ifdef USE_IWMALLOC
#include "iwmalloc.h"
#endif

static int verbose = 0;     // visible in substructure_demerits

static int keep_going_after_rejection = 0;

static int molecules_read = 0;

const char * prog_name = NULL;

static int molecules_receiving_demerits = 0;
static int molecules_rejected = 0;

static Elements_to_Remove elements_to_remove;

static Molecule_Output_Object stream_for_non_rejected_molecules;
static Molecule_Output_Object stream_for_rejected_molecules;

static IWString stem_for_demerit_file;

static int do_hard_coded_substructure_queries = 1;

static int append_demerit_text_to_name = 0;

static int write_rejection_reason_like_tsubstructure = 0;

static int write_rejection_reasons_separately = 0;

static int skip_molecules_with_abnormal_valences;
static int molecules_with_abnormal_valences = 0;

static resizable_array_p<Substructure_Hit_Statistics> queries;

/*
  For molecules that lie between the hard and low atom count cutoffs, we apply demerits
*/

static int soft_lower_atom_count_cutoff = 0;
static int hard_lower_atom_count_cutoff = 0;
static int lower_atom_count_demerit = 100;

static int soft_upper_atom_count_cutoff = 0;
static int hard_upper_atom_count_cutoff = 0;
static int upper_atom_count_demerit = 100;

/*
  We reject if the atom contains only C, H or S
  On 20 Nov, we expanded that to be C, H, S, or halogen
  We reject of we see any strange element types.
*/

static int atom_types_count = 0;
static int csxh_count = 0;

static void
do_atom_count_demerits (Molecule & m,
                        Demerit & demerit)
{
  int nf = m.number_fragments ();

  int matoms;
  if (1 == nf)
    matoms = m.natoms ();
  else
  {
    matoms = 0;
    for (int i = 0; i < nf; i++)
    {
      int f = m.atoms_in_fragment (i);
      if (f > matoms)
        matoms = f;
    }
  }

  if (hard_lower_atom_count_cutoff > 0 && matoms <= hard_lower_atom_count_cutoff)
  {
    demerit.extra (lower_atom_count_demerit, "too_few_atoms");

    return;
  }

// Feb 2005. Heuristic for adding demerits beyond the hard atom count cutoff

  if (hard_upper_atom_count_cutoff > 0 && matoms >= hard_upper_atom_count_cutoff)
  {
    demerit.extra (rejection_threshold() + 6 * (matoms - hard_upper_atom_count_cutoff), "too_many_atoms");

//  cerr << hard_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff, matoms = " << matoms << endl;
    return;
  }

  if (matoms > hard_lower_atom_count_cutoff && matoms < soft_lower_atom_count_cutoff)
  {
    float r = static_cast<float> (soft_lower_atom_count_cutoff - matoms) / static_cast<float> (soft_lower_atom_count_cutoff - hard_lower_atom_count_cutoff);
    int d = static_cast<int> (lower_atom_count_demerit * r);

    if (0 == d)
      d = 1;

    demerit.extra (d, "too_few_atoms");

    return;
  }

  if (matoms > soft_upper_atom_count_cutoff && matoms < hard_upper_atom_count_cutoff)
  {
    float r = static_cast<float> (matoms - soft_upper_atom_count_cutoff) / static_cast<float> (hard_upper_atom_count_cutoff - soft_upper_atom_count_cutoff);
    int d = static_cast<int> (upper_atom_count_demerit * r);
//  cerr << "too_many_atoms: " << soft_upper_atom_count_cutoff << " soft_upper_atom_count_cutoff, matoms = " << matoms << ", d = " << d << endl;
//  cerr << "upper_atom_count_demerit " << upper_atom_count_demerit << endl;
//  cerr << "too_many_atoms: soft_upper_atom_count_cutoff " << soft_upper_atom_count_cutoff << " hard_upper_atom_count_cutoff " << hard_upper_atom_count_cutoff << " matoms " << matoms << " D = " << d << endl;
    if (0 == d)
      d = 1;

    demerit.extra (d, "too_many_atoms");

    return;
  }

  return;
}


static void
run_a_set_of_queries(Molecule_to_Match & target,
                     Demerit & demerit,
                     resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  int nqueries = queries.number_elements ();

  for (int i = 0; i < nqueries; i++)
  {
    Substructure_Hit_Statistics * q = queries[i];

    int nhits = q->substructure_search (target);
    if (0 == nhits)
      continue;

    if (verbose > 1)
      cerr << nhits << " matches to query " << i << ' ' << q->comment () << endl;

    double d;
    (void) q->numeric_value (d);
    
    d = d * nhits;

    int intd = int (d * 1.0001);

    if (intd >= rejection_threshold ())
    {
      demerit.reject (q->comment ());
      if (0 == keep_going_after_rejection)
        return;
    }
    else
      demerit.extra(intd, q->comment ());

    if (demerit.rejected() && 0 == keep_going_after_rejection)
      return;
  }

  return;
}

static void
iwdemerit (Molecule & m,
           resizable_array_p<Substructure_Hit_Statistics> & q1,
           resizable_array_p<Substructure_Hit_Statistics> & q2,
           Demerit & demerit)
{
  if (soft_lower_atom_count_cutoff > 0 || soft_upper_atom_count_cutoff > 0)
  {
    do_atom_count_demerits (m, demerit);
    if (demerit.rejected () && 0 == keep_going_after_rejection)
      return;
  }

  if (skip_molecules_with_abnormal_valences && ! m.valence_ok ())
  {
    demerit.reject ("valence");
    molecules_with_abnormal_valences++;
    if (0 == keep_going_after_rejection)
      return;
  }

  if (do_hard_coded_substructure_queries)
  {
    substructure_demerits::hard_coded_queries (m, demerit);
//  cerr << "After hard coded queries, score is " << demerit.score() << endl;

    if (demerit.rejected () && 0 == keep_going_after_rejection)
      return;
  }

  Molecule_to_Match target (&m);

  if (q1.number_elements())
  {
    run_a_set_of_queries(target, demerit, q1);

//  cerr << "After command line queries, score is " << demerit.score() << " rej? " << demerit.rejected() << endl;
    if (demerit.rejected())
      return;
  }

  if (q2.number_elements())
  {
    m.reduce_to_largest_fragment_carefully();
    Molecule_to_Match target(&m);

    run_a_set_of_queries(target, demerit, q2);

    if (demerit.rejected())
      return;
  }

  return;
}

/*
  In jul 97, people wanted to know how many molecules received multiple
  demerits, and how many were rejected as a result of multiple demerits
*/

static extending_resizable_array<int> demerits_per_molecule;
static extending_resizable_array<int> demerits_per_rejected_molecule;

/*
  Later in Jul 97 they wanted a list of molecules for which there were
  rejected as a result of multiple different types of demerit
*/

static Molecule_Output_Object stream_for_multiple_demerits;

static int
do_append_demerit_text_to_name (Molecule & m,
                                const Demerit & demerit)
{
  IWString new_name (m.molecule_name ());

  if (write_rejection_reasons_separately && write_rejection_reason_like_tsubstructure)
    new_name << ' ' << demerit.types ();
  else if (write_rejection_reason_like_tsubstructure)
    new_name << " (" << demerit.number_different_demerits_applied () << " matches to 'D" << demerit.score () << ' ' << demerit.types () << "')";
  else
    new_name << " : D(" << demerit.score () << ") " << demerit.types ();

  m.set_name (new_name);

  return 1;
}

static int
iwdemerit (Molecule & m,
           resizable_array_p<Substructure_Hit_Statistics> & q1,
           resizable_array_p<Substructure_Hit_Statistics> & q2,
           ofstream & output)
{
  elements_to_remove.process (m);

  Demerit demerit;

  iwdemerit (m, q1, q2, demerit);

  if (demerit.rejected ())
    demerits_per_rejected_molecule[demerit.number_different_demerits_applied ()]++;
  if (demerit.score ())
    demerits_per_molecule[demerit.number_different_demerits_applied ()]++;

  if (demerit.score () && append_demerit_text_to_name)
    do_append_demerit_text_to_name (m, demerit);

//cerr << m.name() << " rejected? " << demerit.rejected () << " stream is " << stream_for_non_rejected_molecules.active () << endl;

  if (demerit.rejected ())
  {
    if (stream_for_rejected_molecules.number_elements ())
      stream_for_rejected_molecules.write (&m);
    if (stream_for_multiple_demerits.number_elements () && ! demerit.rejected_by_single_rule ())
      stream_for_multiple_demerits.write (&m);
  }
  else if (stream_for_non_rejected_molecules.active ())
  {
    stream_for_non_rejected_molecules.write (&m);
  }

  if (verbose > 1 && demerit.score ())
    demerit.debug_print (cerr);

  if (! output.rdbuf ()->is_open ())
    return 1;

  output << "$SMI<" << m.smiles () << ">\n";
  output << "PCN<" << m.name () << ">\n";

  if (demerit.score ())
  {
    molecules_receiving_demerits++;
    if (demerit.rejected ())
      molecules_rejected++;
    demerit.write_in_tdt_form (output);
  }

  output << "|\n";
  output.flush ();

  return output.good ();
}

static void
preprocess (Molecule & m)
{
  m.remove_all (1);

  return;
}

static int
iwdemerit (data_source_and_type<Molecule> & input,
           resizable_array_p<Substructure_Hit_Statistics> & q1,
           resizable_array_p<Substructure_Hit_Statistics> & q2,
           ofstream & output)
{
  assert (input.good ());

  Molecule * m;
  while (NULL != (m = input.next_molecule ()))
  {
    if (verbose > 1)
      cerr << input.molecules_read () <<  " processing '" << m->name () << "'\n";

    molecules_read++;

    auto_ptr<Molecule> free_m (m);

    preprocess (*m);

    if (! iwdemerit (*m, q1, q2, output))
      return 0;

#ifdef USE_IWMALLOC
    terse_malloc_status (stderr);
#endif
  }

  if (verbose)
  {
    cerr << input.molecules_read () << " molecules read\n";
    elements_to_remove.report (cerr);
    cerr << molecules_receiving_demerits << " molecules were assigned demerits\n";
    cerr << molecules_rejected << " molecules were rejected\n";

    cerr << "Details on queries\n";
    cerr << "atom_types_count = " << atom_types_count << endl;
    cerr << "only C S X or H = " << csxh_count << endl;

    if (do_hard_coded_substructure_queries)
      substructure_demerits::hard_coded_queries_statistics (cerr);

    for (int i = 0; i < queries.number_elements (); i++)
    {
      cerr << "Results of query " << i << endl;
      queries[i]->report (cerr, verbose);
    }
  }

  return 1;
}

static int
iwdemerit (const char * fname, int input_type,
           resizable_array_p<Substructure_Hit_Statistics> & q1,
           resizable_array_p<Substructure_Hit_Statistics> & q2)
{
  if (0 == input_type)
  {
    input_type = discern_file_type_from_name (fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input (input_type, fname);
  if (! input.ok ())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  ofstream output;
  
  if (stem_for_demerit_file.length ())
  {
    IWString output_file_name = stem_for_demerit_file;

    change_suffix (output_file_name, "demerit");

    output.open (output_file_name.null_terminated_chars ());
    if (! output.good ())
    {
      cerr << prog_name << ": cannot open '" << output_file_name << "' for output\n";
      return 2;
    }
  }

  return iwdemerit (input, q1, q2, output);
}

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " <<__TIME__ << endl;

  cerr << "Usage : " << prog_name << " options file1 file2 file3 ....\n";
  cerr << "  -q <file>      specify substructure query file\n";
  cerr << "  -R <stem>      write rejected structures to 'stem.otype'\n";
  cerr << "  -G <stem>      write non rejected (good) structures to 'stem.otype'\n";
  cerr << "  -S <stem>      specify stem for .demerit file\n";
  cerr << "  -M <stem>      write molecules rejected by cumulative demerits to <stem>\n";
  cerr << "  -k             check all criteria even if a molecule is already rejected\n";
  cerr << "  -c <number>    lower atom count cutoff\n";
  cerr << "  -O hard        skip all the hard coded substructure queries\n";
  cerr << "  -V             reject molecules containing unusual valences\n";
  cerr << "  -t             append demerit text to molecule names\n";
  cerr << "  -f <n>         molecules are rejected when they have <n> demerits\n";
  cerr << "  -x             atom count demerits remain scaled to 100\n";
  cerr << "  -E <symbol>    create element with symbol\n";
  cerr << "  -X <symbol>    before processing, delete all <symbol> atoms\n";
  cerr << "  -C <file>      control file for what demerits to apply\n";
  cerr << "  -r             only do rejection rules - no demerits\n";
  cerr << "  -N ...         charge assigner specifications, enter '-N help' for info\n";
  cerr << "  -o <type>      file type for structures written\n";
  cerr << "  -i <type>      specify input file type\n";
  display_standard_aromaticity_options (cerr);
#ifdef USE_IWMALLOC
  cerr << "  -d <block>     die when block <block> is allocated\n";
#endif
  cerr << "  -v             verbose output\n";

  exit (rc);
}

static void
separate_depending_on_fragment_match (resizable_array_p<Substructure_Hit_Statistics> & queries,
                                      resizable_array_p<Substructure_Hit_Statistics> & q1,
                                      resizable_array_p<Substructure_Hit_Statistics> & q2)
{
  int n = queries.number_elements();

  for (int i = n - 1; i >= 0; i--)
  {
    Substructure_Hit_Statistics * q = queries.remove_no_delete(i);

    if (q->only_keep_matches_in_largest_fragment())
      q2.add(q);
    else
      q1.add(q);
  }

  if (verbose > 1)
    cerr << "Separated queries into " << q1.number_elements() << " and " << q2.number_elements() << " queries\n";

  return;
}

int
iwdemerit (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "M:VX:tA:S:R:G:O:kd:Dq:E:vi:o:c:C:N:uyf:xrl");

  if (cl.unrecognised_options_encountered ())
    usage (1);

  verbose = cl.option_count ('v');

  if (verbose)
    substructure_demerits::set_verbose (verbose);

#ifdef USE_IWMALLOC
  if (cl.option_present ('d'))
  {
    int i;
    if (! cl.value ('d', i))
    {
      cerr << "Must specify whole number for malloc death '" << cl.option_value ('d') << "'\n";
      return 1;
    }
    malloc_set_die_when_issued (i);
    if (verbose)
      cerr << "Malloc will die when issuing block " << i << endl;
  }

  if (cl.option_present ('D'))
    iwmalloc_set_full_debug (1);
#endif

  if (cl.option_present ('C'))
  {
    IWString fname = cl.string_value ('C');

    if (! substructure_demerits::initialise_hard_coded_queries_to_do (fname))
    {
      cerr << "Cannot initialise hard coded queries subset (-C option) from '" << fname << "'\n";
      return 8;
    }
  }

  if (cl.option_present('r'))
  {
    substructure_demerits::set_only_apply_rejection_rules ();
    if (verbose)
      cerr << "Only rejection rules will be applied\n";
  }

  if (cl.option_present ('N'))
  {
    Charge_Assigner & charge_assigner = substructure_demerits::charge_assigner ();

    if (! charge_assigner.construct_from_command_line (cl, verbose > 0, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      return 3;
    }
  }

  if (cl.option_present ('f'))
  {
    int f;
    if (! cl.value ('f', f) || f < 1)
    {
      cerr << "The rejection threshold (-f) must be a whole +ve number\n";
      usage(4);
    }

    set_rejection_threshold (f);

    if (verbose)
      cerr << "Rejection threshold set to " << f << endl;

    if (! cl.option_present('x'))
    {
      lower_atom_count_demerit = f;
      upper_atom_count_demerit = f;
    }
  }

  int input_type = 0;
  if (! cl.option_present ('i'))
  {
    if (! all_files_recognised_by_suffix (cl))
    {
      cerr << "Cannot discern input type(s)\n";
      return 3;
    }
  }
  else if (! process_input_type (cl, input_type))
  {
    cerr << prog_name << ": cannot discern input type\n";
    usage (2);
  }
  
  if (cl.option_present ('R'))
  {
    const_IWSubstring fname;
    cl.value ('R', fname);

    if (! cl.option_present ('o'))
      stream_for_rejected_molecules.add_output_type (SMI);
    else if (! stream_for_rejected_molecules.determine_output_types (cl))
    {
      cerr << "Cannot discern file types for -R output\n";
      return 4;
    }

    if (! stream_for_rejected_molecules.new_stem (fname))
    {
      cerr << "Cannot set -R output name to '" << fname << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Rejected structures written to '" << fname << "'\n";
  }

  if (cl.option_present ('G'))
  {
    const_IWSubstring fname;
    cl.value ('G', fname);

    if (! cl.option_present ('o'))
      stream_for_non_rejected_molecules.add_output_type (SMI);
    else if (! stream_for_non_rejected_molecules.determine_output_types (cl))
    {
      cerr << "Cannot discern file types for -G output\n";
      return 5;
    }

    if (! stream_for_non_rejected_molecules.new_stem (fname))
    {
      cerr << "Cannot set -G output name to '" << fname << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Non rejected structures written to '" << fname << "'\n";

//  cerr << "stream_for_non_rejected_molecules number types " << file_for_non_rejected_molecules.number_elements () << endl;
  }

  if (cl.option_present ('S'))
  {
    cl.value ('S', stem_for_demerit_file);

    if (verbose)
      cerr << "Demerit file will have name stem '" << stem_for_demerit_file << "'\n";
  }

  (void) process_elements (cl);

  if (cl.option_present ('M'))
  {
    const_IWSubstring fname;
    cl.value ('M', fname);

    if (! stream_for_multiple_demerits.determine_output_types (cl))
    {
      cerr << "Cannot discern output types for -M output\n";
      return 76;
    }

    if (! stream_for_multiple_demerits.new_stem (fname))
    {
      cerr << "Cannot set -M output stem to '" << fname << "'\n";
      return 18;
    }

    if (verbose)
      cerr << "Molecules rejected for multiple purposes written to '" << fname << "'\n";
  }

  if (cl.option_present ('X'))
  {
    if (! elements_to_remove.construct_from_command_line (cl, verbose, 'X'))
    {
      cerr << "Cannot discern elements to be removed\n";
      usage (21);
    }
  }

  if (cl.option_present ('t'))
  {
    append_demerit_text_to_name = 1;
    if (verbose)
      cerr << "Demerit types will be appended to molecule names\n";
  }

  if (cl.option_present ('u'))
  {
    append_demerit_text_to_name = 1;
    write_rejection_reason_like_tsubstructure = 1;
    set_store_demerit_reasons_like_tsubstructure (1);
    if (verbose)
      cerr << "Rejection reasons written like tsubstructure\n";
  }

  if (cl.option_present ('y'))
  {
    set_demerit_reason_contains_individual_demerits (1);
    write_rejection_reasons_separately = 1;
    if (verbose)
      cerr << "Demerit reasons separated\n";
  }

  if (cl.option_present ('c'))
  {
    int i = 0;
    const_IWSubstring c;
    while (cl.value ('c', c, i++))
    {
      const_IWSubstring directive;
      int dvalue;
      if (! c.split_into_directive_and_value (directive, '=', dvalue))
      {
        cerr << "Invalid -c directive '" << c << "'\n";
        usage (5);
      }

      if (dvalue < 0)
      {
        cerr << "INvalid numeric directive, cannot be negative '" << c << "'\n";
        usage (6);
      }

      if ("hmin" == directive)
      {
        hard_lower_atom_count_cutoff = dvalue;
      }
      else if ("smin" == directive)
      {
        soft_lower_atom_count_cutoff = dvalue;
      }
      else if ("mindmrt" == directive)
      {
      }
      else if ("smax" == directive)
      {
        soft_upper_atom_count_cutoff = dvalue;
      }
      else if ("hmax" == directive)
      {
        hard_upper_atom_count_cutoff = dvalue;
      }
      else if ("maxdmrt" == directive)
      {
      }
      else
      {
        cerr << "Unrecognised -c qualifier '" << c << "'\n";
        usage (5);
      }
    }

    if (0 == soft_lower_atom_count_cutoff && 0 == hard_lower_atom_count_cutoff)
      ;
    else if (hard_lower_atom_count_cutoff <= soft_lower_atom_count_cutoff)
      ;
    else
    {
      cerr << "Invalid lower atom count cutoff values soft:" << soft_lower_atom_count_cutoff << " hard:" << hard_lower_atom_count_cutoff << endl;
      return 5;
    }

    if (0 == soft_upper_atom_count_cutoff && 0 == hard_upper_atom_count_cutoff)
      ;
    else if (soft_upper_atom_count_cutoff < hard_upper_atom_count_cutoff)
      ;
    else
    {
      cerr << "Invalid upper atom count cutoff values soft:" << soft_upper_atom_count_cutoff << " hard:" << hard_upper_atom_count_cutoff << endl;
      return 5;
    }

    if (0 == soft_lower_atom_count_cutoff && 0 == soft_upper_atom_count_cutoff)
      ;
    else if (soft_lower_atom_count_cutoff < soft_upper_atom_count_cutoff)
      ;
    else
    {
      cerr << "Soft cutoffs invalid, lower:" << soft_lower_atom_count_cutoff << " upper:" << soft_upper_atom_count_cutoff << endl;
      return 8;
    }

//  Should do more checks here...
  }

  if (! process_standard_smiles_options (cl, verbose))
  {
    usage (7);
  }

  if (! process_standard_aromaticity_options (cl, verbose))
  {
    usage (8);
  }

  if (cl.option_present ('k'))
  {
    keep_going_after_rejection = 1;
    if (verbose)
      cerr << "Will scan all criteria, ignoring previous rejections\n";

    substructure_demerits::set_keep_going_after_rejection (1);
  }

  if (cl.option_present ('O'))
  {
    IWString tmp;
    int i = 0;
    while (cl.value ('O', tmp, i++))
    {
      if ("hard" == tmp)
      {
        do_hard_coded_substructure_queries = 0;
        if (verbose)
          cerr << "Hard coded substructure queries will be skipped\n";
      }
      else
      {
        cerr << "Unrecognised skip specifier '" << tmp << "'\n";
        usage (19);
      }
    }
  }

  if (cl.option_present ('V'))
  {
    skip_molecules_with_abnormal_valences = 1;
    if (verbose)
      cerr << "Molecules containing unusual valences will be rejected\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << prog_name << ": no files specified\n";
    usage (2);
  }

// As we read in queries, we make sure that each component has a numeric
// value - saves checking them later...

  if (cl.option_present ('r'))
    ;
  else if (cl.option_present ('q'))
  {
    if (! process_queries (cl, queries, verbose, 'q'))
    {
      cerr << prog_name << ": processing of -q option failed\n";
      return 8;
    }

    for (int i = 0; i < queries.number_elements (); i++)
    {
      Substructure_Hit_Statistics * q = queries[i];

      q->set_find_unique_embeddings_only (1);

      for (int j = 0; j < q->number_elements (); j++)
      {
        const Single_Substructure_Query * sq = q->item (j);

        double d;
        if (! sq->numeric_value (d))
        {
          cerr << "Yipes, query '" << q->comment () << " has no demerit value\n";
          return 13;
        }
  
        if (d <= 0.0)
        {
          cerr << "Hmmm, query '" << q->comment () << "' has a non-positive demerit " << d << endl;
          return 14;
        }
      }
    }
  }

  set_remove_hits_not_in_largest_fragment_behaviour(1);    // to get reproducible behaviour with multiple instances of largest fragment

  int rc = 0;

  if (cl.option_present('l'))
  {
    resizable_array_p<Substructure_Hit_Statistics> q1, q2;
    separate_depending_on_fragment_match(queries, q1, q2);

    for (int i = 0; i < cl.number_elements(); i++)
    {
      const char *fname = cl[i];

      if (verbose)
        cerr << prog_name << " processing '" << fname << "'\n";

      if (! iwdemerit (fname, input_type, q1, q2))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    resizable_array_p<Substructure_Hit_Statistics> notused;

    for (int i = 0; i < cl.number_elements (); i++)     // each argument is a file
    {
      const char *fname = cl[i];

      if (verbose)
        cerr << prog_name << " processing '" << fname << "'\n";

      if (! iwdemerit (fname, input_type, queries, notused))
      {
        rc = i + 1;
        break;
      }
    }
  }

  if (verbose)
  {
    if (cl.number_elements () > 1)
      cerr << molecules_read << " molecules read\n";

    if (molecules_with_abnormal_valences)
      cerr << "Skipped " << molecules_with_abnormal_valences << " molecules with abnormal valences\n";

    for (int i = 0; i < demerits_per_molecule.number_elements (); i++)
    {
      int j = demerits_per_molecule[i];
      if (j)
        cerr << j << " molecules had " << i << " sources of demerits\n";
    }

    for (int i = 0; i < demerits_per_rejected_molecule.number_elements (); i++)
    {
      int j = demerits_per_rejected_molecule[i];
      if (j)
        cerr << j << " rejected molecules had " << i << " sources of demerits\n";
    }
  }

#ifdef USE_IWMALLOC
  terse_malloc_status (stderr);
#endif

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

#ifdef USE_IWMALLOC
  cerr << "................. Running with debugging malloc ......................\n";
//iwmalloc_set_full_debug (1);
#endif
  
  int rc = iwdemerit (argc, argv);

#ifdef USE_IWMALLOC
  malloc_status (stderr);
#endif

  return rc;
}
