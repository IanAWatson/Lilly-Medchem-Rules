Lilly-Medchem-Rules
===================

This is an implementation of Eli Lilly Medchem Rules.  They were
published under "Rules for Identifying Potentially Reactive or
Promiscuous Compounds" by Robert F. Bruns and Ian W. Watson,
J. Med. Chem. 2012, 55, 9763--9772 as ACS Author choice, i.e. open
access at [doi 10.1021/jm301008n](https://doi.org/10.1021/jm301008n).


To quote the abstract, "[This approach] describes a set of 275 rules,
developed over an 18-year period, used to identify compounds that may
interfere with biological assays, allowing their removal from
screening sets. Reasons for rejection include reactivity (e.g., acyl
halides), interference with assay measurements (fluorescence,
absorbance, quenching), activities that damage proteins (oxidizers,
detergents), instability (e.g., latent aldehydes), and lack of
druggability (e.g., compounds lacking both oxygen and nitrogen). The
structural queries were profiled for frequency of occurrence in
druglike and nondruglike compound sets and were extensively reviewed
by a panel of experienced medicinal chemists. As a means of profiling
the rules and as a filter in its own right, an index of biological
promiscuity was developed. The 584 gene targets with screening data at
Lilly were assigned to 17 subfamilies, and the number of subfamilies
at which a compound was active was used as a promiscuity index."

Scrutinizing the SMILES string of a molecule, the program identifies
pattern which are are a knock-off for a candidate, such as unwanted
elements (e.g., Ag, Hg, Zn), or a too low atom count (less than 7 heavy
atoms).  More importantly, however, most of the rules are _scaled_
in respect of each other.  Thus, the _demerits_ of a butyl, pentyl,
hexyl, heptyl, and cyclohexyl group in a molecule equate 10, 25, 50,
100, and 170, respectively.  Eventually, the demerits of a molecule
are summed up and compared with an arbitrary threshold; the program's
adjustable default cut-off equates to 100.  (For details adjusting the
settings, see the program's documentation.)

After downloading the software, multiple options are offered to
install the program.  These include a compilation with `make` e.g.,
in Cygwin or Linux Ubuntu, or as docker file and are documented in
dedicated `.md` files.  In addition to C++, an installation of Ruby
is required allowing to perform a basic functionality by

`ruby Lilly_Medchem_Rules.rb input.smi > okmedchem.smi`

where `input.smi` is your collection of SMILES strings of structures
to be checked.  In file `okmedchem.smi` the program lists the molecules
whose added demerits are below the critical threshold applied.  If
applicable, molecules with demerits equal or greater than the threshold
will be reported in one or multiple additional files written by the
program, e.g., `bad0.smi`.

The freely accessible publication and its supplementary material at ACS
outline the structural pattern scrutinized and their demerits.  As a
user training, this repository contains multiple reference `.smi` files
in folder `test` to probe the program and its options.  Simulating a
screening, applying the default parameters, the collection of 24986
pubchem molecules in file `example_molecules.smi` yields a set of 4576
acceptable molecules reported in file `okmedchem.correct.smi`.

Test data `table_S3.smi`, retrieved from [table S3](https://pubs.acs.org/doi/suppl/10.1021/jm301008n/suppl_file/jm301008n_si_001.pdf)
of publication, and `200_prescriptions_2011.smi`, retrieved from a
cross-linked [Wikipedia project](https://en.wikipedia.org/wiki/Wikipedia:WikiProject_Pharmacology/Top_200_US_Prescriptions_2011)
are provided to illustrate the outcome of this set of rules among drugs
eventually marketed.  These of course represent a stage of development
much later than the of screening the program targets.

## Adjusting Defaults

By default, molecules with fewer than 7 heavy atoms are rejected. Molecules with
between 25 and 40 heavy atoms are progressively demerited, and after 40, molecules are rejected.
These upper limits are referred to as the soft and hard upper cutoffs. 
Atom count paramters can be adjusted via the `-c` (lower cutoff) and
`-Cs` and `-Ch` (soft and hard upper atom count) options.  For example
if you wanted to filter to molecules that contained between 10 and 40
atoms with no demerits, and then reject at 50 heavy atoms, that could
be

`ruby Lilly_Medchem_Rules.rb -c 10 -Cs 40 -Ch 50 input.smi > okmedchem.smi`

Unfortunately as currently implemented the combination `-Cs 50 -Ch 50` does
not work, so if you want to largely avoid demerits for heavy atom count try
something like `-Cs 50 -Ch 51`.
