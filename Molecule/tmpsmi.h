  private:


    int _find_raw_rings (const atom_number_t previous_atom,
                         const atom_number_t current_atom,
                         resizable_array<Ring *> & rings,
                         resizable_array<atom_number_t> & active_rings,
                         int * already_done);
    int _find_raw_rings_for_fragment (int id, int * already_done);
    int _find_raw_rings (int * already_done);
    int _find_raw_rings ();

//  Various functions for building smiles

    int _construct_smiles_for_fragment (const int * zorder,
                                        int * already_done,
                                        IWString & smiles,
                                        Ring_Number_Manager & rnm,
                                        atom_number_t previous_atom,
                                        atom_number_t a);
    int _construct_smiles (const int *, int *, IWString &);
    int _construct_smiles (const int *, IWString &);

    int _build_smiles_ordering (int * zorder,
                                int frag_order,
                                int (*identify_next_atom) (Molecule *, const int *, atom_number_t, atom_number_t &),
                                const atom_number_t previous_atom,
                                const atom_number_t a,
                                int & icounter);
    int _build_smiles_ordering (int * zorder,
                                int (* identify_start_atom) (Molecule *, const int *, atom_number_t & a),
                                int (* identify_next_atom) (Molecule *, const int *, atom_number_t, atom_number_t &));
    int _include_atom_in_smiles (atom_number_t) const;
    int _mark_atoms_not_in_smiles (int * zorder);
    int _build_smiles_ordering ();

  public:

