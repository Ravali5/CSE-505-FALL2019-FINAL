
:- use_module(library(lib)).
:- lib( bio_db ).

/** Use Case 1 : bio_db_stats gives list of predicates installed  */

bio_db_stats :-
    bio_db_stats( _Trips ),
    halt.

bio_db_stats_csv :-
    bio_db_stats( Trips ),
    findall( row(Pn,Pa,Len), member(Pn/Pa-Len,Trips), Rows ),
    lib( mtx ),
    mtx( 'bio_db_stats.csv', Rows ),
    halt.

bio_db_stats( Trips ) :-
    findall( Pn/Pa-Len, bio_db_pred_stats(Pn,Pa,Len), Trips ),
    length( Trips, NoPreds ),
    findall( Len, member(_-Len,Trips), Lens ),
    sumlist( Lens, NoRecs ).

bio_db_pred_stats( Pn, Pa, Len ) :-
    current_predicate( bio_db:Pn/Pa ),
    member( Pfx, [map_,edge_] ),
    atom_concat( Pfx, _, Pn ),
    functor( G, Pn, Pa ),
    \+ atom_concat( _, info, Pn ),
    findall( 1, G, Ones ),
    length( Ones, Len ),
    debug( bio_db_stats, '~w/~d has ~d records.', [Pn,Pa,Len] ).


%% Use Case 2 : For a specified gene symbol this function finds the corresponding Gene Ontology term

lmtk3_go(Symb1, Oymbs) :-
	map_gont_symb_gont( Symb1, Gont ),
	findall( Symb, map_gont_gont_symb(Gont,Symb), Symbs ),
	%map_gont_gont_gonm( Gont, Gonm ),
	sort( Symbs, Oymbs ).
	%length( Oymbs, Len )


gene_pops( GeneIn ) :-
	upcase_atom( GeneIn, Gene ),
	findall( row(GoT,GoN,GoP), gene_pops(Gene,GoT,GoN,GoP), Rows ),
	atomic_list_concat( [pop,Gene], '_', Stem ),
	file_name_extension( Stem, csv, MtxF ),
	mtx( MtxF, Rows ).

gene_pops( Gene, Gont, Gonm, Len ) :-
	map_gont_symb_gont( Gene, Gont ),
	findall( Symb, map_gont_gont_symb(Gont,Symb), Symbs ),
	map_gont_gont_gonm( Gont, Gonm ),
	sort( Symbs, Oymbs ),
	length( Oymbs, Len ),
	Mess = '~w gene ontology term: ~a-~a, with population of: ~d',
	debug( lmtk3, Mess, [Gene,Gont,Gonm,Len] ).

	go_plot( Gont ) :-
	findall( Symb, map_gont_gont_symb(Gont,Symb), Symb1 ),
	edge_strg_hs_symb( Symb1, Symb2, 941 ),
	Popts = [plotter(igraph),vertex.size=4,orphan_edge_weight(0.001),
	         stem(Gont)],
	wgraph_plot( Symb2, Popts ).


%% Use Case 3 : Finding the edges between the gene symbols for a specific Gene ontology term.

go_term_graph(Gont,Lim,Graph):-
	findall( Symb, map_gont_gont_symb(Gont,Symb), Symbs ),
	findall( Symb1-Symb2:W, (
	member(Symb1,Symbs),
	member(Symb2,Symbs),
	edge_strg_hs_symb(Symb1,Symb2,W),
	Lim < W
	),
	Graph ).

%% Use Case 4 : To find how many ensemble genes are part of hgnc and not part of ense .

gene_map_compare :-
    findall( EnsG, ( map_hgnc_ensg_hgnc(EnsG,Hgnc),
                  \+ map_ense_ensg_hgnc(EnsG,Hgnc)
                            ), NotInEnse ),
    findall( EnsG, ( map_ense_ensg_hgnc(EnsG,Hgnc),
                     \+map_hgnc_ensg_hgnc(EnsG,Hgnc) ),
                        NotInHgnc ),
    findall( EnsG, map_hgnc_ensg_hgnc(EnsG,_HgncH), EnsGsH ),
    sort( EnsGsH, EnsGsHo ),
    length( EnsGsHo, NofGsH ),
    findall( EnsG, map_ense_ensg_hgnc(EnsG,_HgncE), EnsGsE ),
    sort( EnsGsE, EnsGsEo ),
    length( EnsGsEo, NofGsE ),
    length( NotInEnse, NofNIE ),
    length( NotInHgnc, NofNIH ),
    write( not_in_ensemble(NofNIE/NofGsH) ), nl,
    write( not_in_hgnc(NofNIH/NofGsE) ), nl,
    findall( EnsG-(HgncH,HgncE), (
                                    map_hgnc_ensg_hgnc(EnsG,HgncH),
                                    map_ense_ensg_hgnc(EnsG,HgncE),
                                    HgncH \== HgncE
                                 ),
                                    Trips ),
    length( Trips, NofTrips ),
    write( nof_diffs(NofTrips) ), nl,
    write( diffs(Trips) ), nl.



imported_from( Pname/Parity, Mod ) :-
    !,
    functor( Goal, Pname, Parity ),
    imported_from( Goal, Mod ).
imported_from( Goal, Mod ) :-
    predicate_property( Goal, imported_from(Mod) ),
    !.
imported_from( _Goal, user ).

%% Use Case 5 : Display all the predicates available along with the semantics and from where they are imported .

lib_preds :-
    lib_pred( Pred ),
    write( Pred ),  nl,
    fail.
lib_preds.

lib_pred( Pred-From ) :-
    predicate_property( Pred, imported_from(From) ),
    atom_concat( bio_db, _, From ).


