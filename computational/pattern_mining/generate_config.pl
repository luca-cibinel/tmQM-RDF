#!/usr/bin/perl -w
# vi:si:ts=4:sw=4
use strict;

use YAML::Syck qw(Dump LoadFile DumpFile);
use Clone qw(clone);
use MCE::Loop;
use IO::File;
use IO::Uncompress::Gunzip (); # for reading from compressed file
use IO::Compress::Gzip qw(gzip $GzipError); # for writing to compressed file

use Data::Dumper;

my $cfg_dir = "./configs/";

my $session_keyword = "test"; # a keyword used to identify both the graph dataset in intermediate/local_data/ and the results directory in ./results

my $CFG = {
            version => "111125",
            comments => [
                "Test for deterministic sampling (part 2)"
            ],
            
            input_files                 => [],
            input_directory                => "./intermediate/local_data/" . $session_keyword . "/",
            verbose                 => 1,
            incremental_output_file            => "./results/" . $session_keyword . "/output-size-%d",
            incremental_output_file_per_worker    => "./results/" . $session_keyword . "/output-size-%d-worker-%d",
            incremental_input_file_per_worker    => "./results/" . $session_keyword . "/input-size-%d-worker-%d",
            triplepattern_output_file        => "./results/" . $session_keyword . "/triplepatterns.yml",

            output_readable_patterns        => 0,
            output_parseable_patterns        => 1,
            output_list_of_matched_graphs        => 0,
            output_matched_graphs_bitstring        => 0,
            output_number_of_matched_graphs     => 0,
            output_mappings                => 1, # better set to 0, otherwise result files might be very large
            output_number_of_mappings        => 0,

            min_pattern_size            => 1, # this parameter is not necessary anymore, as each round outputs resuls to file
            max_pattern_size            => 6,
            min_matches_per_graph             => 1,
            max_matches_per_graph             => undef,
            min_support_per_graph            => 1,     # TODO: not used yet
            max_support_per_graph            => undef, # TODO: not used yet
            min_number_of_matched_graphs_abs     => 4,
            max_number_of_matched_graphs_abs     => undef,
            min_number_of_matched_graphs_rel     => undef,
            max_number_of_matched_graphs_rel     => undef,

            max_share_of_varnodes            => 1, # value between 0 and 1
            start_patternsize_for_share_of_varnodes_constraint => 4,

            #extension_sample_threshold        => 0, # default is 0. value needs to be >=0 and < 1. the higher the value, the more extensions will be discarded. extensions are discareded randomly (uniformely). TODO sanity checks.

            deterministic_sampling =>  1,
            random_number_file_1 => "./intermediate/random_numbers/" . $session_keyword . "/random_numbers_1.txt",
            random_number_file_2 => "./intermediate/random_numbers/" . $session_keyword . "/random_numbers_2.txt",

            # TODO: might be cool to ask the user during runtime to define the threshold for the next round.
            extension_sample_threshold => {
                #"any" => 0.77,        # "any" means, does not depend on the pattern size. if any is used, then the other thresholds defined here are ignored
                2 => 0.0,
                3 => 0.55,
                4 => 0.55,
                5 => 0.55,
                6 => 0.55,
                7 => 0.66,
                8 => 0.82,
                9 => 0.88,
                10 => 0.9,
                11 => 0.9,
                12 => 0.9
            },

            number_of_workers => 15,

            max_matches_per_graph_per_patternsize => { # TODO sanity checks, compare to min/max_matches_per_graph
                #2 => 10,
                #3 => 15,
                #4 => 20,
            },
    # <resource://integreat/p5/molecule/ABATEC> <resource://integreat/p5/chemical_information/feature_charge> "0"^^<http://www.w3.org/2001/XMLSchema#float> .
            seed_triple_patterns => [
                [
                    '?v0',
                    '<resource://integreat/p5/complex/TMC/hasMetalCentre>',
                    '?v1'
                ],
            ],

            xxx_undesired_subgraphs => [
                [
                    [
                        "?v1",
                        "<http://example.org/sentence#next>",
                        "?v2"
                    ],
                    [
                        "?v2",
                        "<http://example.org/sentence#pos>",
                        "\"PUNCT\""
                    ],
                ],
            ],

            allowed_abstractions => {
                s     => 1,
                p     => 0,
                o     => 1,
                sp    => 0,
                so    => 1,
                po     => 0,
                spo     => 0, # TODO: not supported yet
                spo_xxx => 0, # TODO: not supported yet
                spo_xxy => 0, # TODO: not supported yet
                spo_xyx => 0, # TODO: not supported yet
                spo_xyy => 0, # TODO: not supported yet
                spo_xyz => 0, # TODO: not supported yet
            },
            allowed_join_types => {
                "s-s" => 1,
                "s-p" => 0,
                "s-o" => 1,
                "p-s" => 0,
                "p-p" => 1,
                "p-o" => 0,
                "o-s" => 1,
                "o-p" => 0,
                "o-o" => 1,
            },
            predicate_blacklist => {
                "http://www.w3.org/1999/02/22-rdf-syntax-ns#type" => 1, # Useless information
                # ATOMIC LEVEL
                # "resource://integreat/p5/atomic/structure/b" => 1, # We don't focus on atomic bonds
                #LIGAND LEVEL
                "resource://integreat/p5/ligand/centre/hasAtom" => 1, # We only retrieve centre information at ligand level
                "resource://integreat/p5/ligand/ligand/hasAtom" => 1 # We only focus on binding atoms
            },
            domain_specific_constraints => {
                #argumentation_abs_1 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_2 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_3 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_4 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_5 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_6 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_7 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_8 => 1, # about abstractions of triples to triple patterns
                #argumentation_abs_9 => 1, # about abstractions of triples to triple patterns
                #argumentation_shape_1 => 1, # about building patterns with multiple "next" edges
            },

            triple_abstraction_constraints => {
                p_based => { # TODO <> are inconvenient
                    #"<resource://integreat/p5/ligand/centre/isMetalCentre>" => { s => 1, p => 0, o => 0, sp => 0, po => 0, so => 0 },
                    "<resource://integreat/p5/ligand/ligand/isLigand>" => { s => 1, p => 0, o => 0, sp => 0, po => 0, so => 0 },
                    "<resource://integreat/p5/atomic/atom/isAtom>" => { s => 1, p => 0, o => 0, sp => 0, po => 0, so => 0 }
                },
            },

            prefix_definitions => {
                "resource://integreat/p5/"    => "inp5",
                "resource://integreat/p5/datasets/"    => "inp5ds",
                "resource://integreat/p5/complex/"    => "cm",
                "resource://integreat/p5/complex/TMC/"    => "cmT",
                "resource://integreat/p5/complex/TMC/property/"    => "cmTp",
                "resource://integreat/p5/ligand/"    => "lg",
                "resource://integreat/p5/ligand/centre/"    => "lgC",
                "resource://integreat/p5/ligand/centre/property/"    => "lgCp",
                "resource://integreat/p5/ligand/centre/reference/"    => "lgCr",
                "resource://integreat/p5/ligand/centre/reference/property/"    => "lgCrp",
                "resource://integreat/p5/ligand/ligand/"    => "lgL",
                "resource://integreat/p5/ligand/ligand/property/"    => "lgLp",
                "resource://integreat/p5/ligand/ligand/reference/"    => "lgLr",
                "resource://integreat/p5/ligand/ligand/reference/property/"    => "lgLrp",
                "resource://integreat/p5/ligand/ligand/reference/structure/"    => "lgLrs",
                "resource://integreat/p5/ligand/ligand/reference/reference/"    => "lgLrr",
                "resource://integreat/p5/ligand/bond/"    => "lgB",
                "resource://integreat/p5/ligand/bond/property/"    => "lgBp",
                "resource://integreat/p5/ligand/bond/reference/"    => "lgBr",
                "resource://integreat/p5/ligand/bond/reference/property/"    => "lgBrp",
                "resource://integreat/p5/ligand/structure/"    => "lgS",
                "resource://integreat/p5/atomic/"    => "tm",
                "resource://integreat/p5/atomic/atom/"    => "tmA",
                "resource://integreat/p5/atomic/atom/property/"    => "tmAp",
                "resource://integreat/p5/atomic/atom/reference/"    => "tmAr",
                "resource://integreat/p5/atomic/atom/reference/property/"    => "tmArp",
                "resource://integreat/p5/atomic/bond/"    => "tmB",
                "resource://integreat/p5/atomic/bond/property/"    => "tmBp",
                "resource://integreat/p5/atomic/bond/reference/"    => "tmBr",
                "resource://integreat/p5/atomic/bond/reference/property/"    => "tmBrp",
                "resource://integreat/p5/atomic/structure/"    => "tmS",
                "http://www.w3.org/2001/XMLSchema#"    => "xmls"
            },

            pattern_visualization => {
                visualize_patterns => 1,
                visualization_directory => "./results/" . $session_keyword . "/patterns-%d",
                filetype => "png",
                layout => "dot",
                # TODO have, not in here, an option visualize_only. skip mining, do visualization of what was mined in a previous run.
                # maybe, define a string here that is used in the label of the graph
            },
        };

if(exists $CFG->{input_directory}){
    unless(-d $CFG->{input_directory}){
        die "Local data not found! Please populate the directory " . $CFG->{input_directory} . " with your local data.";
    }
    
    my $cnt = 0; # TODO remove
    foreach my $file (sort glob($CFG->{input_directory} . "*.nt")){
        push(@{$CFG->{input_files}}, $file);
    }
}else{
    die "Local data folder not specified! Specify an input_directory entry."
}

unless(-d ("./results/" . $session_keyword)){
    print(("./results/" . $session_keyword) . " does not exist! Creating it now...");
    mkdir ("./results/");
    mkdir ("./results/" . $session_keyword . "/");
}

unless(-d $cfg_dir){
    print($cfg_dir . " does not exist! Creating it now...");
    mkdir $cfg_dir;
}
DumpFile($cfg_dir . $session_keyword . ".cfg", $CFG);

