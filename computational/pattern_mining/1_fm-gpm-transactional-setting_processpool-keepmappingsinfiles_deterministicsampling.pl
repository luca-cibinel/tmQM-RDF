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

my $CFG;
if($ARGV[0]){
    my $config_file = $ARGV[0];
    $CFG = LoadFile($config_file);
    #print "warning: reduced the number of workers. and number of files. and other things\n"; <STDIN>;
    #$CFG->{number_of_workers} = 2; # TODO remove
    #$CFG->{min_number_of_matched_graphs_abs} = 2; # TODO remove
} else {
    print "Using default configuration generated on the fly\n";
    $CFG = {
                version => "111125",
                comments => [
                    "Test for deterministic sampling"
                ],
                
                input_files                 => [],
                input_directory                => "./local_data/",
                verbose                 => 1,
                incremental_output_file            => "./results/output-size-%d",
                incremental_output_file_per_worker    => "./results/output-size-%d-worker-%d",
                incremental_input_file_per_worker    => "./results/input-size-%d-worker-%d",
                triplepattern_output_file        => "./results/triplepatterns.yml",

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
                random_number_file_1 => "random_numbers_1.txt",
                random_number_file_2 => "random_numbers_2.txt",

                # TODO: might be cool to ask the user during runtime to define the threshold for the next round.
                extension_sample_threshold => {
                    #"any" => 0.77,        # "any" means, does not depend on the pattern size. if any is used, then the other thresholds defined here are ignored
                    2 => 0.0, # Large: ?, small: 0.45 # TODO add warning to sanity checks
                    3 => 0.55, # Large: ?, small: 0.55
                    4 => 0.55, # Large: ?, small: 0.55
                    5 => 0.55, # Large: ?, small: 0.55
                    6 => 0.55, # Large: ?, small: 0.55
                    7 => 0.66, # Large: ?, small: 0.66
                    8 => 0.82, # Large: ?, small: 0.82
                    9 => 0.88, # Large: ?, small: 0.9
                    10 => 0.9, # Large: ?, small: 0.925
                    11 => 0.9, # Large: ?, small: 0.94
                    12 => 0.9, # Large: ?, small: 0.96
                    13 => 0.92,
                    14 => 0.92,
                    15 => 0.92,
                    16 => 0.92,
                    17 => 0.95,
                    18 => 0.95,
                    19 => 0.95,
                    20 => 0.95
                    #21 => 0.9,
                    #22 => 0.9,
                    #23 => 0.9,
                    #24 => 0.9,

                },

                number_of_workers => 4,

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
                    visualize_patterns => 0,
                    visualization_directory => "./results/patterns-%d",
                    filetype => "png",
                    layout => "dot",
                    # TODO have, not in here, an option visualize_only. skip mining, do visualization of what was mined in a previous run.
                    # maybe, define a string here that is used in the label of the graph
                },
            };


    print "VERSION = $CFG->{version}", "\n";
    if (scalar @{$CFG->{comments}} > 0){
        print "COMMENTS:", "\n";
        foreach my $comment (@{$CFG->{comments}}){
            print "\t", "$comment", "\n";
        }
    }
    print "extension_sample_threshold:";
    print Dumper($CFG->{extension_sample_threshold});
    
    if(exists $CFG->{input_directory}){
        my $cnt = 0; # TODO remove
        foreach my $file (sort glob($CFG->{input_directory} . "*.nt")){
            push(@{$CFG->{input_files}}, $file);
            #print "used $file\n";
    #        last if ++$cnt > 10; # TODO remove
        }
    }

    DumpFile("default.cfg", $CFG);
    print " > default.cfg\n"; #<STDIN>;
    #$CFG = LoadFile("default.cfg");
    
}

# TODO: decide where this should be done
#if(exists $CFG->{input_directory}){
#    my $cnt = 0; # TODO remove
#    foreach my $file (sort glob($CFG->{input_directory} . "*.nt")){
#        push(@{$CFG->{input_files}}, $file);
#        last if ++$cnt > 5000; # TODO remove
#    }
#}

if(exists $CFG->{min_number_of_matched_graphs_rel} and defined $CFG->{min_number_of_matched_graphs_rel}){
    $CFG->{min_number_of_matched_graphs_abs} = int ($CFG->{min_number_of_matched_graphs_rel} * scalar @{$CFG->{input_files}});
}

if(exists $CFG->{max_number_of_matched_graphs_rel} and defined $CFG->{max_number_of_matched_graphs_rel}){
    $CFG->{max_number_of_matched_graphs_abs} = int ($CFG->{max_number_of_matched_graphs_rel} * scalar @{$CFG->{input_files}});
}

#print Dump { CFG => $CFG }; <STDIN>;

&config_sanity_checks($CFG);

#my $term_to_termID = {};
#my $termID_to_term = {};
#my $last_term_ID = -1;

my $tpstring_to_tpID = {};
my $tpID_to_tp = {};
my $last_pattern_ID = 0;


#print "warning: changing CFG\n"; <STDIN>;
#$CFG->{number_of_workers} = 3;
#$CFG->{max_pattern_size} = 7;


my $continuation_is_possible = 1;
my $continuation_from_size = 2;

if(not -s $CFG->{triplepattern_output_file}){
    $continuation_is_possible = 0;
} else {

        my $pattern_sh_1 = $CFG->{incremental_output_file} . ".empty";
    $pattern_sh_1 =~ s/\%d/*/;

    if(glob($pattern_sh_1)){
        print "Continuation not possible, nothing to do, end of growth has been reached.\n";
        exit;
    }

    # what if a tar was created in the previous run?
    # solution: untar file for largest pattern size
    my $largest_size = 0;

    my $tar_file_pattern_sh = $CFG->{incremental_input_file_per_worker} . ".dat.tar.gz";
    $tar_file_pattern_sh =~ s/\%d/*/;
    $tar_file_pattern_sh =~ s/\%d/all/;

    my $tar_file_pattern_pl = $CFG->{incremental_input_file_per_worker} . ".dat.tar.gz";
        $tar_file_pattern_pl =~ s/\%d/(\\d+)/;
        $tar_file_pattern_pl =~ s/\%d/all/;

    my @tar_files = glob($tar_file_pattern_sh);

    my @numbers0;
    foreach my $tar_file (@tar_files){
        if($tar_file =~ m/$tar_file_pattern_pl/){
            push (@numbers0, $1);
        }
    }

    if(scalar @numbers0){
        @numbers0 = reverse sort { $a <=> $b } @numbers0;
        my $number = $numbers0[0];
        if($number > $CFG->{max_pattern_size}){
            print "Continuation not possible with currently set max_pattern_size.\n";
            exit;
        }
        my $tar_file_name = $CFG->{incremental_input_file_per_worker} . ".dat.tar.gz";
                $tar_file_name =~ s/\%d/$number/;
        $tar_file_name =~ s/\%d/all/;
        print "untar $tar_file_name\n";
        system("tar -xvf $tar_file_name");
    }
    
    #print "wait\n"; <STDIN>;

    my $pattern_sh_2 = $CFG->{incremental_input_file_per_worker} . ".dat.gz.Z";
        my $pattern_pl = $pattern_sh_2;
    
    $pattern_sh_2 =~ s/\%d/*/;
    $pattern_sh_2 =~ s/\%d/1/;
    
    $pattern_pl =~ s/\%d/(.*)/;
        $pattern_pl =~ s/\%d/1/;

    my @numbers = ();
    foreach my $file (glob($pattern_sh_2)){
        if($file =~ m/$pattern_pl/){
            push (@numbers, $1);
        }
    }

    @numbers = reverse sort @numbers;
    #print Dump { numbers => \@numbers }; <STDIN>;

    if(scalar @numbers){
        $continuation_from_size = $numbers[0];
        print "continuation_from_size: $continuation_from_size\n"; <STDIN>;
    } else {
        $continuation_is_possible = 0;
    }


}

#print "wait (size: $continuation_is_possible, $continuation_from_size)\n"; <STDIN>; # thus, try to create patterns of this size by using patterns of size -1.

if($continuation_is_possible){
    print "Continuation is possible. Continuing to grow patterns of size $continuation_from_size\n";
}

my $can_skip_initial_phase = $continuation_is_possible && $continuation_from_size > 2 ? 1 : 0;
print "can_skip_initial_phase: $can_skip_initial_phase\n"; #<STDIN>;

if(not $can_skip_initial_phase){
    # initial phase - part 1/6: create list of all admissible triple patterns
    my $T = {};
    my $cnt = 1;
    foreach my $input_file (@{$CFG->{input_files}}){
        print &get_timestamp() . " read $cnt files (phase 1)\n" if $cnt % 1000 == 0;
        $cnt++;
        #last if $cnt > 5000; # TODO remove
        open(DAT,"<$input_file") or die "$! <$input_file>"; # TODO remove die
        WHL: while(defined(my $line=<DAT>)){ # WHL
            next if $line =~ m/\A# / or $line =~ m/\A\n\Z/;
            my $obj = &parse_NT_into_obj($line);
            next if not defined $obj;
            $obj->{p}->{rep} =~ m/\A<(.*)>\Z/;
            my $uri_value = $1;
            if(defined $CFG->{predicate_blacklist} and scalar keys %{$CFG->{predicate_blacklist}}){
                next if exists $CFG->{predicate_blacklist}->{$uri_value} and $CFG->{predicate_blacklist}->{$uri_value};
            }
            if(defined $CFG->{predicate_whitelist} and scalar keys %{$CFG->{predicate_whitelist}}){
                next if not exists $CFG->{predicate_whitelist}->{$uri_value};
            }

            my ($s, $p, $o) = ($obj->{s}->{rep}, $obj->{p}->{rep}, $obj->{o}->{rep});

            foreach my $tpID_mu (&AdmissibleAbstractions($s,$p,$o)){
                my ($tpID, $mu) = @{$tpID_mu};
                
                $T->{$tpID}->{$input_file} = 1;
            }

        } # WHL
        close DAT;
    }
    my $cnt_triple_patterns = scalar keys %{$T};
    my $log = int (log($cnt_triple_patterns)/log(10));

    # initial phase - part 2/6: reduce list of triple patterns to those that are frequent and add tpID to each tp
    foreach my $tpID (keys %{$tpID_to_tp}){
        if(scalar keys %{$T->{$tpID}} < $CFG->{min_number_of_matched_graphs_abs}){
            delete $tpID_to_tp->{$tpID};
        }
    }
    $T = {}; # i.e., clean up

    # initial phase - part 3/6: find tpID for seed triple patterns
    if(defined $CFG->{seed_triple_patterns} and scalar @{$CFG->{seed_triple_patterns}}){
        foreach my $seed_triple_pattern (@{$CFG->{seed_triple_patterns}}){

            #if(
            #    $seed_triple_pattern->[0] !~ m/\A\?/ or
            #    $seed_triple_pattern->[1] !~ m/\A\?/ or
            #    $seed_triple_pattern->[2] !~ m/\A\?/
            #){
            #    print "Warning: a seed triple pattern is specified that contains a term that has not been seen in the data:\n";
            #    print Dump { seed_triple_pattern => $seed_triple_pattern };
            #    print "Hit Return to ignore and continue"; <STDIN>;
            #    next;
            #}

            my $can = &canonicalize_triple_pattern($seed_triple_pattern);
            my $tpstring = join("-", $can->[0], $can->[1], $can->[2]);
            if(exists $tpstring_to_tpID->{$tpstring}){
                my $tpID = $tpstring_to_tpID->{$tpstring};
                $CFG->{seed_triple_pattern_IDs}->{$tpID} = 1;
            } else{
                print "Warning: a seed triple pattern is specified that is infrequent:\n";
                print Dump { seed_triple_pattern => $seed_triple_pattern };
                print "Hit Return to ignore and continue"; <STDIN>;
            }
        }
    }

    # initial phase - part 4/6: find tpIDs for tps in undesired subgraphs
    if(defined $CFG->{undesired_subgraphs} and scalar @{$CFG->{undesired_subgraphs}}){ print "taking care of undesired subgraphs\n";
        my $undesired_subgraph_ID = -1;
        FEUS: foreach my $undesired_subgraph (@{$CFG->{undesired_subgraphs}}){
            $undesired_subgraph_ID++;
            foreach my $tp (@{$undesired_subgraph}){
                my $can = &canonicalize_triple_pattern($tp);
                my $tpstring = join("-", $can->[0], $can->[1], $can->[2]);
                if(exists $tpstring_to_tpID->{$tpstring}){
                    my $tpID = $tpstring_to_tpID->{$tpstring};
                    $tp->[3] = $tpID;
                    $CFG->{undesired_subgraphs_tpcount}->{$undesired_subgraph_ID}->{$tpID}++;
                } else {
                    print "Warning: an undesired subgraph is specified that contains a triple pattern that is infrequent:\n";
                    print Dump { undesired_subgraph => $undesired_subgraph, tp => $tp, can => $can };
                    print "Hit Return to ignore and continue"; <STDIN>;
                    next FEUS;
                }
            }
        }
    }

    # initial phase - part 5/6: create Omega for each frequent triple pattern
    my $number_of_mappings = {};
    $cnt = 1;
    foreach my $input_file (@{$CFG->{input_files}}){
        print &get_timestamp() . " read $cnt files (phase 2)\n" if $cnt % 1000 == 0;
        $cnt++;
        #last if $cnt > 5000; #TODO remove
        open(DAT,"<$input_file");
        WHL: while(defined(my $line=<DAT>)){ # WHL
            next if $line =~ m/\A# / or $line =~ m/\A\n\Z/;
            my $obj = &parse_NT_into_obj($line);
            next if not defined $obj;
            $obj->{p}->{rep} =~ m/\A<(.*)>\Z/;
            my $uri_value = $1;
            if(defined $CFG->{predicate_blacklist} and scalar keys %{$CFG->{predicate_blacklist}}){
                next if exists $CFG->{predicate_blacklist}->{$uri_value} and $CFG->{predicate_blacklist}->{$uri_value};
            }
            if(defined $CFG->{predicate_whitelist} and scalar keys %{$CFG->{predicate_whitelist}}){
                next if not exists $CFG->{predicate_whitelist}->{$uri_value};
            }

            my ($s, $p, $o) = ($obj->{s}->{rep}, $obj->{p}->{rep}, $obj->{o}->{rep});

            foreach my $tpID_mu (&AdmissibleAbstractions($s,$p,$o)){
                my ($tpID, $mu) = @{$tpID_mu};
                if(exists $tpID_to_tp->{$tpID}){
                    push(@{$tpID_to_tp->{$tpID}->{omega}->{$input_file}}, $mu);
                    $number_of_mappings->{$tpID}++;
                }
            }

        } # WHL
        close DAT;
    }

    # initial phase - part 6/6: store tpID_to_tp, so that recontinuation of the growth becomes possible
    # TODO YAML might not be the best format
    DumpFile($CFG->{triplepattern_output_file}, $tpID_to_tp);

    # growth phase - part 1/2: preparations

    my $TP = {};
    foreach my $tpID (sort keys %{$tpID_to_tp}){
        my $tp = $tpID_to_tp->{$tpID};

        if(not exists $CFG->{seed_triple_pattern_IDs} or exists $CFG->{seed_triple_pattern_IDs}->{$tpID}){

            my $list_of_matched_graphs = [];
            if($CFG->{output_list_of_matched_graphs}){
                foreach my $graphID (sort { $a <=> $b } keys %{$tp->{omega}}){
                    push (@{$list_of_matched_graphs}, $graphID);
                }
            }
            my $matched_graphs_bitstring = q{};
            if($CFG->{output_matched_graphs_bitstring}){
                foreach my $file (@{$CFG->{input_files}}){
                    my $graphID = $file;
                    if(exists $tp->{omega}->{$graphID}){
                        $matched_graphs_bitstring .= "1";
                    } else {
                        $matched_graphs_bitstring .= "0";
                    }
                }
            }
            my $omega = {};
            if($CFG->{output_mappings}){
                foreach my $graphID (keys %{$tp->{omega}}){
                    my $graph_name = $graphID;
                    foreach my $mu (@{$tp->{omega}->{$graphID}}){
                        my $new_mu = {};
                        foreach my $var (keys %{$mu}){
                            $new_mu->{$var} = $mu->{$var};
                        }
                        push(@{$omega->{$graph_name}}, $new_mu);
                    }
                }
            }

            my $block = {};
            $block->{pattern_readable}         = &print_pattern_to_string($tp)                 if $CFG->{output_readable_patterns};
            $block->{pattern_parseable}         = $tp->{pattern_parseable}                     if $CFG->{output_parseable_patterns};
            $block->{list_of_matched_graphs}     = $list_of_matched_graphs                     if $CFG->{output_list_of_matched_graphs};
            $block->{matched_graphs_bitstring}     = $matched_graphs_bitstring                     if $CFG->{output_matched_graphs_bitstring};
            $block->{number_of_matched_graphs}     = scalar keys %{$tp->{omega}}                     if $CFG->{output_number_of_matched_graphs};
            $block->{omega}             = $omega                             if $CFG->{output_mappings};
            $block->{number_of_mappings}         = $number_of_mappings->{$tpID}                     if $CFG->{output_number_of_mappings};
            $TP->{$tpID} = $block;
        }
    }

    print "#patterns of size 1: " . (scalar keys %{$TP}) . "\n";

    # create an input file for each worker, distribute the patterns

    my $num_workers = $CFG->{number_of_workers};
    my $num_patterns = scalar keys %{$TP};

    my $cnt_patterns_per_worker = {};
    my $patterns_per_worker = {};

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        $cnt_patterns_per_worker->{$worker_ID} = $worker_ID == 1 ? int($num_patterns / $num_workers) + $num_patterns % $num_workers : int($num_patterns / $num_workers);
    }

    my $current_worker_ID = 1;
    foreach my $pattern_ID (sort keys %{$TP}){
        $current_worker_ID++ if scalar keys %{$patterns_per_worker->{$current_worker_ID}} == $cnt_patterns_per_worker->{$current_worker_ID};
        $patterns_per_worker->{$current_worker_ID}->{$pattern_ID} = 1;
    }


    my $single_output_file = sprintf($CFG->{incremental_output_file}, 1) . ".dat.gz.Z";
    my $output_all = new IO::Compress::Gzip($single_output_file)
        or die &get_timestamp() . " gzip failed: $GzipError\n";

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        my $output_file_name = sprintf($CFG->{incremental_input_file_per_worker}, 2, $worker_ID) . ".dat.gz.Z";

        print " > $output_file_name\n";
        my $output_single = new IO::Compress::Gzip($output_file_name)
            or die &get_timestamp() . " gzip failed: $GzipError\n";

        foreach my $pattern_ID (sort keys %{$patterns_per_worker->{$worker_ID}}){
            $output_single->write(&pattern_to_lines($TP->{$pattern_ID}, $pattern_ID));
            $output_all->write(&pattern_to_lines($TP->{$pattern_ID}, $pattern_ID));
        }

        close $output_single;
    }
    close $output_all;

    print "Written " . (scalar keys %{$TP}) . " patterns of size 1 to incremental output file <$single_output_file>\n";

    if(not scalar keys %{$TP}){
        system("touch " . sprintf($CFG->{incremental_output_file}, 1) . ".empty");
    }

    #print "wait\n"; <STDIN>;

    # TP can now be emptied. Unsure whether this happens automatically anyways.
    $TP = {};

} else {
    # TODO YAML might not be the best format
        $tpID_to_tp = LoadFile($CFG->{triplepattern_output_file});
}

#print "wait\n"; <STDIN>;



# TODO Warning about continuation. it should be continued with the same number of workers. otherwise, if previously more worker have been used, then some files will be omitted. if previously less workers have been used, then some workers will be idle in the first round.


# growth phase - part 2/2: iterative growth
my $pattern_size = 1;
if($continuation_is_possible){
    $pattern_size = $continuation_from_size -1;
}


my $patterns_were_created = 0;

while(1){

    print "\n\n\n[ v$CFG->{version} ] MAIN LOOP: Extending patterns of size $pattern_size\n";
    #print "  wait\n"; <STDIN>;

    last if defined $CFG->{max_pattern_size} and $pattern_size >= $CFG->{max_pattern_size};

    MCE::Loop->init(
        max_workers => $CFG->{number_of_workers},
        chunk_size => 1,
        user_begin => sub {
            print &get_timestamp() . " " . MCE->wid, " started\n";
        },
        user_end => sub {
            print &get_timestamp() . " " . MCE->wid, " completed\n";
        }
    );

    my %RES = mce_loop {
        my ($mce, $chunk_ref, $worker_id) = @_;
        #MCE->say(&get_timestamp() . " worker $worker_id started");

        my $input_file_name = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+1, $worker_id) . ".dat.gz.Z";
        my $output_file_name = sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_id) . ".dat.gz.Z";

	my $index = 0;

        #MCE->say("worker #$worker_id reads from compressed $input_file_name");
        #MCE->say("worker #$worker_id writes to compressed $output_file_name");

        #MCE->say("worker $worker_id reads from $input_file_name");
        my $input_fh = IO::Uncompress::Gunzip->new($input_file_name)
            or MCE->say(&get_timestamp() . " worker $worker_id cannot open input file $input_file_name: gunzip-1 failed: $!\n");

        #MCE->say(&get_timestamp() . " worker $worker_id writes to $output_file_name");
            
        my $output_fh = new IO::Compress::Gzip($output_file_name)
            or MCE->say(&get_timestamp() . " worker $worker_id cannot open output file $output_file_name: gzip failed: $GzipError\n");

        while(defined(my $tuple = &load_next_pattern($input_fh, 1))){
            my ($gpID, $gp) = @{$tuple};

            #MCE->say("worker $worker_id has loaded graph pattern $gpID");
            
            FETP: foreach my $tpID (sort keys %{$tpID_to_tp}){ # { $number_of_mappings->{$a} <=> $number_of_mappings->{$b} || $a <=> $b }
    
                my $tp = $tpID_to_tp->{$tpID}; #print Dump { tp => $tp->{pattern} };
                
                my $cnt_extensions = 0;

                FEEXT: foreach my $mj_mo (&ExtendPattern($gp, $tp, $index, $CFG->{random_numbers_1})){
                    my ($m_j, $m_o) = @{$mj_mo};
                    my $tp_dash = &substitute_variables_in_pattern($tp->{pattern_parseable}, $m_j, $m_o);
                    my $gp_ext = {
                        pattern_parseable => clone($gp->{pattern_parseable}),
                        omega => {},
                        number_of_mappings => 0
                    };
                    push(@{$gp_ext->{pattern_parseable}}, clone($tp_dash->[0]));

                    if(scalar keys %{$CFG->{undesired_subgraphs_tpcount}} and &contains_undesired_subgraph($gp_ext->{pattern_parseable})){
                        MCE->say("Ignored an extension that contains an undesired subgraph");
                        next FEEXT;
                    }

                    foreach my $input_file (sort keys %{$gp->{omega}}){
                        if(exists $tpID_to_tp->{$tpID}->{omega}->{$input_file}){
                            foreach my $mu_1 (@{$gp->{omega}->{$input_file}}){
                                foreach my $mu_2 (@{$tpID_to_tp->{$tpID}->{omega}->{$input_file}}){
                                    my $mu_2_dash = &substitute_variables_in_mapping($mu_2, $m_j, $m_o);
                                    
                                    if(&compatible($mu_1, $mu_2_dash) and &embedding_is_overlap_free($gp_ext->{pattern_parseable}, $mu_1, $mu_2_dash)){ # TODO: use this for p-confidentiality code, too
                                        my $mu_3 = {};
                                        $mu_3->{$_} = $mu_1->{$_} foreach keys %{$mu_1};
                                        $mu_3->{$_} = $mu_2_dash->{$_} foreach keys %{$mu_2_dash};

                                        push(@{$gp_ext->{omega}->{$input_file}}, $mu_3);
                                        $gp_ext->{number_of_mappings}++;

                                        if(
                                            defined $CFG->{max_matches_per_graph} and
                                            scalar @{$gp_ext->{omega}->{$input_file}} >= $CFG->{max_matches_per_graph}
                                        ){
                                            MCE->say("next-1");
                                            next FEEXT;
                                        }

                                        if(
                                            exists $CFG->{max_matches_per_graph_per_patternsize}->{$pattern_size + 1} and
                                            defined    $CFG->{max_matches_per_graph_per_patternsize}->{$pattern_size + 1} and
                                            scalar @{$gp_ext->{omega}->{$input_file}} >= $CFG->{max_matches_per_graph_per_patternsize}->{$pattern_size + 1}
                                        ){
                                            #print "skipped pattern because of max_matches_per_graph_per_patternsize\n"; <somoSTDIN>;
                                            MCE->say("next-2");
                                            next FEEXT;
                                        }
                                    }
                                }
                            }
                            if(    exists $gp_ext->{omega}->{$input_file} and
                                defined $CFG->{min_matches_per_graph} and
                                scalar @{$gp_ext->{omega}->{$input_file}} < $CFG->{min_matches_per_graph}
                            ){
                                MCE->say("next-3");
                                next FEEXT;
                            }
                        }

                    }

                    if(
                        (defined $CFG->{min_number_of_matched_graphs_abs} and scalar keys %{$gp_ext->{omega}} >= $CFG->{min_number_of_matched_graphs_abs})
                            or not defined $CFG->{min_number_of_matched_graphs_abs}
                    ){

                        if(
                            (defined $CFG->{max_number_of_matched_graphs_abs} and scalar keys %{$gp_ext->{omega}} < $CFG->{max_number_of_matched_graphs_abs})
                                or not defined $CFG->{max_number_of_matched_graphs_abs}
                        ){
                            my $fingerprint = &create_fingerprint($gp_ext);
                            $gp_ext->{fingerprint} = $fingerprint;

                            $cnt_extensions++;
                            my $gp_ext_ID = "(w$worker_id-gp$gpID-tp$tpID-num$cnt_extensions)";

                            #MCE->say("worker $worker_id writes a pattern to incremental output_file $output_file_name");
                            &write_to_output_file_handle($output_fh, $gp_ext_ID, $gp_ext);

                            #if($CFG->{verbose}){
                            #MCE->say("new pattern:\n");
                            #    MCE->say(&print_pattern_to_string($gp_ext));
                            #}
                            
                            delete $gp_ext->{omega};
                            MCE->gather($gp_ext_ID, $gp_ext);
                        }
                    }
                }
            }
        } # while
        $input_fh->close()
            or MCE->say("worker $worker_id cannot close input file $input_file_name: $!");
        $output_fh->close()
            or MCE->say("worker $worker_id cannot close output file $output_file_name: $!");

    } (1..$CFG->{number_of_workers}); # the MCE loop

    MCE::Loop->finish;

    print &get_timestamp() . " finished loop\n";

    # now, identify isomorphic patterns, rebalance output files while creating input files, remove isomorphic patterns at the same time

    print &get_timestamp() . " BEG iso, balance, shuffle\n";

    # first create fingerprint_to_pattern
    my $fingerprint_to_patternIDs;
    foreach my $pattern_ID (keys %RES){
        my $pattern = $RES{$pattern_ID};
        if(not exists $fingerprint_to_patternIDs->{$pattern->{fingerprint}}){
            $fingerprint_to_patternIDs->{$pattern->{fingerprint}} = [$pattern_ID];
        } else {
            push (@{$fingerprint_to_patternIDs->{$pattern->{fingerprint}}}, $pattern_ID);
        }
    }

    my $isomorphic_patterns = {};
    my $non_isomorphic_patterns = {};
    my $keep_isomorphic_patterns = {};

        print &get_timestamp() . " BEG iso\n";

    #print Dump { fingerprint_to_patternIDs => $fingerprint_to_patternIDs };

    my $cnt_iso_test = 0;

    foreach my $fingerprint (sort keys %{$fingerprint_to_patternIDs}){
        if(scalar @{$fingerprint_to_patternIDs->{$fingerprint}} == 1){
            my $pattern_ID = $fingerprint_to_patternIDs->{$fingerprint}->[0];
            $non_isomorphic_patterns->{$pattern_ID} = 1;
        } else {
            #print &get_timestamp() . " with fingerprint $fingerprint\n";
            foreach my $pattern_1_ID (sort @{$fingerprint_to_patternIDs->{$fingerprint}}){
                next if exists $isomorphic_patterns->{$pattern_1_ID};

                #print &get_timestamp() . "  with pattern_1_ID $pattern_1_ID\n";

                my $pattern_1 = $RES{$pattern_1_ID};

                my $gp1_var_char = {};
                foreach my $tp (@{$pattern_1->{pattern_parseable}}){
                    my $tpID = $tp->[3];
                    $gp1_var_char->{$tp->[0]}->{$tpID}->{s}++ if $tp->[0] =~ m/\A\?/;
                    $gp1_var_char->{$tp->[1]}->{$tpID}->{p}++ if $tp->[1] =~ m/\A\?/;
                    $gp1_var_char->{$tp->[2]}->{$tpID}->{o}++ if $tp->[2] =~ m/\A\?/;
                }

                my $an_isomorphic_pattern_exists = 0;

                foreach my $pattern_2_ID (sort @{$fingerprint_to_patternIDs->{$fingerprint}}){
                    next if $pattern_1_ID eq $pattern_2_ID;
                    next if $isomorphic_patterns->{$pattern_1_ID}->{$pattern_2_ID};
                    #print &get_timestamp() . "   with pattern_2_ID $pattern_2_ID\n";
                    my $pattern_2 = $RES{$pattern_2_ID};

                    $cnt_iso_test++;
                    print &get_timestamp() . " cnt_iso_test $cnt_iso_test\n" if $cnt_iso_test % 100000 == 0;
                    #print &get_timestamp() . "     <iso-test ...\n";
                    if(&patterns_are_isomorphic($pattern_1, $gp1_var_char, $pattern_2)){
                        $isomorphic_patterns->{$pattern_1_ID}->{$pattern_2_ID} = 1;
                        $isomorphic_patterns->{$pattern_2_ID}->{$pattern_1_ID} = 1;
                        $an_isomorphic_pattern_exists = 1;
                        $keep_isomorphic_patterns->{$pattern_1_ID} = 1;
                    }
                    #print &get_timestamp() . "     ... iso-test>\n";
                }

                if(not $an_isomorphic_pattern_exists){
                    $non_isomorphic_patterns->{$pattern_1_ID} = 1;
                }
            }
        }
    }

    print Dump {
        timestamp => &get_timestamp(),
        cnt_initially => scalar keys %RES,
        cnt_iso_test => $cnt_iso_test,
        cnt_non_isomorphic_patterns => scalar keys %{$non_isomorphic_patterns},
        cnt_keep_isomorphic_patterns => scalar keys %{$keep_isomorphic_patterns},
        cnt_all => scalar keys %RES,
    }; #<STDIN>;

    print &get_timestamp() . " END iso\n";

    my $num_workers = $CFG->{number_of_workers};
        my $num_patterns = (scalar keys %{$non_isomorphic_patterns}) + (scalar keys %{$keep_isomorphic_patterns});

    $patterns_were_created = $num_patterns > 0 ? 1 : 0;

    print &get_timestamp() . " #patterns of size " . ($pattern_size + 1) . ": $num_patterns\n";
    
        my $target_cnt_patterns_per_worker = {};

        for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
                $target_cnt_patterns_per_worker->{$worker_ID} = $worker_ID == 1 ? int($num_patterns / $num_workers) + $num_patterns % $num_workers : int($num_patterns / $num_workers);
        }

        my $cnt_patterns_per_worker = {};
    
    my $current_output_ID = 1;


    # shuffle while rebalancing the files
    my $number_of_random_numbers = scalar(@{$CFG->{random_numbers_2}});
    #print "number_of_random_numbers $number_of_random_numbers\n";
    #print Dump { random_numbers_2 => $CFG->{random_numbers_2} };
    my $index = 0;

    my $cnt_patterns = 0;
    my $output_fhs = {};
        for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        my $filename = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat.gz.Z";
        my $fh = new IO::Compress::Gzip($filename) or die &get_timestamp() .  " gzip failed: $GzipError\n";
        $output_fhs->{$worker_ID} = $fh;
    }

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        my $filename = sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat.gz.Z";
        my $input_fh = IO::Uncompress::Gunzip->new($filename)
                        or die &get_timestamp() . " gunzip-1 failed: $!\n(open $filename)";

        while(defined(my $tuple = &load_next_pattern($input_fh, 0))){
            my ($gpID, $gp) = @{$tuple};
            $cnt_patterns++;
            next if not (exists $non_isomorphic_patterns->{$gpID} or exists $keep_isomorphic_patterns->{$gpID});


	    # added: deterministic shuffling
            my @numbers = (1..$CFG->{number_of_workers});
	    my $number;
	    if($CFG->{deterministic_sampling}){
           	$number = 1 + $CFG->{random_numbers_2}->[$index++ % $number_of_random_numbers];
		#print "number = $number ($index)\n"; 
	    } else {
            	$number = $numbers[rand @numbers];
	    }
            my $random_output_fh = $output_fhs->{$number};

            &write_to_output_file_handle($random_output_fh, $gpID, $gp);
        }

        close $input_fh;
    }

    close($output_fhs->{$_}) foreach keys %{$output_fhs};
    
    if(not $cnt_patterns){
                system("touch " . sprintf($CFG->{incremental_output_file}, $pattern_size + 1) . ".empty");
        for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
                    my $filename = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat.gz.Z";
                    unlink $filename;
        }
    }

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
                my $filename1 = sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat.gz.Z";
        unlink $filename1;

        my $filename2 = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat.gz.Z";
                unlink $filename2;
    }

    if($cnt_patterns){
        my $output_file = sprintf($CFG->{incremental_output_file}, $pattern_size+1) . ".dat";
        for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
            if(not -s sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat.gz.Z"){
                my $name = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat.gz.Z";
                system("touch $name");
                if($? != 0){
                    print "touch-1 $! [$?]\n(touch $name)\n";
                }
            }
            my $part_file = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat.gz.Z";
            system("zcat $part_file >> $output_file");
            if($? != 0){
                print "zcat-1 $!\n(zcat $part_file >> $output_file)\n";
            }

            unlink sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat.gz.Z";
        }

        system("gzip -f $output_file");
        if($? != 0){
            print "gzip-2 $! [$?]\n(gzip -f $output_file)\n";
        }

        for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
            close $output_fhs->{$worker_ID};
        }
    }

        print &get_timestamp() . " END iso, balance, shuffle\n";



    #open(my $output_fh, '>' . sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $current_output_ID)  . ".dat");
    #print " > " . sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $current_output_ID)  . ".dat\n";
    
    #for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
    #    #print "with input worker $worker_ID\n";
    #    #print " < " . sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat\n";
    #    open(my $input_fh, '<' . sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat");

    #    while(defined(my $tuple = &load_next_pattern($input_fh))){
    #        my ($gpID, $gp) = @{$tuple};
    #        next if not (exists $non_isomorphic_patterns->{$gpID} or exists $keep_isomorphic_patterns->{$gpID});

    #        #print "writing gp from $worker_ID to $current_output_ID using write_to_output_file_handle\n";
    #        &write_to_output_file_handle($output_fh, $gpID, $gp);
    #        $cnt_patterns_per_worker->{$current_output_ID}++;
            
    #        if($cnt_patterns_per_worker->{$current_output_ID} == $target_cnt_patterns_per_worker->{$current_output_ID}){
    #            $current_output_ID++;
    #            close($output_fh);
    #            if($current_output_ID <= $CFG->{number_of_workers}){
    #                open($output_fh, '>' . sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $current_output_ID)  . ".dat");
    #                print " > " . sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $current_output_ID)  . ".dat\n";
    #            }
    #        }
    #    }
    #    if(not -s sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat"){
    #        system("touch " . sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat");
    #    }
    #}


    #my $output_file = sprintf($CFG->{incremental_output_file}, $pattern_size+1) . ".dat";
    #for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
    #    my $part_file = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+2, $worker_ID)  . ".dat";
    #    system("cat $part_file >> $output_file");
    #}
    #system("zip $output_file.zip $output_file");
    #unlink $output_file;

    #for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
    #    unlink sprintf($CFG->{incremental_output_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat";
    #    my $filename = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+1, $worker_ID)  . ".dat";
    #    system("zip $filename.zip $filename");
    #    unlink $filename;
    #}
    
    $pattern_size++;

    last if $num_patterns == 0;

} # while


if($patterns_were_created){
    my $file_pattern = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+1, "123456789")  . ".dat.gz.Z"; # TODO bad hack
    $file_pattern =~ s/123456789/*/; # TODO bad hack
    my $tar_file_name = sprintf($CFG->{incremental_input_file_per_worker}, $pattern_size+1, "123456789")  . ".dat.tar.gz"; # TODO bad hack
    $tar_file_name =~ s/123456789/all/; # TODO bad hack
    system("gtar czf $tar_file_name $file_pattern --remove-files 2>/dev/null");
    if($? != 0){
        print "tar-3 command: $! [$?]\n(tar czf -C / $tar_file_name $file_pattern)\n";
    }
}


if(exists $CFG->{pattern_visualization} and $CFG->{pattern_visualization}->{visualize_patterns}){
    print "pattern_visualization\n";
    #print "<STDIN>\n"; <STDIN>;
    for(my $i=$CFG->{min_pattern_size}; $i<=$CFG->{max_pattern_size}; $i++){
        print " pattern_size $i\n";
        #next if $i < 25; # TODO remove
        my $file = sprintf($CFG->{incremental_output_file}, $i) . ".dat.gz";
        my $dir = sprintf($CFG->{pattern_visualization}->{visualization_directory}, $i);
        
        print "  file: $file\n";
        # Warning: this might visualize patterns obtained in previous runs
        if(-e $file){

            print "Visualizing patterns of size $i ...";
            
            mkdir $dir if not -d $dir;# and scalar keys %{$DAT->{results}->{$i}};
        
            DumpFile("$dir/config.cfg", {
                config => $CFG,
            });

            my $input_fh = IO::Uncompress::Gunzip->new($file)
                        or die &get_timestamp() . " gunzip-1 failed: $!\n(open $file)";

            while(defined(my $tuple = &load_next_pattern($input_fh, 0))){
                my ($patternID, $pattern) = @{$tuple};

                my $number_of_mappings = 0;
                my $number_of_matched_graphs = scalar keys %{$pattern->{omega}};
                foreach my $graphname (keys %{$pattern->{omega}}){
                    $number_of_mappings += scalar @{$pattern->{omega}->{$graphname}};
                }

                my $patternname = $patternID;
                $patternname =~ s/\(//g;
                $patternname =~ s/\)//g;

                my $introduced_nodes = {};
                my $last_node_id = -1;
                open (OUT,">$dir/pattern-$patternname.dot");
                print OUT "digraph G {\n";
                print OUT "\t node [ shape=\"none\" ];\n";
                foreach my $tp (@{$pattern->{pattern_parseable}}){
                    if(not exists $introduced_nodes->{$tp->[0]}){
                        $last_node_id++;
                        $introduced_nodes->{$tp->[0]} = "node$last_node_id";
                        print OUT "\tnode$last_node_id [label=\"" . &shorten($tp->[0]) . "\"];\n";
                    }

                    if(not exists $introduced_nodes->{$tp->[2]}){
                        $last_node_id++;
                        $introduced_nodes->{$tp->[2]} = "node$last_node_id";
                        my $nodelabel;
                        if($tp->[2] =~ m/\A"(.*)"\^\^(<.*>)\Z/){
                            $nodelabel = "\"\\\"" . $1 . "\\\"" . "\\\^\\\^" . &shorten($2) . "\"";
                        } elsif($tp->[2] =~ m/\A"(.*)"\Z/){
                            $nodelabel = "\"\\\"$1\\\"\"";
                        } elsif($tp->[2] =~ m/\A"(.*)"\@(.*)\Z/){ # TODO: might not be safe
                            $nodelabel = "\"\\\"$1\\\"\\\@$2\"";
                        } else {
                            $nodelabel = "\"" . &shorten($tp->[2]) . "\"";
                        }
                        print OUT "\tnode$last_node_id [label=$nodelabel];\n";
                    }

                    print OUT "\t" . $introduced_nodes->{$tp->[0]} . " -> " . $introduced_nodes->{$tp->[2]} . " [label=\"" . &shorten($tp->[1]). "\"];\n";
                }
                print OUT "\t labelloc=\"b\";\n";
                print OUT "\t label=\"patternID: $patternID, size: $i, number_of_mappings: "
                . $number_of_mappings
                . ", number_of_matched_graphs: "
                . $number_of_matched_graphs
                . "\";\n";

                print OUT "}\n";
                close OUT;

                system("dot -T" . $CFG->{pattern_visualization}->{filetype}. " -K" . $CFG->{pattern_visualization}->{layout} . " $dir/pattern-$patternname.dot -o $dir/pattern-$patternname.png");
                unlink "$dir/pattern-$patternname.dot";
                #print "wait\n"; <STDIN>;
            }
            close $input_fh;

            print " done\n";
        } else { print "file <$file> does not exist\n"; }
    }
}

###############################################################################

sub load_next_pattern {
        my ($p1_fh, $in_MCE) = @_;

        my $pattern_ID = undef;
        my $p = {};
        my $g = undef;
        my $m = {};

        while(defined(my $line = <$p1_fh>)){
                if($line =~ m/\A0 (.*)\n/){
                        $p = {};
                        $pattern_ID = $1;
                }
                elsif($line =~ m/\A1 (.*)\n/){
                        $p->{number_of_mappings} = $1;
                }
        #elsif($line =~ m/\A1 (.*)\n/){
        #       $p->{number_of_mappings} = $1;
        #}
                elsif($line =~ m/\A2 (.*)\n/){
                        $p->{number_of_matched_graphs} = $1;
                }
                elsif($line =~ m/\A3 (.*)\n/){
                        push(@{$p->{list_of_matched_graphs}}, $1);
                }
                elsif($line =~ m/\A4 (.*)\n/){
                        $p->{matched_graphs_bitstring} = $1;
                }
                elsif($line =~ m/\A5 (.*)\n/){
                        $p->{pattern_readable} = $1;
                }
                elsif($line =~ m/\A6 (.*)\t(.*)\t(.*)\t(.*)\n/){
                        push(@{$p->{pattern_parseable}}, [$1, $2, $3, $4]);
                }
                elsif($line =~ m/\A7 (.*)\n/){
                        $g = $1;
                        $m = {};
                }
                elsif($line =~ m/\A8 (.*)\t(.*)\n/){
                        $m->{$1} = $2;
                }
                elsif($line =~ m/\A9\n/){
                        push(@{$p->{omega}->{$g}}, clone($m));
                        $m = {};
                }
                elsif($line =~ m/\A10\n/){
            #print Dump { loaded => { pID => $pattern_ID, p => $p } }; <STDIN>;
                        return [$pattern_ID, $p];
                }
                else {
            if($in_MCE){
                MCE->say("line format unexpected: \n$line");
            } else {
                            print "line format unexpected: \n$line\n"; <STDIN>;
            }
                }
        }
        return undef;
}

sub pattern_to_lines {
    my ($p, $pattern_ID) = @_;

    my $lines = q{};

    $lines .= "0 $pattern_ID\n";
    $lines .= "1 " . $p->{number_of_mappings} . "\n"
        if exists $p->{number_of_mappings};

    $lines .= "2 " . $p->{number_of_matched_graphs} . "\n"
        if exists $p->{number_of_matched_graphs};

    if(exists $p->{list_of_matched_graphs}){
        foreach my $graph_ID (@{$p->{list_of_matched_graphs}}){
            $lines .= "3 $graph_ID\n";
        }
    }

    if(exists $p->{matched_graphs_bitstring}){
        $lines .= "4 " . $p->{matched_graphs_bitstring} . "\n";
    }

    if(exists $p->{pattern_readable}){
        $lines .= "5 " . $p->{pattern_readable} . "\n";
    }

    if(exists $p->{pattern_parseable}){
        foreach my $tp (@{$p->{pattern_parseable}}){
            $lines .= "6 " . join("\t", @{$tp}) . "\n";
        }
    }

    if(exists $p->{omega}){
        foreach my $graph_ID (sort keys %{$p->{omega}}){
            $lines .= "7 $graph_ID\n";
            foreach my $m (@{$p->{omega}->{$graph_ID}}){
                foreach my $v (sort keys %{$m}){
                    $lines .= "8 $v\t" . $m->{$v} . "\n";
                }
                $lines .= "9\n";
            }
        }
    }
    $lines .= "10\n";

    return $lines;
}

sub write_to_n_input_files {
    my ($next_size, $CFG, $patterns) = @_;

    my $num_workers = $CFG->{number_of_workers};
    my $num_patterns = scalar keys %{$patterns};
    
    my $cnt_patterns_per_worker = {};
    my $patterns_per_worker = {};

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        $cnt_patterns_per_worker->{$worker_ID} = $worker_ID == 1 ? int($num_patterns / $num_workers) + $num_patterns % $num_workers : int($num_patterns / $num_workers);
    }

    #print Dump { cnt_patterns_per_worker => $cnt_patterns_per_worker, num_workers => $num_workers, num_patterns => $num_patterns }; <STDIN>;

    my $current_worker_ID = 1;
    foreach my $pattern_ID (sort keys %{$patterns}){
        $current_worker_ID++ if scalar keys %{$patterns_per_worker->{$current_worker_ID}} == $cnt_patterns_per_worker->{$current_worker_ID};
        $patterns_per_worker->{$current_worker_ID}->{$pattern_ID} = 1;
    }

    #print Dump { patterns_per_worker => $patterns_per_worker }; <STDIN>;

    for(my $worker_ID = 1; $worker_ID <= $CFG->{number_of_workers}; $worker_ID++){
        my $output_file_name = sprintf($CFG->{incremental_input_file_per_worker}, $next_size, $worker_ID) . ".dat";

        print " > $output_file_name\n";
        open(OUT, '>' . $output_file_name);

        foreach my $pattern_ID (sort keys %{$patterns_per_worker->{$worker_ID}}){
            #print "file $worker_ID: add $pattern_ID\n";

            my $p = $patterns->{$pattern_ID};

            print OUT "0 $pattern_ID\n";
            print OUT "1 " . $p->{number_of_mappings} . "\n"
                if exists $p->{number_of_mappings};

            print OUT "2 " . $p->{number_of_matched_graphs} . "\n"
                if exists $p->{number_of_matched_graphs};

            if(exists $p->{list_of_matched_graphs}){
                foreach my $graph_ID (@{$p->{list_of_matched_graphs}}){
                    print OUT "3 $graph_ID\n";
                }
            }

            if(exists $p->{matched_graphs_bitstring}){
                print OUT "4 " . $p->{matched_graphs_bitstring} . "\n";
            }

            if(exists $p->{pattern_readable}){
                print OUT "5 " . $p->{pattern_readable} . "\n";
            }

            if(exists $p->{pattern_parseable}){
                foreach my $tp (@{$p->{pattern_parseable}}){
                    print OUT "6 " . join("\t", @{$tp}) . "\n";
                }
            }

            if(exists $p->{omega}){
                foreach my $graph_ID (sort keys %{$p->{omega}}){
                    print OUT "7 $graph_ID\n";
                    foreach my $m (@{$p->{omega}->{$graph_ID}}){
                        foreach my $v (sort keys %{$m}){
                            print OUT "8 $v\t" . $m->{$v} . "\n";
                        }
                        print OUT "9\n";
                    }
                }
            }
            print OUT "10\n";
        }

        close OUT;
    }
}

sub write_to_output_file {
    my ($output_file_name, $CFG, $patterns) = @_;

    #print " > $output_file_name.yml\n";
    #DumpFile("${output_file_name}.yml", $CFG); # TODO. unclear why this causes an error

    print " > $output_file_name.dat\n";
    open(OUT,">$output_file_name.dat");
    foreach my $pattern_ID (sort { $a <=> $b } keys %{$patterns}){
        my $p = $patterns->{$pattern_ID};
        print OUT "0 $pattern_ID\n";
        print OUT "1 " . $p->{number_of_mappings} . "\n"
            if exists $p->{number_of_mappings};

        print OUT "2 " . $p->{number_of_matched_graphs} . "\n"
            if exists $p->{number_of_matched_graphs};

        if(exists $p->{list_of_matched_graphs}){
            foreach my $graph_ID (@{$p->{list_of_matched_graphs}}){
                print OUT "3 $graph_ID\n";
            }
        }

        if(exists $p->{matched_graphs_bitstring}){
            print OUT "4 " . $p->{matched_graphs_bitstring} . "\n";
        }

        if(exists $p->{pattern_readable}){
            print OUT "5 " . $p->{pattern_readable} . "\n";
        }

        if(exists $p->{pattern_parseable}){
            foreach my $tp (@{$p->{pattern_parseable}}){
                print OUT "6 " . join("\t", @{$tp}) . "\n";
            }
        }

        if(exists $p->{omega}){
            foreach my $graph_ID (sort keys %{$p->{omega}}){
                print OUT "7 $graph_ID\n";
                foreach my $m (@{$p->{omega}->{$graph_ID}}){
                    foreach my $v (sort keys %{$m}){
                        print OUT "8 $v\t" . $m->{$v} . "\n";
                    }
                    print OUT "9\n";
                }
            }
        }
        print OUT "10\n";
    }
    close OUT;
}

sub write_to_output_file_handle {
    my ($output_fh, $pattern_ID, $pattern) = @_;

    $output_fh->write("0 $pattern_ID\n");
    
    $output_fh->write("1 " . $pattern->{number_of_mappings} . "\n")
        if exists $pattern->{number_of_mappings};

    $output_fh->write("2 " . $pattern->{number_of_matched_graphs} . "\n")
        if exists $pattern->{number_of_matched_graphs};

    if(exists $pattern->{list_of_matched_graphs}){
        foreach my $graph_ID (@{$pattern->{list_of_matched_graphs}}){
            $output_fh->write("3 $graph_ID\n");
        }
    }

    if(exists $pattern->{matched_graphs_bitstring}){
        $output_fh->write("4 " . $pattern->{matched_graphs_bitstring} . "\n");
    }

    if(exists $pattern->{pattern_readable}){
        $output_fh->write("5 " . $pattern->{pattern_readable} . "\n");
    }

    if(exists $pattern->{pattern_parseable}){
        foreach my $tp (@{$pattern->{pattern_parseable}}){
            $output_fh->write("6 " . join("\t", @{$tp}) . "\n");
        }
    }

    if(exists $pattern->{omega}){
        foreach my $graph_ID (sort keys %{$pattern->{omega}}){
            $output_fh->write("7 $graph_ID\n");
            foreach my $m (@{$pattern->{omega}->{$graph_ID}}){
                foreach my $v (sort keys %{$m}){
                    $output_fh->write("8 $v\t" . $m->{$v} . "\n");
                }
                $output_fh->write("9\n");
            }
        }
    }
    $output_fh->write("10\n");
}


sub patterns_are_isomorphic { #print "patterns_are_isomorphic?\n";
    my ($gp, $gp_var_char, $can) = @_;

    my $can_var_char = {};
    foreach my $tp (@{$can->{pattern_parseable}}){
        my $tpID = $tp->[3];
        $can_var_char->{$tp->[0]}->{$tpID}->{s}++ if $tp->[0] =~ m/\A\?/;
        $can_var_char->{$tp->[1]}->{$tpID}->{p}++ if $tp->[1] =~ m/\A\?/;
        $can_var_char->{$tp->[2]}->{$tpID}->{o}++ if $tp->[2] =~ m/\A\?/;
    }

    my $gp_var_char_flat = {}; # TODO this should be done elsewhere, so that the same computation is not repeated
    foreach my $v (sort keys %{$gp_var_char}){
        $gp_var_char_flat->{$v} = q{};
        foreach my $tpID (sort { $a <=> $b } keys %{$gp_var_char->{$v}}){
            foreach my $pos (qw(s p o)){
                if(exists $gp_var_char->{$v}->{$tpID}->{$pos}){
                    $gp_var_char_flat->{$v} .= "$tpID $pos " . $gp_var_char->{$v}->{$tpID}->{$pos} . " ";
                }
            }
        }
        $gp_var_char_flat->{$v} .= ".";
    }

    my $can_var_char_flat = {};
        foreach my $v (sort keys %{$can_var_char}){
                $can_var_char_flat->{$v} = q{};
                foreach my $tpID (sort { $a <=> $b } keys %{$can_var_char->{$v}}){
                        foreach my $pos (qw(s p o)){
                                if(exists $can_var_char->{$v}->{$tpID}->{$pos}){
                                        $can_var_char_flat->{$v} .= "$tpID $pos " . $can_var_char->{$v}->{$tpID}->{$pos} . " ";
                                }
                        }
                }
                $can_var_char_flat->{$v} .= ".";
        }


    my $posmap_t = {};
    my $size = scalar @{$gp->{pattern_parseable}} -1;

    foreach my $tp1_num (0..$size){
        my $tp1 = $gp->{pattern_parseable}->[$tp1_num];
        my $tp1_ID = $tp1->[3];
        FETP2: foreach my $tp2_num (0..$size){
            my $tp2 = $can->{pattern_parseable}->[$tp2_num];
            my $tp2_ID = $tp2->[3];
            next if $tp1_ID ne $tp2_ID;

            foreach my $posnum (qw(0 1 2)){
                if($tp1->[$posnum] =~ m/\A\?/){
                    my ($v1, $v2) = ($tp1->[$posnum], $tp2->[$posnum]);
                    next FETP2 if $gp_var_char_flat->{$v1} ne $can_var_char_flat->{$v2};

                    #foreach my $tpID (keys %{$gp_var_char->{$v1}}){
                    #    next FETP2 if not exists $can_var_char->{$v2}->{$tpID};
                    #    foreach my $pos (qw(s p o)){
                    #        next FETP2 if exists $gp_var_char->{$v1}->{$tpID}->{$pos} and not exists $can_var_char->{$v2}->{$tpID}->{$pos};
                    #        next FETP2 if exists $gp_var_char->{$v1}->{$tpID}->{$pos} and $gp_var_char->{$v1}->{$tpID}->{$pos} ne $can_var_char->{$v2}->{$tpID}->{$pos};
                    #    }
                    #}
                }
            }

            $posmap_t->{$tp1_num}->{$tp2_num} = 1;
        }
    }

    return 0 if scalar keys %{$posmap_t} ne scalar @{$gp->{pattern_parseable}};

    return 0 if not scalar keys %{$posmap_t};

    return 1 if &suitable_injection_exists($posmap_t, $gp->{pattern_parseable}, $can->{pattern_parseable});

    return 0;
}

sub subgraph_isomorphism_exists {
    my ($small, $large, $small_var_char, $large_var_char) = @_;

    my $posmap_t = {};
    for(my $pos1=0; $pos1<scalar @{$small}; $pos1++){
        for(my $pos2=0; $pos2<scalar @{$large}; $pos2++){
            if($small->[$pos1]->[3] == $large->[$pos2]->[3]){

                # TODO: look at how it is done in isomorphic_pattern_exists
                if($small->[$pos1]->[0] =~ m/\A\?/){
                    # compare characteristics of this var against corresponding var in corresponding tp
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{s} and not exists $large_var_char->{$large->[$pos2]->[0]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{p} and not exists $large_var_char->{$large->[$pos2]->[0]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{o} and not exists $large_var_char->{$large->[$pos2]->[0]}->{o};
                    
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{s} and $small_var_char->{$small->[$pos1]->[0]}->{s} > $large_var_char->{$large->[$pos2]->[0]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{p} and $small_var_char->{$small->[$pos1]->[0]}->{p} > $large_var_char->{$large->[$pos2]->[0]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[0]}->{o} and $small_var_char->{$small->[$pos1]->[0]}->{o} > $large_var_char->{$large->[$pos2]->[0]}->{o};
                }
                if($small->[$pos1]->[1] =~ m/\A\?/){
                    # compare characteristics of this var against corresponding var in corresponding tp
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{s} and not exists $large_var_char->{$large->[$pos2]->[1]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{p} and not exists $large_var_char->{$large->[$pos2]->[1]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{o} and not exists $large_var_char->{$large->[$pos2]->[1]}->{o};
                    
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{s} and $small_var_char->{$small->[$pos1]->[1]}->{s} > $large_var_char->{$large->[$pos2]->[1]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{p} and $small_var_char->{$small->[$pos1]->[1]}->{p} > $large_var_char->{$large->[$pos2]->[1]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[1]}->{o} and $small_var_char->{$small->[$pos1]->[1]}->{o} > $large_var_char->{$large->[$pos2]->[1]}->{o};
                }
                if($small->[$pos1]->[2] =~ m/\A\?/){
                    # compare characteristics of this var against corresponding var in corresponding tp
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{s} and not exists $large_var_char->{$large->[$pos2]->[2]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{p} and not exists $large_var_char->{$large->[$pos2]->[2]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{o} and not exists $large_var_char->{$large->[$pos2]->[2]}->{o};
                    
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{s} and $small_var_char->{$small->[$pos1]->[2]}->{s} ne $large_var_char->{$large->[$pos2]->[2]}->{s};
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{p} and $small_var_char->{$small->[$pos1]->[2]}->{p} ne $large_var_char->{$large->[$pos2]->[2]}->{p};
                    next if exists $small_var_char->{$small->[$pos1]->[2]}->{o} and $small_var_char->{$small->[$pos1]->[2]}->{o} ne $large_var_char->{$large->[$pos2]->[2]}->{o};
                }

                $posmap_t->{$pos1}->{$pos2}++;
            }
        }
    }

    return &suitable_injection_exists($posmap_t, $small, $large);
}

sub suitable_injection_exists { #print "suitable_injection_exists\n";
    my ($posmap_t, $gp1, $gp2) = @_;

    my $F = [
        {
            t_gp1_to_t_gp2 => {},
            t_gp2_to_t_gp1 => {},
            v_gp1_to_v_gp2 => {},
            v_gp2_to_v_gp1 => {},
        }
    ];

    #print Dump {
    #    gp1 => $gp1,
    #    gp2 => $gp2,
    #    posmap_t => $posmap_t
    #}; <STDIN>;

    foreach my $key1 (sort { $a <=> $b } keys %{$posmap_t}){
        my $F_new = [];
        foreach my $key2 (sort { $a <=> $b } keys %{$posmap_t->{$key1}}){
            FEF: foreach my $f (@{$F}){
                next FEF if exists $f->{t_gp2_to_t_gp1}->{$key2};

                my $tp_gp1 = $gp1->[$key1];
                my $tp_gp2 = $gp2->[$key2];

                my $f_new = clone($f);
                $f_new->{t_gp1_to_t_gp2}->{$key1} = $key2;
                $f_new->{t_gp2_to_t_gp1}->{$key2} = $key1;

                foreach my $pos (qw(0 1 2)){
                    if($tp_gp1->[$pos] =~ m/\A\?/){
                        next FEF if exists $f_new->{v_gp1_to_v_gp2}->{$tp_gp1->[$pos]} and $f_new->{v_gp1_to_v_gp2}->{$tp_gp1->[$pos]} ne $tp_gp2->[$pos];
                        next FEF if exists $f_new->{v_gp2_to_v_gp1}->{$tp_gp2->[$pos]} and $f_new->{v_gp2_to_v_gp1}->{$tp_gp2->[$pos]} ne $tp_gp1->[$pos];
                        $f_new->{v_gp1_to_v_gp2}->{$tp_gp1->[$pos]} = $tp_gp2->[$pos];
                        $f_new->{v_gp2_to_v_gp1}->{$tp_gp2->[$pos]} = $tp_gp1->[$pos];
                    }
                }

                #print Dump { f_new => $f_new }; <STDIN>;

                push(@{$F_new}, $f_new);
                #if(
                #    (exists $f_new->{t_gp1_to_t_gp2}->{2} and not exists $f_new->{t_gp1_to_t_gp2}->{1}) or
                #    (exists $f_new->{t_gp1_to_t_gp2}->{3} and not exists $f_new->{t_gp1_to_t_gp2}->{2}) or
                #    (exists $f_new->{t_gp1_to_t_gp2}->{3} and not exists $f_new->{t_gp1_to_t_gp2}->{1})
                #){
                #    print Dump {
                #        f_new => $f_new,
                #        posmap_t => $posmap_t,
                #    }; <STDIN>;
                #}
            }
        
            if(not scalar @{$F_new}){
                return 0;
            }
        }
        $F = $F_new;
    }

    if(scalar @{$F} > 0){
        #print "suitable injection exists!\n";
        #print Dump { posmap_t => $posmap_t, gp1 => $gp1, gp2 => $gp2, F => $F }; <STDIN>;
        return 1;
    }

    return 0;
}

sub ExtendPattern {
    my ($gp, $tp, $index, $random_numbers) = @_;
    my @X = ();

    #print Dump { ExtendPattern => { gp => $gp->{pattern}, tp => $tp->{pattern} } };
    
    #print "gp:\n";
    #&print_pattern($gp);
    #print "tp:\n";
    #&print_pattern($tp);

    #MCE->say(Dump {random_numbers => $random_numbers }); die;
    my $number_of_random_numbers = scalar(@{$random_numbers});

    my $varnodes = {};
    my $nonvarnodes = {};

    my $V_gp = {};
    my $V_tp = {};
    my $Omega_gp = $gp->{omega};
    my $Omega_tp = $tp->{omega};

    my $gp_varpos = {};
    my $tp_varpos = {};

    foreach (@{$gp->{pattern_parseable}}){
        if($_->[0] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_gp->{$_->[0]} = 1;
            $gp_varpos->{$_->[0]}->{s} = 1;
            $varnodes->{$_->[0]} = 1;
        } else {
            $nonvarnodes->{$_->[0]} = 1;
        }
        if($_->[1] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_gp->{$_->[1]} = 1;
            $gp_varpos->{$_->[1]}->{p} = 1;
        }
        if($_->[2] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_gp->{$_->[2]} = 1;
            $gp_varpos->{$_->[2]}->{o} = 1;
            $varnodes->{$_->[2]} = 1;
        } else {
            $nonvarnodes->{$_->[2]} = 1;
        }
    }
    
    foreach (@{$tp->{pattern_parseable}}){
        if($_->[0] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_tp->{$_->[0]} = 1;
            $tp_varpos->{$_->[0]}->{s} = 1;
        }
        if($_->[1] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_tp->{$_->[1]} = 1;
            $tp_varpos->{$_->[1]}->{p} = 1;
        }
        if($_->[2] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
            $V_tp->{$_->[2]} = 1;
            $tp_varpos->{$_->[2]}->{o} = 1;
        }
    }
    
    my @tp_varlist = sort keys %{$V_tp};
    my @gp_varlist = sort keys %{$V_gp};

    my @tplists; # this could be done via Data::PowerSet, but the less dependencies the better.
    if(scalar keys %{$V_tp} == 1){
        # if the triple pattern contains exactly one variable, then this variable needs to be mapped to some gp variable.
        @tplists = ([0]);
    } elsif(scalar keys %{$V_tp} == 2){
        # if the triple pattern contains exactly two variables, then either the first or the second or both can be mapped.
        @tplists = ([0], [1], [0,1]);
    } elsif(keys %{$V_tp} == 3){
        @tplists = ([0], [1], [2], [0,1], [0,2], [1,2], [0,1,2]);
        # Explanation. [0,2] means that the tp-variable in first variable and the third variable
        # will be mapped to some gp-variable. In this case, as there are three variables,
        # it means that the first variable occurs in subject position and the third variable
        # occurs in object position of the triple pattern.
    } else {
        die "Something is really wrong. A triple pattern cannot have more than 3 variables.\n";
    }
    

    my @possibilities = ();
    foreach my $tplist (@tplists){
        my $mj_mo_set = [{m_j => {}, m_o => {}, used_gpv => {}}];
        foreach my $tp_n (@{$tplist}){
            my $new_mj_mo_set = [];
            my $tp_v = $tp_varlist[$tp_n];
            FEGPV: foreach my $gp_v (sort keys %{$V_gp}){
                foreach my $tp_v_pos (sort keys %{$tp_varpos->{$tp_v}}){
                    foreach my $gp_v_pos (sort keys %{$gp_varpos->{$gp_v}}){
                        my $join_type = $tp_v_pos . "-" . $gp_v_pos;
                        next FEGPV if not exists $CFG->{allowed_join_types}->{$join_type} or not $CFG->{allowed_join_types}->{$join_type};
                    }
                }
                FEMJMO: foreach my $mj_mo (@{$mj_mo_set}){
                    next FEMJMO if exists $mj_mo->{used_gpv}->{$gp_v};
                    my $new_mj_mo = clone($mj_mo);
                    $new_mj_mo->{m_j}->{$tp_v} = $gp_v;
                    $new_mj_mo->{used_gpv}->{$gp_v} = 1;
                    push (@{$new_mj_mo_set}, $new_mj_mo);
                }

            }
            $mj_mo_set = $new_mj_mo_set;
        }
        FEMJMO: foreach my $mj_mo (@{$mj_mo_set}){
            if(scalar keys %{$mj_mo->{m_j}} == scalar @{$tplist}){
                foreach my $tp_v (sort keys %{$V_tp}){
                                        if(not exists $mj_mo->{m_j}->{$tp_v}){
                                                $mj_mo->{m_o}->{$tp_v} = "?v" . ((scalar keys %{$mj_mo->{m_o}}) + (scalar keys %{$V_gp}));
                                        }
                                }

                my $tp_dash = &substitute_variables_in_pattern($tp->{pattern_parseable}, $mj_mo->{m_j}, $mj_mo->{m_o});

                if($tp_dash->[0]->[0] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
                    $varnodes->{$tp_dash->[0]->[0]} = 1;
                } else {
                    $nonvarnodes->{$tp_dash->[0]->[0]} = 1;
                }
                if($tp_dash->[0]->[2] =~ m/\A\?/){ # not very elegant, but probably fast. and memory-efficient.
                    $varnodes->{$tp_dash->[0]->[2]} = 1;
                } else {
                    $nonvarnodes->{$tp_dash->[0]->[2]} = 1;
                }


                # make sure that tp_dash does not already exist in gp:
                foreach my $gp_tp (@{$gp->{pattern_parseable}}){
                    if(
                        $gp_tp->[0] eq $tp_dash->[0]->[0] and
                        $gp_tp->[1] eq $tp_dash->[0]->[1] and
                        $gp_tp->[2] eq $tp_dash->[0]->[2]
                    ){
                        next FEMJMO;
                    }
                }


                my $cnt_varnodes = scalar keys %{$varnodes};
                my $cnt_nonvarnodes = scalar keys %{$nonvarnodes};

                if(exists $CFG->{max_share_of_varnodes} and defined $CFG->{max_share_of_varnodes}){
                    if(
                        exists $CFG->{start_patternsize_for_share_of_varnodes_constraint} and
                        defined $CFG->{start_patternsize_for_share_of_varnodes_constraint} and
                        (scalar @{$gp->{pattern_parseable}}) + 1 >= $CFG->{start_patternsize_for_share_of_varnodes_constraint}
                    ){
                        my $share_of_varnodes = $cnt_varnodes/($cnt_varnodes + $cnt_nonvarnodes);
                        if($share_of_varnodes >= $CFG->{max_share_of_varnodes}){
                            #print "skipped an extension due to the constraint max_share_of_varnodes\n";
                            #print Dump { varnodes => $varnodes, nonvarnodes => $nonvarnodes, cnt_varnodes => $cnt_varnodes, cnt_nonvarnodes => $cnt_nonvarnodes, share_of_varnodes => $share_of_varnodes }; <STDIN>;
                            next FEMJMO;
                        }
                    }
                }



                if(
                    (not defined $CFG->{extension_sample_threshold})
                    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        exists $CFG->{extension_sample_threshold}->{"any"} and
                        $CFG->{extension_sample_threshold}->{"any"} != 0 and
			not $CFG->{deterministic_sampling} and 
                        rand(1) >= $CFG->{extension_sample_threshold}->{"any"}
                    )
		    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        exists $CFG->{extension_sample_threshold}->{"any"} and
                        $CFG->{extension_sample_threshold}->{"any"} != 0 and
			not $CFG->{deterministic_sampling} and 
                        $CFG->{random_numbers_1}->[$index++ % $number_of_random_numbers] >= $CFG->{extension_sample_threshold}->{"any"}
                    )
                    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        not defined $CFG->{extension_sample_threshold}->{"any"} and
                        exists $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} and
                        $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} != 0 and
			not $CFG->{deterministic_sampling} and 			
                        rand(1) >= $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1}
                    )
		    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        not defined $CFG->{extension_sample_threshold}->{"any"} and
                        exists $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} and
                        $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} != 0 and
			$CFG->{deterministic_sampling} and 			
                        $CFG->{random_numbers_1}->[$index++ % $number_of_random_numbers] >= $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1}
                    )
                    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        not defined $CFG->{extension_sample_threshold}->{"any"} and
                        exists $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} and
                        $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1} == 0
                    )
                    or
                    (
                        defined $CFG->{extension_sample_threshold} and
                        not defined $CFG->{extension_sample_threshold}->{"any"} and
                        not exists $CFG->{extension_sample_threshold}->{scalar @{$gp->{pattern_parseable}} + 1}
                    )


                ){
                    push (@possibilities, [$mj_mo->{m_j}, $mj_mo->{m_o}]);
                }
                
            }
        }
        #print Dump { tplist => $tplist, mj_mo_set => $mj_mo_set }; <STDIN>;
        
    }

    #print Dump { possibilities => \@possibilities }; <STDIN>;
    return @possibilities;
}

sub compatible {
    my ($mu_1, $mu_2) = @_;

    #print Dump {
    #    mu_1 => $mu_1,
    #    mu_2 => $mu_2,
    #}; #<STDIN>;

    foreach my $v1 (sort keys %{$mu_1}){
        my $t1 = $mu_1->{$v1};
        foreach my $v2 (sort keys %{$mu_2}){
            my $t2 = $mu_2->{$v2};
            return 0 if $v1 eq $v2 and $t1 ne $t2;
            return 0 if $v1 ne $v2 and $t1 eq $t2;
        }
    }

    #print "return 1\n"; <STDIN>;
    return 1;
}

sub embedding_is_overlap_free {
    my ($pattern, $mu_1, $mu_2_dash) = @_;

    #print Dump { embedding_is_overlap_free => { pattern => $pattern, mu_1 => $mu_1, mu_2_dash => $mu_2_dash } }; <STDIN>;
    my $triples = {};
    foreach my $tp (@{$pattern}){
        my @parts = ();
        foreach my $x ($tp->[0], $tp->[1], $tp->[2]){
            if($x =~ m/\A\?/){
                if(exists $mu_1->{$x}){
                    push(@parts, $mu_1->{$x});
                } elsif(exists $mu_2_dash->{$x}){
                    push(@parts, $mu_2_dash->{$x});
                }
            } else {
                push (@parts, $x);
            }
        }
        my $triple = join(" ", @parts);
        #print Dump { triple => $triple, triples => $triples };
        return 0 if exists $triples->{$triple};
        $triples->{$triple} = 1;
    }
    #print "embedding_is_overlap_free\n"; <STDIN>;
    return 1;
}

sub substitute_variables_in_pattern {
    my ($p, $m_j, $m_o) = @_;
    my $new_p = [];
    foreach my $tp (@{$p}){
        my $new_tp = [];
        foreach my $x ($tp->[0], $tp->[1], $tp->[2]){
            if(exists $m_o->{$x}){
                push(@{$new_tp}, $m_o->{$x});
            } elsif(exists $m_j->{$x}){
                push(@{$new_tp}, $m_j->{$x});
            } else {
                push(@{$new_tp}, $x);
            }
        }
        push(@{$new_tp}, $tp->[3]); # add the tpID
        push (@{$new_p}, $new_tp);
    }

    #print Dump { subtitute_variables_in_pattern => { p => $p, m_j => $m_j, m_o => $m_o, new_p => $new_p }};
    return $new_p;

}

sub substitute_variables_in_mapping {
    my ($mu, $m_j, $m_o) = @_;
    my $new_mu = {};
    
    foreach my $v (keys %{$mu}){
        if(exists $m_o->{$v}){
            $new_mu->{$m_o->{$v}} = $mu->{$v};
        } elsif(exists $m_j->{$v}){
            $new_mu->{$m_j->{$v}} = $mu->{$v} ;
        } else {
            $new_mu->{$v} = $mu->{$v};
        }
    }

    #print Dump { subtitute_variables_in_mapping => { mu => $mu, m_j => $m_j, m_o => $m_o, new_mu => $new_mu }};
    return $new_mu;
}


sub AdmissibleAbstractions {
    my ($s, $p, $o) = @_;

    my @AA = ();

    foreach my $at (qw(s p o sp so po)){
        my $all_constraints_satisfied = 1;

        $all_constraints_satisfied = 0 if (
            exists $CFG->{triple_abstraction_constraints}->{s_based}->{$s} and
            exists $CFG->{triple_abstraction_constraints}->{s_based}->{$s}->{$at} and
            not $CFG->{triple_abstraction_constraints}->{s_based}->{$s}->{$at}
        );

        $all_constraints_satisfied = 0 if (
            exists $CFG->{triple_abstraction_constraints}->{p_based}->{$p} and
            exists $CFG->{triple_abstraction_constraints}->{p_based}->{$p}->{$at} and
            not $CFG->{triple_abstraction_constraints}->{p_based}->{$p}->{$at}
        );

        $all_constraints_satisfied = 0 if (
            exists $CFG->{triple_abstraction_constraints}->{o_based}->{$o} and
            exists $CFG->{triple_abstraction_constraints}->{o_based}->{$o}->{$at} and
            not $CFG->{triple_abstraction_constraints}->{o_based}->{$o}->{$at}
        );

        if($all_constraints_satisfied){
            if($at eq "s" and $CFG->{allowed_abstractions}->{s}){
                my $v = "?v0";
                if($s ne $p and $s ne $o){
                    my $tpID = &get_tpID($v, $p, $o, [$v]);
                    push(@AA, [$tpID, {$v => $s}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-1 $v $p $o\n"; <STDIN>; }
            } elsif($at eq "p" and $CFG->{allowed_abstractions}->{p}){
                my $v = "?v0";
                if($p ne $s and $p ne $o){
                    my $tpID = &get_tpID($s, $v, $o, [$v]);
                    push(@AA, [$tpID, {$v => $p}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-2 $s $v $o\n"; <STDIN>; }
            } elsif($at eq "o" and $CFG->{allowed_abstractions}->{o}){
                my $v = "?v0";
                if($o ne $s and $o ne $p){
                    my $tpID = &get_tpID($s, $p, $v, [$v]);
                    push(@AA, [$tpID, {$v => $o}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-3 $s $p $v\n"; <STDIN>; }
            } elsif($at eq "sp" and $CFG->{allowed_abstractions}->{sp}){
                my $v1 = "?v0";
                my $v2 = $s eq $p ? "?v0" : "?v1";
                my $varlist = $v1 eq $v2 ? [$v1] : [$v1, $v2];
                if($s ne $o and $p ne $o){
                    my $tpID = &get_tpID($v1, $v2, $o, $varlist);
                    push(@AA, [$tpID, {$v1 => $s, $v2 => $p}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-4 $v1 $v2 $o\n"; <STDIN>; }
            } elsif($at eq "so" and $CFG->{allowed_abstractions}->{so}){
                my $v1 = "?v0";
                my $v2 = $s eq $o ? "?v0" : "?v1";
                my $varlist = $v1 eq $v2 ? [$v1] : [$v1, $v2];
                if($s ne $p and $o ne $p){
                    my $tpID = &get_tpID($v1, $p, $v2, $varlist);
                    push(@AA, [$tpID, {$v1 => $s, $v2 => $o}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-5 $v1 $p $v2\n"; <STDIN>; }
            } elsif($at eq "po" and $CFG->{allowed_abstractions}->{po}){
                my $v1 = "?v0";
                my $v2 = $p eq $o ? "?v0" : "?v1";
                my $varlist = $v1 eq $v2 ? [$v1] : [$v1, $v2];
                if($p ne $s and $o ne $s){
                    my $tpID = &get_tpID($s, $v1, $v2, $varlist);
                    push(@AA, [$tpID, {$v1 => $p, $v2 => $o}]);
                } #else { print " s $s\n p $p\n o $o\n"; print "prevented-6 $s $v1 $v2\n"; <STDIN>; }
            }
        }

        if($CFG->{allowed_abstractions}->{spo_xxx}){
            # xxx
            if($s eq $p and $p eq $o){
                my $v1 = "?v1"; # v0
                my $varlist = [$v1];
                my $tpID = &get_tpID($v1, $v1, $v1, $varlist);
                push(@AA, [$tpID, {$v1 => $s}]);
            }
        }

        if($CFG->{allowed_abstractions}->{spo_xxy}){
            # xxy
            if($s eq $p and $p ne $o){
                my $v1 = "?v1"; # v0
                my $v2 = "?v2"; # v1
                my $varlist = [$v1, $v2];
                my $tpID = &get_tpID($v1, $v1, $v2, $varlist);
                push(@AA, [$tpID, {$v1 => $s, $v2 => $o}]);
            }
        }

        if($CFG->{allowed_abstractions}->{spo_xyx}){
            #xyx
            if($s eq $o and $s ne $p){
                my $v1 = "?v1"; # v0
                my $v2 = "?v2"; # v1
                my $varlist = [$v1, $v2];
                my $tpID = &get_tpID($v1, $v2, $v1, $varlist);
                push(@AA, [$tpID, {$v1 => $s, $v2 => $p}]);
            }
        }

        if($CFG->{allowed_abstractions}->{spo_xyy}){
            # xyy
            if($s ne $p and $p eq $o){
                my $v1 = "?v1"; # v0
                my $v2 = "?v2"; # v1
                my $varlist = [$v1, $v2];
                my $tpID = &get_tpID($v1, $v2, $v2, $varlist);
                push(@AA, [$tpID, {$v1 => $s, $v2 => $p}]);
            }
        }

        if($CFG->{allowed_abstractions}->{spo_xyz}){
            #xyz
            if($s ne $p and $s ne $o and $p ne $o){
                my $v1 = "?v1"; # v0
                my $v2 = "?v2"; # v1
                my $v3 = "?v3"; # v2
                my $varlist = [$v1, $v2, $v3];
                my $tpID = &get_tpID($v1, $v2, $v3, $varlist);
                push(@AA, [$tpID, {$v1 => $s, $v2 => $p, $v3 => $o}]);
            }
        }
    }

    # TODO's:
    # 1. allow specification of more complex constraints, e.g., sp_based, so_based, ...

    return @AA;
}

sub create_fingerprint {
    my $p = shift;
    my $tpID_count = {};
    foreach my $tp (@{$p->{pattern_parseable}}){
        $tpID_count->{$tp->[3]}++;
    }
    my @parts = ();
    foreach my $tpID (sort { $tpID_count->{$a} <=> $tpID_count->{$b} || $a <=> $b } keys %{$tpID_count}){
        push(@parts, $tpID . ":" . $tpID_count->{$tpID});
    }
    my $fingerprint = join("-", @parts) . "=" . (scalar keys %{$p->{omega}}) . "=" . $p->{number_of_mappings};
    #print Dump { create_fingerprint => { pattern_parseable => $p->{pattern_parseable}, fingerprint => $fingerprint }}; <STDIN>;
    return $fingerprint;
}

sub get_tpID {
    my ($s, $p, $o, $varlist) = @_;
    my $tpstring = join("-", $s, $p, $o);
    return $tpstring_to_tpID->{$tpstring} if exists $tpstring_to_tpID->{$tpstring};
    my $new_tpID = $last_pattern_ID++; #scalar keys %{$tpstring_to_tpID};
    $tpstring_to_tpID->{$tpstring} = $new_tpID;
    $tpID_to_tp->{$new_tpID} = { pattern_parseable => [[$s, $p, $o, $new_tpID]], omega => {}, fingerprint => "$new_tpID:1", varlist => $varlist };
    return $new_tpID;
}

sub contains_undesired_subgraph {
    my $gp = shift;
    my $gp_tpcount = {};
    
    foreach my $tp (@{$gp}){
        $gp_tpcount->{$tp->[3]}++;
    }
    
    FEUS: foreach my $usID (sort keys %{$CFG->{undesired_subgraphs_tpcount}}){
        foreach my $tpID (sort keys %{$CFG->{undesired_subgraphs_tpcount}->{$usID}}){
            next FEUS if not exists $gp_tpcount->{$tpID};
            next FEUS if $CFG->{undesired_subgraphs_tpcount}->{$usID}->{$tpID} > $gp_tpcount->{$tpID};
        }

        return 1 if &subgraph_isomorphism_exists($CFG->{undesired_subgraphs}->[$usID], $gp);
    }

    return 0;
}

sub canonicalize_triple_pattern { # TODO: take care of blank nodes
    my $tp = shift;
    my $can = [
        $tp->[0] =~ m/\A\?/ ? $tp->[0] : $tp->[0], # TODO does this make sense?
        $tp->[1] =~ m/\A\?/ ? $tp->[1] : $tp->[1], # TODO does this make sense?
        $tp->[2] =~ m/\A\?/ ? $tp->[2] : $tp->[2], # TODO does this make sense?
    ];

    if($tp->[0] =~ m/\A\?/){
        $can->[0] = "?v0";
        if($tp->[1] =~ m/\A\?/){
            if($tp->[1] eq $tp->[0]){
                $can->[1] = "?v0";
            } else {
                $can->[1] = "?v1";
            }
            if($tp->[2] =~ m/\A\?/){
                if($tp->[2] eq $tp->[0]){
                    $can->[2] = "?v0";
                } elsif($tp->[2] eq $tp->[1]){
                    $can->[2] = "?v1";
                } else {
                    $can->[2] = "?v2";
                }
            }
        } else {
            if($tp->[2] =~ m/\A\?/){
                if($tp->[2] eq $tp->[0]){
                    $can->[2] = "?v0";
                } else {
                    $can->[2] = "?v1";
                }
            }

        }
    } else {
        if($tp->[1] =~ m/A\?/){
            $can->[1] = "?v0";
            if($tp->[2] =~ m/\A\?/){
                if($tp->[2] eq $tp->[1]){
                    $can->[2] = "?v0";
                } else {
                    $can->[2] = "?v1";
                }
            }
        } else {
            if($tp->[2] =~ m/\A\?/){
                $can->[2] = "?v0";
            }
        }
    }
    #print Dump { tp => $tp, can => $can }; <STDIN>;
    return $can;
}

sub print_pattern {
    my $pattern = shift;
    print &print_pattern_to_string($pattern) . "\n";
}

sub print_pattern_to_string {
    my $p = shift;
    my $string = q{};
    foreach my $tp (@{$p->{pattern_parseable}}){
        my ($s, $p, $o, $tpID) = (
            $tp->[0] =~ m/\A\?/ ? $tp->[0] : &shorten($tp->[0]),
            $tp->[1] =~ m/\A\?/ ? $tp->[1] : &shorten($tp->[1]),
            $tp->[2] =~ m/\A\?/ ? $tp->[2] : &shorten($tp->[2]),
            $tp->[3]
        );
        #$string .= "  [" . sprintf("%0${log}d", $tpID) . "] $s $p $o .\\n";
        $string .= "  [$tpID] $s $p $o .\\n";

        
    }
    return $string;
}

sub shorten {
    my $term = shift;

    # this is a hack: how are URIs distinguished from variables and literals
    if($term =~ m/\A\?/){
        return $term;
    } elsif($term =~ m/\A</){ # OLD REGEXP: m/\A<http/
        foreach my $part (keys %{$CFG->{prefix_definitions}}){
            my $part_qm = quotemeta $part;
            if($term =~ m/\A<$part_qm(.*)>\Z/){
                #print "shorten ($term) to (" . $CFG->{prefix_definitions}->{$part} . ":" . $1 . ")\n";
                return $CFG->{prefix_definitions}->{$part} . ":" . $1;
            }
        }
        return $term;
    } elsif($term =~ m/"(.*)\^\^(<.*>)"\Z/){
        my ($value, $datatype) = ($1, $2);
        foreach my $part (keys %{$CFG->{prefix_definitions}}){
            my $part_qm = quotemeta $part;
            if($datatype =~ m/\A<$part_qm(.*)>\Z/){
                return "\"$value\"^^" . $CFG->{prefix_definitions}->{$part} . ":" . $1;
            }
    
        }
            return $term;
    } else {
        return "$term";
    }
}

sub get_timestamp {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub parse_NT_into_obj {
    my $string = shift;

          # TODO: this can be done with less code
    # Warning: this code does not handle blank nodes
    
        # URI URI URI
        if($string =~ m/(<.+>) (<.+>) (<.+>) .\n\Z/){

        my ($sID, $pID, $oID) = ($1, $2, $3);
        
        return {
                        s => {
                type     => "uri",
                rep     => $sID,
            },
                        p => {
                type     => "uri",
                rep     => $pID,
            },
                        o => {
                type     => "uri",
                rep     => $oID,
            },
                };
        }
        
    # URI URI LIT-LANG
        elsif($string =~ m/(<.+>) (<.+>) (\"(.+)\"\@(.+)) .\n\Z/){

        my ($sID, $pID, $oID) = ($1, $2, $3);
        
        return {
                        s => {
                type     => "uri",
                rep     => $sID,
            },
                        p => {
                type     => "uri",
                rep     => $pID,
            },
                        o => {
                type => "literal",
                rep => $oID,
            },
                };
        }
        
    # URI URI LIT-DAT
        elsif($string =~ m/(<.+>) (<.+>) (\"(.+)\"\^\^<(.+)>) .\n\Z/){

        my ($sID, $pID, $oID) = ($1, $2, $3);
        
        return {
                        s => {
                type     => "uri",
                rep     => $sID,
            },
                        p => {
                type     => "uri",
                rep     => $pID,
            },
                        o => {
                type     => "typed-literal",
                rep     => $oID
            },
                };
        }

    # URI URI LIT
        elsif($string =~ m/(<.+>) (<.+>) (\"(.+)\") .\n\Z/){

        my ($sID, $pID, $oID) = ($1, $2, $3);
        
        return {
                        s => {
                type     => "uri",
                rep     => $sID,
            },
                        p => {
                type     => "uri",
                rep    => $pID,
            },
                        o => {
                type     => "literal",
                rep => $oID,
            },
                };
        } else {
        print "Warning: line cannot be parsed!\n$string"; <STDIN>;
    }
}


sub config_sanity_checks {
    my $CFG = shift;

    # TODO: check number_of_workers_workers

    # TODO predicate blacklist and white list may not be used together

    if(exists $CFG->{deterministic_sampling}){

	if(not defined $CFG->{deterministic_sampling}){
		$CFG->{deterministic_sampling} = 0;
	}

	if(
		not exists $CFG->{random_number_file_1} 
		or not defined $CFG->{random_number_file_1}
		or not exists $CFG->{random_number_file_2} 
		or not defined $CFG->{random_number_file_2}
	){
		print "Error: random_number_file_1 or random_number_file_2 not provided, but required for deterministic sampling";
	} else {
		if(not -s $CFG->{random_number_file_1}){
			# try to create it
			my $count = 10_000;
			open(DAT, '>' . $CFG->{random_number_file_1}) or die "Cannot create " . $CFG->{random_number_file_1} . ": $!";
			for (1 .. $count) {
				my $num = sprintf("%.10f", rand());
				print DAT "$num\n";
				push(@{$CFG->{random_numbers_1}}, $num);
			}
			close DAT;
			print "Populated random number 1 file with $count numbers.\n";
		} else {
			# load it from file
			open(DAT,'<' . $CFG->{random_number_file_1}) or die $!;
			while(defined(my $number = <DAT>)){
				chomp $number;
				push(@{$CFG->{random_numbers_1}}, $number);
			}
			close DAT;
			print "Loaded random number file 1 with " . (scalar @{$CFG->{random_numbers_1}}) . " values.\n";
		}

		if(not -s $CFG->{random_number_file_2}){
			# try to create it
			my $count = 10_000;
			open(DAT,'>' . $CFG->{random_number_file_2}) or die "Cannot create " . $CFG->{random_number_file_2} . ": $!";
			for (1 .. $count) {
				my $num = int(rand($CFG->{number_of_workers}));
				print DAT "$num\n";
				push(@{$CFG->{random_numbers_2}}, $num);
			}
			close DAT;
			print "Populated random number 2 file with $count numbers.\n";
		} else {
			# load it from file
			open(DAT, '<' . $CFG->{random_number_file_2}) or die $!;
			while(defined(my $number = <DAT>)){
				chomp $number;
				push(@{$CFG->{random_numbers_2}}, $number);
			}
			close DAT;
			print "Loaded random number file 2 with " . (scalar @{$CFG->{random_numbers_2}}) . " values.\n";
		}

		#print Dump { random_numbers_1 => $CFG->{random_numbers_1} };
		#print Dump { random_numbers_2 => $CFG->{random_numbers_2} };
		#<STDIN>;
		#print Dump { CFG => $CFG };	
	}
    } else {
	$CFG->{deterministic_sampling} = 0;	    
	$CFG->{random_numbers_1} = [];
    }


    if(exists $CFG->{incremental_output} and $CFG->{incrememental_output} and (not exists $CFG->{incremental_output_file} or not defined $CFG->{incremental_output_file})){
        print "Error: incremental_output_file needs to be specified when incremental output is activated\n"; die;
    }

    # TODO: there should not ne more than one %d.
    if(exists $CFG->{incremental_output} and $CFG->{incrememental_output} and $CFG->{incremental_output_file} !~ m/\%d/){
        print "Error: incremental_output_file needs to contain a \"%d\" (the position of the pattern size to be inserted into the filename)\n"; die;
    }



    if(
        (not exists $CFG->{min_number_of_matched_graphs_rel} or not defined $CFG->{min_number_of_matched_graphs_rel}) and
        (not exists $CFG->{min_number_of_matched_graphs_abs} or not defined $CFG->{min_number_of_matched_graphs_abs})
    ){
        print "Error: either min_number_of_matched_graphs_rel or min_number_of_matched_graphs_abs needs to be specified.\n"; die;
    }

    if(exists $CFG->{allowed_abstractions}->{spo} and $CFG->{allowed_abstractions}->{spo}){
        print "Warning: allowing to abstact all three part of a triple can lead to a very large seach space. Only allow this abstration if you really know what you are doing. Usually, you would not want to allow this kind of abstraction. (Thus, better set allowed_abstractions -> spo to 0.)\n Hit enter to continue anyways.";  <STDIN>;
    }

    if(exists $CFG->{max_matches_per_graph} and defined $CFG->{max_matches_per_graph}){
        print "Warning: using max_matches_per_graph is usually not a good idea, because tis is not an anti-monotonic measure. Usuallz, you would not want to use this and instead use max_support_per_graph.\n Hit enter to continue anyways.";  <STDIN>;

    }

    if(exists $CFG->{min_number_of_matched_graphs_abs} and defined $CFG->{min_number_of_matched_graphs_abs} and $CFG->{min_number_of_matched_graphs_abs} > scalar @{$CFG->{input_files}}){
        print "Error: the value for the parameter min_number_of_matched_graphs_abs is larger than the number of graphs given, but should be less or equal to the number of graphs given. There cannot be results.\n"; die;
    }

    if(exists $CFG->{max_number_of_matched_graphs_abs} and defined $CFG->{max_number_of_matched_graphs_abs} and $CFG->{max_number_of_matched_graphs_abs} > scalar @{$CFG->{input_files}}){
        print "Warning: the value for the parameter max_number_of_matched_graphs is larger than the number of graphs given, but should be less or equal to the number of graphs given. This is not really a problem, the results will be the same as setting max_number_of_matched_graphs to the number of graphs given. Hit enter to continue.\n"; <STDIN>;
    }

    # TODO min/max_number_of_matched_graphs_rel between 0 and 1
    # TODO min_support should be leq max_support
    # same for min_pattern_size max_pattern_size,
    # same for min_number_of_matched_graphs max_number_of_matched_graphs
}

