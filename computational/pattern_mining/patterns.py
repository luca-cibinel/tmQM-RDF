"""
This file defines a series of utility classes used to read and match patterns
"""

from tqdm import tqdm
from tqdm.contrib.telegram import tqdm as tqdmtg
from collections import defaultdict
from itertools import product, combinations
from rdflib.plugins import sparql

import os
import re
import gzip
import rdflib
import pandas as pd

class progressbar:
    """
    Helper class.
    Automatically handles different logging options via tqdm.
    """
    
    def __init__(self, obj, desc = "", log = False, log_token = None, log_chat_id = None):
        """
        Arguments:
            - obj: the object to decorate with tqdm
            - desc: an optional description
            - log: whether log should be displayed
            - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
        """
        if log:
            if log_token is not None and log_chat_id is not None:
                self.obj = tqdmtg(obj, token = log_token, chat_id = log_chat_id, desc = desc)
            else:
                self.obj = tqdm(obj, desc = desc)
        else:       
            self.obj = obj
        
    def get(self):
        """
        Returns:
            - the appropriately decorated object
        """
        return self.obj

class Pattern:
    """
    This class represents a pattern, intended a set of triples of (T x V)^3 where
    T is the set of terms and V is the set of variables. It does not support Blank Nodes nor Literals.
    
    This class acts as a wrapper of rdflib.Graph, where the only allowed term types are
    rdflib.term.URIRef and rdflib.term.Variable.
    
    It supports the following operations:
        - sum (adds a list of triples to the pattern)
        - length (returns the number of triples in the pattern)
        - get item (same as rdflib.Graph)
        - contains (check if a triple is in the pattern; 
                    uses a SPARQL query to account for variable relabeling)
        - iter (same as rdflib.Graph)
        - next (same as rdflib.Graph)
        - <= (check if pattern is 'dominated by' another pattern;
              domination means that if dominating pattern does not match, neither does the dominated one; 
              uses a SPARQL query to account for variable relabeling)
        - < (check if a pattern is a subpattern of another pattern; 
             uses a SPARQL query to account for variable relabeling)
        - call (calling pattern(rdf_graph) on an rdflib.Graph or TMCGraph object will run a SPARQL query on rdf_graph
                that computes the possible matches of the pattern inside the given graph under the NRA semantics.
                The result of the call is a list of dictionaries, where every dictionary corresponds 
                to a match and maps each variable to its grounding)
        
    When building or expanding the pattern, triples must be supplied in the form of a list of lists/tuples:
        [(s1, p1, o1), (s2, p2, o2), ...]
    The terms s, p, o must either be strings ('<...>' for URIs and '?...' for variables) or rdflib.term.URIRef/rdflib.term.Variable objects.
    """
    
    def __init__(self, triples = []):
        """
        Clas constructor.
        
        Arguments:
            - triples: the list of triples in the pattern
        """
        
        self._pattern = rdflib.Graph()
        
        self._variables = set()
        self._terms = set()
        
        self._query_representation = None
        
        self._parse_triples(triples)
        
    def __add__(self, other):
        new_pattern = Pattern()
        new_pattern += self
        new_pattern += other
        
        return new_pattern
    
    def __iadd__(self, other):
        self._parse_triples(other)
        
        return self
        
    def __getitem__(self, key):
        return self._pattern[key]
    
    def __contains__(self, item):
        return Pattern([item]) <= self
    
    def __len__(self):
        return len(self._pattern)
    
    def __iter__(self):
        return iter(self._pattern)
    
    def __next__(self):
        return next(self._pattern)
    
    def __str__(self):
        return " ".join(" ".join(t.n3() for t in triple) + " ." for triple in self._pattern)
    
    def __repr__(self):
        return "\n".join(" ".join(t.n3() for t in triple) + " ." for triple in self._pattern)
    
    def __le__(self, other):
        return other._dominates(self, strict = False)
    
    def __lt__(self, other):
        return other._dominates(self, strict = True)
    
    def __call__(self, graph):
        """
        Tries to match this pattern object to the provided graph under the No Repeated Anything semantics.
        It first converts the pattern into a SPARQL query and then runs the query onto the provided graph.
        
        Arguments:
            - graph: the graph to query. Either an rdflib.Graph or a TMCGraph object
            
        Returns:
            - A list of matches, represented as dictionaries in which every rdflib.term.Variable object is
                matched to its grounding
        """
        try:
            query_result = graph.query(self._as_sparql_query())
        except Exception as e:
            print("Attempted to parse the following query:")
            print(self._as_sparql_query())
            
            raise e
            
        return [{v: row[v] for v in self._variables} for row in query_result]
    
    def _as_sparql_query(self, return_string = False, assume_named_graphs = False):
        """
        Converts this pattern object into a viable SPARQL query that will try to match the pattern
        under the No Repeated Anything semantics.
        
        This method tries to recover the precomputed query representation, if available, before computing a new one.
        Recomputation is necessary if new triples are added to the pattern (see _parse_triples).
        
        Arguments:
            - return_string: boolean; whether the string representation should also be returned. Default = False
            - assume_named_graphs: boolean; whether the query should encompass the pattern in a GRAPH statement (assuming
                                    that the pattern will be used in a graph dataset). Default = False
        
        Returns:
            - The prepared query representation
            - The string representation of the query (only if return_string is True)
        """
        
        where_clause = "\n".join([" ".join([t.n3() for t in triple]) + " ." for triple in self._pattern])
        
        # Prescribe that no two variables can match to the same term
        filter_clause = [f"\t{v1.n3()} != {v2.n3()}" for v1, v2 in combinations(self._variables, 2)]
        
        # Prescribe that no variable can match to a term already present in the pattern
        filter_clause += [f"\t{v.n3()} != {t.n3()}" for v, t in product(self._variables, self._terms)]
        
        # Join the filters together
        filter_clause = " &&\n".join(filter_clause)
        
        if assume_named_graphs:
            self._query_string = "SELECT *\nWHERE {\nGRAPH ?graph {" + where_clause + "}\nFILTER (\n" + filter_clause + "\n)\n}"
        else:
            self._query_string = "SELECT *\nWHERE {\n" + where_clause + "\nFILTER (\n" + filter_clause + "\n)\n}"
        
        self._query_representation = sparql.prepareQuery(self._query_string)
        
        if return_string:
            return self._query_representation, self._query_string
        
        return self._query_representation
    
    def _parse_triples(self, triples):
        """
        Parses a set of triples that needs to be added to the pattern. 'Parsing' means checking that the new triples
        are in the correct format (see _sanitise_triple), adding them to the internal rdflib.Graph representation and updating the set of terms and
        variables used by the pattern.
        
        This method also resets the cached query representation, so as to force the computation of the correct representation in
        _as_sparql_query.
        
        Arguments:
            - triples: the list of triples to be added    
        """
        
        if len(triples) > 0:
            # Set of triples is being updated, query representation will have to be recomputed.
            self._query_representation = None
        
        for triple in triples:
            # Sanitise the triple (check format and cast to rdflib.term.URIRef/rdflib.term.Variable objects)
            sanitised = type(self)._sanitise_triple(triple)
            
            # Save the triple
            self._pattern.add(sanitised)
            
            # Extract the variables and the terms added by the triple
            self._variables = self._variables.union(v for v in sanitised if isinstance(v, rdflib.term.Variable))
            self._terms = self._terms.union(t for t in sanitised if isinstance(t, rdflib.term.URIRef))
    
    def _dominates(self, other, strict = False):
        """
        Verifies whether this pattern object dominates another pattern, while accounting for the possible
        relabeling of the variables, in the sense that if the subpattern has no matches in a graph, neither has the superpattern.
        This is checked by trying to match this pattern to the rdflib.Graph representation of the candidate superpattern.
        A strict comparison (checking whether the subpattern relation also holds) requires checking 
        whether there exists a match where all the variables of this pattern are matched to variables in the candidate superpattern.
        
        Arguments:
            - other: the Pattern object that represents the candidate superpattern
            - strict: boolean, whether the comparison should also require that variables are matched to variables
        
        Returns:
            - whether this pattern object dominates the other
        """
        
        if not isinstance(other, Pattern):
            raise TypeError("Cannot compare a pattern to an object of class " + str(type(other)) + " !")
            
        matches = self(other._pattern)
        
        if strict:
            for m in matches:
                if all(isinstance(grounding, rdflib.term.Variable) for grounding in m.values()):
                    # Found a match in which all the groundings are variables, return True
                    return True
        
            # Couldn't find a match in which all groundings are variables, return False
            return False
        
        return len(matches) > 0
    
    @staticmethod
    def _sanitise_triple(triple):
        """
        Sanitises a candidate triple: checks that the format is correct (a list/tuple of three elements which are either strings
        or rdflib.term.URIRef/rdflib.term.Variable objects) and converts all the elements to rdflib.term.URIRef/rdflib.term.Variable objects.
          
        If the elements of the triple are provided as strings, "<...>" identifies URIs and "?..." identifies variables.
        
        Arguments:
            - triple: the triple to be sanitised
        
        Returns:
            - the sanitised triple (as a tuple)
        """
        if len(triple) != 3:
            raise Exception("Error! Triples must be provided in a [(subject, predicate, object)] format!")
            
        s_triple = []
        
        for term in triple:
            if not (isinstance(term, str) or isinstance(term, rdflib.term.URIRef) or isinstance(term, rdflib.term.Variable)):
                raise Exception("Error! Triples must be sequences of strings or rdflib.term.URIRef/Variable objects!")
            
            if not isinstance(term, rdflib.term.Identifier): # This is a string
                if term.startswith("<"): # Found a URI: convert to rdflib.term.URIRef
                    term = rdflib.term.URIRef(term[1:-1])
                elif term.startswith("?"): # Found a variable: convert to rdflib.term.Variable
                    term = rdflib.term.Variable(term)
                else:
                    raise Exception(f"Unrecognised term: {term}. Terms must be URIs (<...>) or variables (?...)!")
                
            s_triple += [term]
        
        return tuple(s_triple)
    
    def query(self, query_object, *args, **kwargs):
        """
        Wrapper for the query method of the internal rdflib.Graph pattern representation. See the documentation
        for rdflib.Graph.query.
        
        Arguments:
            - query_object: the query object
            - *args: further arguments to pass to rdflib.Graph.query
            - **kwargs: further arguments to pass to rdflib.Graph.query
        
        Returns:
            - an rdflib.query.Result object
        """
        return self._pattern.query(query_object, *args, **kwargs)
    
class PatternSet:
    """
    This class reads the .dat.gz file for a specific pattern size, within a specified directory, and stores each pattern as a
    Pattern object. It acts as a wrapper for a dictionary where each Pattern object is stored under its pattern id.
    
    This class supports the following operations:
        - get item (access Pattern objects via pattern id)
        - contains (as in dictionary, checks if a pattern id is stored)
        - length (return the number of patterns)
        - iter (as in dictionary)
        - next (as in dictionary)
        - keys (as in dictionary)
        - values (as in dictionary)
        - items (as in dictionary)
        - call (given a list of rdflib.Graph/TMCGraph objects, returns a pandas.DataFrame where the columns are indexed by the pattern ids and the
                rows correspond to the graphs in the list. Each entry contains the number of matches of the corresponding pattern in the
                corresponding graph)
    """
    
    def __init__(self, directory, size, output_file_pattern = "output-size-%d.dat", log = False, log_token = None, log_chat_id = None):
        """
        Class constructor.
        
        Arguments:
            - directory: the directory from which to retrieve the patterns.
            - size: the size of the patterns to retrieve.
            - output_file_pattern: the pattern of the name of the .dat file to retrieve. 
                It must contain a '%d' token in place of the pattern size.
            - log: optional boolean parameter. Whether a detailed log of the retrieving operations should be provided.
            - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
        """
        #self._multiproc_manager = Manager() # multiprocessing.Manager object (used to manage dictionary access)
        
        self._directory = directory
        self._size = size
        
        if log:
            print(f"Reading .dat file for size {self._size}...")
        
        with gzip.open(self._directory + (output_file_pattern % self._size) + ".gz", "rt") as f:
            pattern_lines = f.read().split('\n')
        
        if log:
            print("Done!")
        
        """
        In the .dat files, lines start with a number, with each number having the following meaning:
            - 0: beginning of pattern description
            - 1: number of matches
            - 6: pattern description
            - 7: beginning of list of match description (containes matched graph name)
            - 8: variable grounding within matched graph
            - 9: end of match description (if multiple matches are available for the same graph
                                           more 8-type lines will follow a 9-type line)
        """
        
        pattern_lines = progressbar(
                pattern_lines,
                f"Processing pattern data (size: {self._size})",
                log,
                log_token,
                log_chat_id
            )
            
        self._patterns = defaultdict(Pattern)
        for line in pattern_lines.get():
            if line.startswith('0'):
                pattern_id = line[2:].strip()
        
            if line.startswith('6'):
                self._patterns[pattern_id] += [re.match(r"6 (.*)\t(.*)\t(.*)\t\d*", line).groups()]
                
        self._patterns = dict(self._patterns)
        self._keys = list(self._patterns.keys())
        
    def __getitem__(self, key):
        try:
            return self._patterns[key]
        except KeyError:
            return self._patterns[self._keys[key]]
    
    def __contains__(self, item):
        return item in self._patterns
    
    def __len__(self):
        return len(self._patterns)
    
    def __iter__(self):
        return iter(self._patterns)
    
    def __next__(self):
        return next(self._patterns)
    
    def keys(self):
        return self._patterns.keys()
    
    def values(self):
        return self._patterns.values()
    
    def items(self):
        return self._patterns.items()
    
    # def __call__(self, rdf_graphs, inactive_graphs = [], key = None, desc = None, log = False, log_token = None, log_chat_id = None, n_workers = 1):
    #     """
    #     Computes a feature matrix containing the number of matches of the stored patterns in each provided graph.
        
    #     Arguments:
    #         - rdf_graphs: a list of rdflib.Graph/TMCGraph objects
    #         - inactive_graphs: optional, a list of (pattern_id, graphs) for which it is known that the given pattern does not match the given graph.
    #             The possible matches will not be computed and a 0 will be automatically placed in the feature matrix.
    #         - key: an optional function that assigns a short identifier to each RDF graph
    #         - desc: an optional description for the progress bar
    #         - log: optional boolean parameter. Whether a detailed log should be provided.
    #         - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
    #         - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            
    #     Returns:
    #         - A pandas.DataFrame containing the described feature matrix
    #     """
        
    #     multiproc_manager = Manager()
        
    #     def partition(l, n):
    #         """
    #         Internal utility function. Partitions a list l into n sublists. It tries to allocate elements to the sublists
    #         as evenly as possible, assigning leftover elements to the first sublist.
    #         If n < len(l), then it just partitions the list into singletons.
            
    #         Arguments:
    #             - l: the list to partition
    #             - n: the desired number of sublists
            
    #         Returns:
    #             - a list of sublists
    #         """
    #         if len(l) < n:
    #             return [[i] for i in l]
            
    #         chunk_size = len(l)//n
    #         remainder = len(l) % n
            
    #         parts = [l[:(chunk_size + remainder)]]
            
    #         for i in range(1,n):
    #             parts += [ l[ (chunk_size*i + remainder):(chunk_size*(i+1) + remainder) ] ]
            
    #         return parts
        
    #     desc = f"Matching patterns of size {self._size}..." if desc is None else desc
        
    #     raw_df_data = multiproc_manager.dict()
        
    #     workers = []
        
    #     for i, chunk in enumerate(partition(list(self.items()), n_workers)):
    #         workers += [
    #                 Process(
    #                         target = type(self)._count_pattern_matches,
    #                         args = (chunk, rdf_graphs, raw_df_data, inactive_graphs, i, desc, log, log_token, log_chat_id)
    #                     )
    #             ]
            
    #     for w in workers:
    #         w.start()
            
    #     for w in workers:
    #         w.join()
        
    #     if key is not None:
    #         return pd.DataFrame(dict(raw_df_data), index = [key(g) for g in rdf_graphs])
        
    #     return pd.DataFrame(dict(raw_df_data))
    
    # @staticmethod
    # def _count_pattern_matches(pattern_items, rdf_graphs, res, inactive_graphs, worker_id, desc, log = False, log_token = None, log_chat_id = None):
    #     """
    #     Utility method used during parallelised matches computation.
        
    #     Computes a feature matrix containing the number of matches of the stored patterns in each provided graph.
        
    #     Arguments:
    #         - pattern_items: a list of (pattern_id, pattern) pairs
    #         - rdf_graphs: a list of rdflib.Graph/TMCGraph objects
    #         - res: a multiprocessing.managers.DictProxy which will contain the resulting match counts
    #         - inactive_graphs: optional, a list of (pattern_id, graphs) for which it is known that the given pattern does not match the given graph.
    #             The possible matches will not be computed and a 0 will be automatically placed in the feature matrix.
    #         - worker_id: the id of the current worker
    #         - desc: an optional description for the progress bar
    #         - log: optional boolean parameter. Whether a detailed log should be provided.
    #         - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
    #         - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            
    #     """
    #     for pattern_id, pattern in progressbar(pattern_items, f"W{worker_id} | {desc}", log, log_token, log_chat_id).get():
    #         res[pattern_id] = [ (len(pattern(graph)) if (pattern_id, graph) not in inactive_graphs else 0) for graph in rdf_graphs ]
    
class PatternBatch:
    """
    This class reads an entire folder and stores all the patterns (of every size) into appropriate PatternSet objects.
    It acts as a wrapper for a dictionary which stores each PatternSet object under its corresponding size.
    
    It supports the following operations:
        - get item (access the PatternSet object via the pattern size)
        - contains (same as dictionary)
        - length (returns the number of sizes found)
        - iter (same as dictionary)
        - next (same as dictionary)
        - keys (same as dictionary)
        - values (same as dictioanry)
        - items (same as dictionary)
        - call (given a list of rdflib.Graph/TMCGraph objects, returns a pandas.DataFrame where the columns are indexed by the pattern ids of all the
                patterns in the retrieved PatternSet objects and the rows correspond to the graphs in the list. 
                Each entry contains the number of matches of the corresponding pattern in the corresponding graph)
    """
    
    def __init__(
                self, 
                directory, 
                detect_domination = False, 
                output_file_pattern = "output-size-%d.dat", 
                log = False,
                log_token = None,
                log_chat_id = None
            ):
        """
        Class constructor.
        
        Arguments:
            - directory: the directory from which to retrieve the patterns.
            - detect_domination: boolean, whether dominating patterns should be identified
                (see Pattern._dominates)
            - output_file_pattern: the pattern of the name of the .dat file to retrieve. 
                It must contain a '%d' token in place of the pattern size.
            - log: optional boolean parameter. Whether a detailed log of the retrieving operations should be provided.
            - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
        """
        
        self._directory = directory
        
        def partition(l, n):
            """
            Internal utility function. Partitions a list l into n sublists. It tries to allocate elements to the sublists
            as evenly as possible, assigning leftover elements to the first sublist.
            If n < len(l), then it just partitions the list into singletons.
            
            Arguments:
                - l: the list to partition
                - n: the desired number of sublists
            
            Returns:
                - a list of sublists
            """
            if len(l) < n:
                return [[i] for i in l]
            
            chunk_size = len(l)//n
            remainder = len(l) % n
            
            parts = [l[:(chunk_size + remainder)]]
            
            for i in range(1,n):
                parts += [ l[ (chunk_size*i + remainder):(chunk_size*(i+1) + remainder) ] ]
            
            return parts
        
        # Retrieve the available pattern sizes by scanning for all the available .dat files using a regexp
        regexp = output_file_pattern.replace("%d", "(\\d+)") + ".gz$"
        sizes = [int(m.group(1)) for m in [re.match(regexp, f) for f in os.listdir(self._directory)] if m is not None]
        self._sizes = sorted(sizes)
        
        # Start extracting pattern sets
        self._pattern_sets = dict()
        
        self._read_pattern_sets(output_file_pattern, log, log_token, log_chat_id)
            
    def __getitem__(self, key):
        return self._pattern_sets[key]
    
    def __contains__(self, item):
        return item in self._pattern_sets
    
    def __len__(self):
        return len(self._pattern_sets)
    
    def __iter__(self):
        return iter(self._pattern_sets)
    
    def __next__(self):
        return next(self._pattern_sets)
    
    def keys(self):
        return self._pattern_sets.keys()
    
    def values(self):
        return self._pattern_sets.values()
    
    def items(self):
        return self._pattern_sets.items()
    
    def _read_pattern_sets(self, output_file_pattern, log, log_token, log_chat_id):
        """
        Utility function used during object initialisation.
        Given a set of pattern sizes, reads the pattern set associated with that size.
        
        Arguments:
            - sizes_to_process: set of sizes to process
            - directory: the directory from which to retrieve the patterns.
            - output_file_pattern: the pattern of the name of the .dat file to retrieve. 
                It must contain a '%d' token in place of the pattern size.
            - res: a multiprocessing.managers.DictProxy which will contain the resulting pattern sets
            - log: boolean, whether detailed log should be provided
            - log_token: optional, the token of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
            - log_chat_id: optional, the chat id of the telegram bot to which tqdm updates should be sent (see https://tqdm.github.io/docs/contrib.telegram/)
        """
        
        for size in self._sizes:
            self._pattern_sets[size] = PatternSet(self._directory, size, output_file_pattern, log, log_token, log_chat_id)
    
if __name__ == "__main__":
    log_token = "7919012232:AAHvPKh8hYBSB0TsD30T7nvM5eG4LVM7WQk"
    log_chat_id = "452834227"
    
    print("----------------- DEBUGGING...")
    
    import sys
    sys.path.append("../../data/utils/")
    sys.path.append("../../data/preprocessing/")
    from TMCGraph import TMCGraph
    TMCGraph.path_to_tmQM_RDF = "../../data/tmQM-RDF/"
    
    target_directory = "./../results/for_paper_1/ligand_behaviour_intensive/"
    
    dataset = [TMCGraph("ABEVAH.ttl"), TMCGraph("WAJJOH.ttl")]
    
    print("\n> Pattern\n")
    
    p1 = Pattern()
    p2 = Pattern()
    
    p1 += [("<a>", "<b>", "?v1")]
    p2 += [("<a>", "<b>", "?v2"), ("?v1", "<c>", "<d>")]
    
    print("p1:", p1)
    print("p2:", p2)
    print("Pattern():", Pattern())
    
    print("p2 <= p1:", p2 <= p1, "(expected: T)")
    
    p3 = Pattern([("<a>", "<b>", "?v1"), ("?v1", "<c>", "<d>")])
    p4 = Pattern([("<a>", "<b>", "?v2"), ("?v2", "<c>", "<d>"), ("?v2", "<e>", "<f>")])
    p4_1 = Pattern([("<a>", "<b>", "<g>"), ("<g>", "<c>", "<d>"), ("<g>", "<e>", "<f>")])
    p5 = Pattern([("<a>", "<b>", "?v3"), ("<d>", "<c>", "?v3"), ("?v3", "<e>", "<f>")])
    
    print("p3:", p3)
    print("p4:", p4)
    print("p4_1:", p4_1)
    print("p5:", p5)
    print("p4 <= p3:", p4 <= p3, "(expected: T)")
    print("p4_1 <= p3:", p4_1 <= p3, "(expected: T)")
    print("p4_1 < p3:", p4_1 < p3, "(expected: F)")
    print("p5 <= p3:", p5 <= p3, "(expected: F)")
    
    
    print("\n> PatternSet\n")
    
    ps = PatternSet(target_directory, 2, log = True)
    print("\n\n")
    X = ps(dataset)
    
    for i in range(len(dataset)):
        for j in range(len(ps)):
            print(X.iloc[i, j], end = " ")
        print("")
        
    print("\n Expected:\n 1 4 1 1 0\n 0 6 1 1 0")
    
    
    print("\n> PatternBatch\n")
    
    pb = PatternBatch(target_directory, detect_domination = True, log = True, 
                      output_file_pattern = "TEMP_output-size-%d.dat",
                      log_token = log_token, log_chat_id = log_chat_id,
                      n_workers = 4)
    
    old_pattern_sets = pb._pattern_sets
    pb._pattern_sets = {s: p for s, p in pb.items() if s <= 3}
    
    Y = pb(dataset, log = True,
                      log_token = log_token, log_chat_id = log_chat_id, 
                      n_workers = 4)
    
    print("")
    print("")
    for i in range(len(dataset)):
        for j in range(Y.shape[1]):
            print(Y.iloc[i, j], end = " ")
        print("")
        
    print("\n Expected:\n 0 2 0 0 4 4 4 2 0 1 1 4 1 1 0\n 0 0 0 0 0 6 6 6 0 0 0 6 1 1 0 \n ...")