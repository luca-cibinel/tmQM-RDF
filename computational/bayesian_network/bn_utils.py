"""
Define a utility class that wraps the BN functionalities from pgmpy
"""

# %% Header
from tqdm import tqdm

import numpy as np
import networkx as nx
import pgmpy.models as bnmod
import pgmpy.estimators as bnest
import pgmpy.factors.discrete.CPD as bncpd

import logging
logging.basicConfig(level = logging.WARNING, force = True) # Only enables logging for messages of level WARNING or higher
                                                            # Prevents pgmpy to flood the console with logs

# %% BayesianNetwork (classes)

class BayesianNetwork(bnmod.DiscreteBayesianNetwork):
    
    def __init__(self, gml_network):
        
        self.nxgraph = nx.read_gml(gml_network)
        #self.nxgraph = nx.relabel_nodes(self.nxgraph, {str(i + 1): str(i) for i in range(len(self.nxgraph))})
        
        # Complete edge orientation (if needed)
        # non_oriented = []
        # oriented = []
        
        # for e in self.nxgraph.edges(data = True):
        #     if e[2]["mult"] > 1:
        #         non_oriented += [e[:2]]
        #     else:
        #         oriented += [e[:2]]
        
        # pdag = bnbase.PDAG(directed_ebunch = oriented, undirected_ebunch = non_oriented)
        # pdag.add_nodes_from(nx.isolates(self.nxgraph))
        
        # Call super constructor
        super().__init__(self.nxgraph)
    
    def fit(self, data, *args, **kwargs):
        """
        Fits the parameters of the BN.
        
        Arguments:
            - data: a pd.DataFrame
            - args, kwargs: further arguments to pgmpy.estimators.BayesianEstimator.get_parameters
        """
        # Standard parameter estimation via pgmpy
        estimator = bnest.BayesianEstimator(
            self, 
            data
        )
        
        theta = estimator.get_parameters(*args, **kwargs)
        self.add_cpds(*theta)
        
        # Add marginal distribution of isolated nodes
        theta_is = []
        for v in nx.isolates(self.nxgraph):
            prob = np.mean(data.loc[:, v])
            
            theta_is += [
                    bncpd.TabularCPD(v, variable_card = 2, values = [[1 - prob], [prob]], state_names = {v: [0, 1]})
                ]
            
        self.add_cpds(*theta_is)
    
    def get_complete_state_probability(self, states, log = False, unsafe = False):
        """
        Implements computation of log probabilities (based on the original pgmpy code of BayesianNetwork.get_state_probability).
        This function is implemented ONLY for complete states
        
        Arguments:
            - states: a dictionary {var: state}
            - log: boolean, should log probabilities be returned instead?
            - unsafe: boolean, should sanity checks be skipped?
        """
        if log:
            aggregate = lambda x, y: x + y
            transform = lambda x : np.log(x)
        else:
            aggregate = lambda x, y: x * y
            transform = lambda x : x
        
        # Step 1: Check that all variables and states are in the model.
        if not unsafe:
            self.check_model()
            for var, state in states.items():
                if var not in self.nodes():
                    raise ValueError(f"{var} not in the model.")
                if state not in self.states[var]:
                    raise ValueError(f"State: {state} not define for {var}")
            
            missing_vars = set(self.nodes()) - set(states.keys())
            if len(missing_vars) > 0:
                raise ValueError(f"The function get_complete_state_probability is defined only for complete cases! Missing variables: {missing_vars}")
        

        # Step 2: Compute the probability
        prob = 0 if log else 1
        for cpd in self.cpds:
            index = []
            for var in cpd.variables:
                index.append(cpd.name_to_no[var][states[var]])
            prob = aggregate(prob, transform(cpd.values[tuple(index)]))
        
        return prob
    
    def score(self, data, normalise = False, log = False, unsafe = False, progress = False):
        """
        Compute the probability of each observation within a dataset.
        
        Arguments:
            - data: a pd.DataFrame
            - normalise: should the probabilities be normalised? Ignored if log is True
            - log: boolean, should log-probabilities be computed?
            - unsafe: boolea, should sanity checks be skipped?
            - progress: boolean, should a progress bar be displayed
            
        Returns:
            - a dictionary mapping the index of each observation to its (log-)probability
        """
        
        it = data.index
        if progress:
            it = tqdm(data.index)
        
        probs = {i: None for i in data.index}
        for i in it:
            probs[i] = self.get_complete_state_probability(dict(data.loc[i,:]), log = log, unsafe = unsafe)
        
        if normalise and not log:
            Z = sum(probs.values())
            probs = {i: p/Z for i, p in probs.items()}
        
        return probs