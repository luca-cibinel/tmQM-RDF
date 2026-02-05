"""
This file defines helper functions that represent turtle prefixes.
Each function has the same name as the prefix it represents. Calling the function with no argument
returns the raw path of the prefix. Calling it with an argument decorates that argument with the prefix.

Prefix logic:
    - inp5: root path
    - xxYwu:
        xx: denotes the level (nm = numerical, ds = datasets, cm = complex, lg = ligand, tm = atom; first two letters after removing vowels)
        Y: denotes the first sublevel (C = centre, L = Ligand, T = TMC, A = Atom, B = Bond, S = Structure; ... first letter capitalised) (if available)
        w, u: denote, respectively, the second and third sublevel (p = properties, r = reference; ... first letter) (if available)
"""

def inp5(item = None):
	if item is None:
		 return 'resource://integreat/p5/'

	return 'inp5:' + item


def nm(item = None):
	if item is None:
		 return 'resource://integreat/p5/numerical/'

	return 'nm:' + item


def ds(item = None):
	if item is None:
		 return 'resource://integreat/p5/datasets/'

	return 'ds:' + item


def dsC(item = None):
	if item is None:
		 return 'resource://integreat/p5/datasets/complexes/'

	return 'dsC:' + item


def dsG(item = None):
	if item is None:
		 return 'resource://integreat/p5/datasets/graphs/'

	return 'dsG:' + item


def dsL(item = None):
	if item is None:
		 return 'resource://integreat/p5/datasets/ligands/'

	return 'dsL:' + item


def cm(item = None):
	if item is None:
		 return 'resource://integreat/p5/complex/'

	return 'cm:' + item


def cmT(item = None):
	if item is None:
		 return 'resource://integreat/p5/complex/TMC/'

	return 'cmT:' + item


def cmTp(item = None):
	if item is None:
		 return 'resource://integreat/p5/complex/TMC/property/'

	return 'cmTp:' + item


def lg(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/'

	return 'lg:' + item


def lgC(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/centre/'

	return 'lgC:' + item


def lgCp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/centre/property/'

	return 'lgCp:' + item


def lgCr(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/centre/reference/'

	return 'lgCr:' + item


def lgCrp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/centre/reference/property/'

	return 'lgCrp:' + item


def lgL(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/'

	return 'lgL:' + item


def lgLp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/property/'

	return 'lgLp:' + item


def lgLr(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/reference/'

	return 'lgLr:' + item


def lgLrp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/reference/property/'

	return 'lgLrp:' + item


def lgLrm(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/reference/motif/'

	return 'lgLrs:' + item


def lgLro(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/ligand/reference/occurrence/'

	return 'lgLrr:' + item


def lgB(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/bond/'

	return 'lgB:' + item


def lgBp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/bond/property/'

	return 'lgBp:' + item


def lgBr(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/bond/reference/'

	return 'lgBr:' + item


def lgBrp(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/bond/reference/property/'

	return 'lgBrp:' + item


def lgS(item = None):
	if item is None:
		 return 'resource://integreat/p5/ligand/structure/'

	return 'lgS:' + item


def tm(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/'

	return 'tm:' + item


def tmA(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/atom/'

	return 'tmA:' + item


def tmAp(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/atom/property/'

	return 'tmAp:' + item


def tmAr(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/atom/reference/'

	return 'tmAr:' + item


def tmArp(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/atom/reference/property/'

	return 'tmArp:' + item


def tmB(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/bond/'

	return 'tmB:' + item


def tmBp(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/bond/property/'

	return 'tmBp:' + item


def tmBr(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/bond/reference/'

	return 'tmBr:' + item


def tmBrp(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/bond/reference/property/'

	return 'tmBrp:' + item


def tmS(item = None):
	if item is None:
		 return 'resource://integreat/p5/atomic/structure/'

	return 'tmS:' + item


def xmls(item = None):
	if item is None:
		 return 'http://www.w3.org/2001/XMLSchema#'

	return 'xmls:' + item


def rdf(item = None):
	if item is None:
		 return 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'

	return 'rdf:' + item


def rdfs(item = None):
	if item is None:
		 return 'http://www.w3.org/2000/01/rdf-schema#'

	return 'rdfs:' + item