import numpy as np
import matplotlib.pyplot as plt
from math import factorial 
from ase import Atoms
from ase.build import make_supercell
from ase.spacegroup import crystal
from ase.io.vasp import write_vasp
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs, 
                                             generate_sqs_from_supercells, 
                                             generate_sqs_by_enumeration, 
                                             generate_target_structure)
from icet.input_output.logging_tools import set_log_config
set_log_config(level='INFO')

## FUNCTION DEFINITIONS ##
# try to write docstrings for each of these functions (if i can be bothered)

# want to be able to compare resultant SQS cluster vector with cluster vector of a perfectly random alloy with same composition & cluster space
# formula for ensemble-averaged correlation functions in perfectly random case is <P_(k,m)>_R = (2*x-1)^k, so only dependent on comp. x & cluster order k
def rand_corr_func(comp: float, order: int):
    return (2*comp-1)**order


# also want to compare to the analytical expression for binomial probability mass function
# this gives normalised measure of prob of getting each no. Br atoms 
def binomial(succ: np.array[int], trials: int, prob: float):
    arr = np.zeros(len(succ))
    for i in range(len(succ)):
        nCk = factorial(trials)/(factorial(succ[i])*factorial(trials-succ[i]))
        arr[i] = nCk*(prob**succ[i])*(1-prob)**(trials-succ[i])
        #print(arr[i])
        
    return arr

# sorts atoms according to desired order of chemical symbols
def atom_sort(atoms: Atoms, sym_order: dict[str, int]):
    """
    Sorts the chemical symbol labels of an `Atoms` object, using order specified by `sym_order` dictionary.
    
    `sym_order` should contain all chemical symbols in structure of interest, with item values corresponding to indices in a sorted list.
    
    E.g. let symbol `'Cs'` have index `0` to order it first in the output file.
    """
    tags = atoms.get_chemical_symbols()
    sorted_tags = sorted([(tag, i) for i, tag in enumerate(tags)], key = lambda symbol : sym_order[symbol[0]])
    sorted_indices = [i for tag, i in sorted_tags]
    return atoms[sorted_indices]

# some function that searches through every Pb atom, then finds species of the 6 nearest neighbours and counts no. Br atoms
# in order to find distribution of octahedra with every possible no. Br atoms
def octahedra_distribution(atoms: Atoms):

    # need to find index of the first and last Pb atoms to loop only through these
    # always going to have same no. Cs & Pb atoms, index of the first Pb atom tells us no. Cs t.f. Pb index = no. Pb
    for i in range(len(atoms)):
        if list(atoms.symbols)[i] == 'Pb':
            N_Pb = i
            break
    #print('N_Pb =', N_Pb)

    # know that there will be Pb atoms from indices [N_Pb, N_Pb*2], with halide atoms filling out the indices from N_Pb*2 to the end
    # so can provide list of all halide atoms to calculate distances from each Pb to each halide 
    idx_halides = list(range(N_Pb*2, len(atoms)))
    
    # each octahedron contains central Pb atom and 6 surrounding halide atoms
    # for each Pb atom (symbols in range [N_Pb, N_Pb+N_Pb]), find the 6 nearest neighbour atoms (for cubic will all be equidistant bc of cell symmetry)
    # by searching within atom position list for the current Pb atom, and then returning the indices of the 6 closest atoms
    # can then use these indices to find corresponding chemical symbol (either I or Br) and count how many Br  
    N_Br = np.zeros(N_Pb, dtype=int)
    for i in range(N_Pb, N_Pb*2):
        #pos_Pb = sqs_cubic.positions[i]
        #print(pos_Pb)
        dists = atoms.get_distances(i, idx_halides, mic=True)
        
        # make a list of tuples, pairing up the distance with halide index
        # need to sort dists by min dist but keep track of the halide indices too
        # sorted() sorts by first element of tuple by default
        dists_t = list(zip(dists, idx_halides))
        dists_t = sorted(dists_t)[0:6]
        halides_sorted = [j[1] for j in dists_t]
        #print(halides_sorted)
        
        # finally, get chemical symbol for each index of the 6 closest atoms to each Pb
        syms = [0]*6
        count_Br = 0
        for j in range(6):
            syms[j] = atoms.symbols[halides_sorted[j]]
            
            # count up no. Br symbols for each octahedron and store in array
            # only need info on the Br occupations s(Br) as s(I) = 1-s(Br)
            if syms[j] == 'Br':
                count_Br = count_Br + 1    
                
        N_Br[i-N_Pb] = count_Br        
        #print('Octahedron no.', i-N_Pb+1, ':', syms, ', no. Br atoms :', count)
        
    #print('No. Br atoms in each octahedra:', N_Br)
    #print(np.array(np.unique(N_Br, return_counts=True)).T)
    
    # outputs the unique values of no. Br in octahedra and their occurrence/counts within the full list
    return (np.array(np.unique(N_Br, return_counts=True))[:][0], np.array(np.unique(N_Br, return_counts=True))[:][1])


# actual SQS generation and evaluation of accuracy compared to ideal random alloy
def generate_sqs(cs: ClusterSpace, x: float, supercell: Atoms, N_steps: int):

    # target concentrations: should correspond to final SQS compositions A_1-x, B_x
    tc = {'A': {'I': 1-x, 'Br': x}}
    
    # order of elements for CsPb(I_xBr_1-x)3
    # if using different material, change this!
    CsPbIBr_order = {
        'Cs': 0,
        'Pb': 1,
        'I': 2,
        'Br': 3,
        }

    # this is the bad way to do it
    #sqs = generate_sqs(cluster_space=cs, max_size=8, target_concentrations=tc, include_smaller_cells=False)
    # this doesn't help as the way in which atoms are ordered from atoms.repeat is already in cell-major format,
    # applying matrix transformation w atom-major order doesn't change that how the data is structured
    #sqs = make_supercell(sqs,[[1,0,0],[0,1,0],[0,0,1]],order="atom-major")

    # generate SQS from the specified (atom-major) supercells, will give cell-major output so then sort according to desired order of ch. symbols
    sqs = generate_sqs_from_supercells(cluster_space=cs, supercells=supercell, target_concentrations=tc, n_steps=N_steps)
    sqs = atom_sort(sqs, CsPbIBr_order)
    write_vasp("sqs.vasp",sqs,direct=True)
    #view(sqs_cubic)

    # cluster vector; correlation functions of each cluster in sqs 
    cv_sqs = cs.get_cluster_vector(sqs)
    print('Cluster vector of generated SQS:', np.around(cv_sqs,4))

    # will need to use info from the cluster space to gen the perfectly random corr. funcs.
    # accessing cluster 'order' (k) from each dictionary corresponding to elements of cs
    # find indices of each new k, to find how many clusters of each order are in cs
    # also need radii of each cluster to find radial distribution of correlations/errors, get these in same loop
    cv_radii = np.empty(len(cs))
    cv_rand = np.empty(len(cs))
    for i in range(len(cs)):
        cv_radii[i] = cs.as_list[i]['radius']
        k = cs.as_list[i]['order']
        cv_rand[i] = rand_corr_func(x,k) 
        
        if k == 1:
            k2_idx = i + 1
        
        if k == 2:
            k3_idx = i + 1
    
        if k == 3:
            k4_idx = i+1
            
    #print('Cluster vector of perfectly random alloy:', np.around(cv_rand,4))

    # absolute percentage error in the sqs correlation functions compared to perfect random alloy
    cv_err = np.empty(len(cs))
    cv_err = np.abs((np.abs(cv_rand)-np.abs(cv_sqs))/cv_rand)
    
    # should get residuals as well, more useful than % errors
    cv_residuals = np.abs(cv_sqs) - np.abs(cv_rand)
    #print('Residuals:', cv_residuals)

    # RMSE of the SQS correlations
    cv_rmse = np.sqrt((cv_residuals**2).mean())
    print('% Error in SQS clusters:', np.around(cv_err,3))
    print('RMSE of SQS cluster vector:', round(cv_rmse,3))
    
    # also find how many correlations match "perfectly" (below small tolerance)
    count_match = 0
    tol_match = 10**(-6)
    for i in range(len(cs)):
        if np.abs(cv_residuals[i]) < tol_match:
            count_match = count_match + 1
        
    print('No. exactly matching SQS correlations:', count_match)
    
    # plot the correlations of the SQS cluster vector, compared with ideal random alloy
    # don't include the zerolet or singlet correlations (always = 1, 0.333), only one data pt for each
    # distinguish 2nd from 3rd order correlations
    plt.figure(1)
    #plt.plot(np.abs(cv_rand[1]), cv_sqs[1], 's', color='tab:blue', label="k=1")
    plt.plot(np.abs(cv_rand[k2_idx:k3_idx-1]), cv_sqs[k2_idx:k3_idx-1], 'o', color='tab:blue', label="k=2")
    plt.plot(np.abs(cv_rand[k3_idx:k4_idx-1]), cv_sqs[k3_idx:k4_idx-1], '^', color='tab:orange', label="k=3")
    #plt.plot(np.abs(cv_rand[k4_idx:]), cv_sqs[k4_idx:], 'v', color='tab:green', label="k=4")
    plt.plot(np.abs(cv_rand[k2_idx:]), np.abs(cv_rand[k2_idx:]), '--', color='k', label="Error=0")
    plt.xticks([np.abs(cv_rand[k2_idx]), np.abs(cv_rand[k3_idx])])
    #plt.xticks([np.abs(cv_rand[2]), np.abs(cv_rand[k3_idx]), np.abs(cv_rand[k4_idx])])
    plt.xlabel("Ideal random alloy correlations")
    plt.ylabel("SQS correlations")
    plt.legend()
    plt.savefig("sqs_correlations_ideal.png")
    
    # plot the residuals of the SQS correlations vs cluster radius
    # shows the preferential minimisation of short-range correlation errors
    # distinguish 2nd from 3rd order correlations
    plt.figure(2)
    #plt.plot(cv_radii[1], cv_residuals[1], 's', color='tab:blue', label="k=1")
    plt.plot(cv_radii[k2_idx:k3_idx-1], cv_residuals[k2_idx:k3_idx-1], 'o', color='tab:blue', label="k=2")
    plt.plot(cv_radii[k3_idx:k4_idx-1], cv_residuals[k3_idx:k4_idx-1], '^', color='tab:orange', label="k=3")
    #plt.plot(cv_radii[k4_idx:], cv_residuals[k4_idx:], 'v', color='tab:green', label="k=4")
    plt.plot([cv_radii[k2_idx], np.max(cv_radii)], [0, 0], '--', color='k', label="Error=0")
    #plt.xticks([2, np.max(cv_radii)+0.5])
    #plt.yticks(cv_sqs)
    plt.xlabel("Cluster radius (Å)")
    plt.ylabel("SQS residuals")
    plt.legend()
    plt.savefig("sqs_correlation_residuals.png")

    # some function that searches through every Pb atom, then finds species of the 6 nearest neighbours and counts no. Br atoms
    # in order to find distribution of octahedra with every possible no. Br atoms 
    oct_distr = octahedra_distribution(sqs)
    occ_sqs = oct_distr[0]
    counts_sqs = oct_distr[1]
    print('Br per octahedra & counts for SQS:', oct_distr)

    # normalising by total octahedra count (e.g. 27 for cubic) to convert count -> prob
    prob_sqs = counts_sqs/(np.sum(counts_sqs))

    # find RMSE in the difference between SQS distribution and binomial distribution!
    # gives me another measure of error to compare with that of the correlation functions
    N_halides = 6
    prob_residuals = prob_sqs - binomial(occ_sqs, N_halides, x)
    prob_rmse = np.sqrt(((prob_residuals)**2).mean())
    print('RMSE of SQS prob distribution:', round(prob_rmse,3))
    
    # all possible values of Br occuptation in each octahedron, 0->6
    k_vals = np.arange(N_halides+1)

    # plot frequency of the Br counts across all the octahedra,
    # compare the SQS octahedral Br distribution with expected binomial distribution in random alloy
    plt.figure(3)
    plt.plot(occ_sqs, prob_sqs, 'o', label = "SQS")
    plt.plot(k_vals, binomial(k_vals, N_halides, x), '-', label = "Binomial")
    plt.xticks(k_vals)
    plt.xlabel("No. Br atoms in octahedron")
    plt.ylabel("Probability")
    plt.legend()
    plt.savefig("sqs_distributions.png")
    
    plt.figure(4)
    plt.plot(occ_sqs, prob_residuals, 'o', color='tab:blue', label = "Residuals")
    plt.plot(occ_sqs, np.zeros(len(occ_sqs)), '--', color='k', label = "Error=0")
    plt.xticks(occ_sqs)
    plt.xlabel("No. Br atoms in octahedron")
    plt.ylabel("SQS residuals")
    plt.legend()
    plt.savefig("sqs_distribution_residuals.png")
    