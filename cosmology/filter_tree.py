import numpy as np
import yt
import ytree
import sys




def select_root_halos(a, criteria, fields):
    """
    Select desired halos given a set of criteria

    Parameters
    -----------
    a        : ytree object
    criteria : sequence of conditional strings that can be
               applied to the tree nodes
    fields   : fields used in the criteria

    Returns
    --------
    roots    : List of root halos in tree node that satisfy the
               criteria
    """
    selection = a.select_halos(criteria, fields = fields)
    roots = [x.find_root() for x in selection]

    return roots


def get_root_halos(a, criteria, fields,
                      unique = True,
                      mutually_exclusive = ''):
    """
    Loop over a list of criteria and return root halos
    that satisfy each criteria.

    Parameters
    -----------
    a        : ytree object
    criteria : sequence of conditional strings that can be
               applied to the tree nodes
    fields   : fields used in the criteria

    unique   : bool, optional
               ensure that the root nodes within each criteria selection
               are unique. Default : True

    mutually_exlcusive : bool, optional
        Method to ensure that halos across criteria are mutually exclusive.
        By default this is '', and this is not used. If set to 'sequential',
        this is meant to be used with criteria that are (in effect) nested
        such that the first is the most restrictive and the last the most
        general, such that it is likely that all roots in each subsequent
        selection are also included in the previous. Default : 'sequential'

    Returns
    ---------
    final_roots : list
        list of the list of root nodes that satisfy the
        criteria
    """

    num_criteria = len(criteria)

    all_roots = []
    for c in criteria:
        all_roots.append(select_root_halos(a, c, fields))

    # now make sure these are unique

    final_roots = [None]*num_criteria
    if mutually_exclusive == 'sequential':

        final_roots[0] = all_roots[0] # always in this mode

        for i in np.arange(1, num_criteria):
            final_roots[i] = [all_roots[i][ii] for ii in np.arange(np.size(all_roots[i])) if (not(all_roots[i][ii] in all_roots[i-1]))]

    elif mutually_exclusive == '':
        final_roots = all_roots
    else:
        print("Mutually exclusive selection type " + mutually_exclusive + " not yet supported")
        raise(ValueError)

    if unique:
        #
        # make sure roots halos are unique within each list
        # np.unique does not work on TreeNode objects (unfortunately)
        # so doing this in a really dumb way
        #
        final_roots_ids = [None]*num_criteria
        for i in np.arange(num_criteria):
            final_roots_ids[i] = [x['id'] for x in final_roots[i]]

        final_roots = [None]*num_criteria
        for i in np.arange(num_criteria):
            final_roots[i] = [ a.query(ii) for ii in np.arange(a.size) if a.query(ii)['id'] in final_roots_ids[i]]


    return final_roots


def print_halo_info(a, treenodes):
    """
    Print some information about the list of treenodes

    Parameters
    -----------
    a         :   ytree object from which the nodes come from
    treenodes : list of tree nodes
    """

    for tn in treenodes:
        print("%4i %5.5E %5.5E %.6f %.6f %.6f"%(tn['id'],tn['mass'],tn['virial_radius'],tn['x']/a.box_size,tn['y']/a.box_size,tn['z']/a.box_size))

    return

if __name__ == "__main__":

    #
    #
    #
    all_criteria = [ '(tree["tree", "redshift"] > 8.0) *'
                     '(tree["tree", "redshift"] < 12.5) *'
                     '(tree["tree", "mass"] > 1.0E9)',

                     '(tree["tree", "redshift"] > 7.0) *'
                     '(tree["tree", "redshift"] < 8.1) *'
                     '(tree["tree", "mass"] > 1.0E9)',

                     '(tree["tree", "redshift"] > 6.0) *'
                     '(tree["tree", "redshift"] < 7.1) *'
                     '(tree["tree", "mass"] > 1.0E9)']

    fields = ['redshift','mass']


    a = ytree.load(str(sys.argv[1]))

    all_halos = get_root_halos(a, all_criteria, fields, unique=True, mutually_exclusive='sequential')

    for i in np.arange(len(all_criteria)):
        print_halo_info(a,all_halos[i])
        print("----------------------------")
