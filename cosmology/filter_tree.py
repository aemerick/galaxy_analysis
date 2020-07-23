import numpy as np
import yt
import ytree




def select_root_halos(a, criteria, fields):
    """
    Select desired halos
    """
    selection = a.select_halos(criteria, fields = fields)
    roots = [x.find_root() for x in selection]

    return roots


def get_unique_root_halos(a, criteria, fields,
                          selection = 'sequential'):
    """
    Loop over a list of criteria and return root halos
    that satisfy each criteria that are unique across all criteria.

    Parameters
    -----------
    selection : bool, optional
        The way to apply the uniqueness selection. 'sequential' is
        meant to be used with criteria that are (in effect) nested
        such that the first is the most restrictive and the last the most
        general, such that it is likely that all roots in each subsequent
        selection are also included in the previous.
    """

    all_roots = []
    for c in criteria:
        all_roots.append(select_root_halos(a, c, fields))

    # now make sure these are unique

    if selection == 'sequential':
        for i in np.arange(1,np.size(all_roots)):
            removed = [all_roots[i].pop(ii) for ii in np.arange(np.size(all_roos[i])) if all_roots[i][ii] in all_roots[i-1]]
            print("From %i "%(i) + " removed ", removed)
    else:
        printf("Selection type " + selection + " not yet supported")
        raise(ValueError)

    return all_roots


if __name__ == "__main__":

    #
    #
    #
    criteria = '(tree["tree", "redshift"] > 8.0) *'
               '(tree["tree", "redshift"] < 12.5) *'
               '(tree["tree", "mass"] > 1.0E9)'
    fields = ['redshift','mass']


    a = ytree.load(treefile)
