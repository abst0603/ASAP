import dionysus as d
import gudhi as gd
import numpy as np
from tqdm import tqdm


class cycle:
    def __init__(self, birthtime, deathtime, idx):
        """
        initialize cycle
        """
        self.birthtime = birthtime
        self.deathtime = deathtime
        self.persistence = deathtime - birthtime
        self.idx = idx
        self.data = None

    def idx2data(self, data):
        """
        index to data point
        """
        self.data = np.zeros((len(self.idx), len(data[0])))
        for cnt2, index in enumerate(self.idx):
            self.data[cnt2, ] = data[index, ]
        del self.idx


def vertex2data(alpha_complex,datashape):
    """
    extract new order based on verticies in alpha complex
    :param alpha_complex: alpha_complex
    :param datashape: datashape
    :return: reordered data
    """
    v2d = np.zeros((datashape[0],datashape[1]))
    for i in range(datashape[0]):
        v2d[i,] = np.array(alpha_complex.get_point(i))
    return v2d


def sort_cycles(cycle_list):
    """
    sort cycles based on how persistent they are
    :param cycle_list: list of cycles
    :return: sorted list of cycles
    """
    parray = np.array([cyl.persistence for cyl in cycle_list])
    sortedarg = np.argsort(parray)
    # descending order
    sorted_list = [cycle_list[sortedarg[-1-i]] for i in range(len(cycle_list))]
    return sorted_list


def extract_cycles(samples, dim, min_persistence):
    """
    extract cycles from the sampled data.
    :param samples: sampled data
    :param dim: represent the dimension of cycles we want to recover
    :param min_persistence: minimum persistency of the cycles to recover
    :return: a list of list of cycles in each sample
    """
    sc_list = []
    for cnt in tqdm(range(len(samples))):
        data = samples[cnt]
        alpha_complex = gd.AlphaComplex(points=data)
        simplex_tree = alpha_complex.create_simplex_tree()
        # Note that Gudhi changes the order of the points, thus we should extract new order based on verticies in alpha
        # complex
        v2d = vertex2data(alpha_complex, data.shape)
        # empty up the memory as soon as we don't need the object
        del alpha_complex
        simplicies = [filtered_value for filtered_value in simplex_tree.get_filtration()]
        del simplex_tree

        # Make the filtration also in Dionysus
        f = d.Filtration(simplicies)
        m = d.homology_persistence(f)
        cycle_list = []
        for i, c in enumerate(m):
            # First condition : remove empty chains
            # Second condition : Only if the dimension is what we want
            # Third condition : only if the cycle is not 0-persistent; This is literally " Death time - Birth time"
            if c and len(f[i])-1 == dim and f[i].data - f[c[len(c)-1].index].data > min_persistence:
                vertices = {v for x in c for v in f[x.index]}
                # we should make a shallow copy so the clear method will not change the results
                obj = cycle(f[c[len(c)-1].index].data, f[i].data, vertices.copy())
                obj.idx2data(v2d)  # change idx to data points
                cycle_list.append(obj)
                vertices.clear()

        cycle_list = sort_cycles(cycle_list)
        sc_list.append(cycle_list)

    return sc_list
