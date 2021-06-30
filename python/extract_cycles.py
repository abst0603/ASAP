import dionysus as d
import gudhi as gd
import scipy.io
import numpy as np


class cycle:
  def __init__(self, birthtime, deathtime, idx):
    self.birthtime = birthtime
    self.deathtime = deathtime
    self.persistence = deathtime - birthtime
    self.idx = idx

  def idx2data(self,data):
    self.data = np.zeros((len(self.idx), len(data[0])))
    for cnt2, index in enumerate(self.idx):
        self.data[cnt2,] = data[index,]
    del self.idx

def vertex2data(alpha_complex,datashape):
    v2d = np.zeros((datashape[0],datashape[1]))
    for i in range(datashape[0]):
        v2d[i,] = np.array(alpha_complex.get_point(i))
    return v2d

def sort_cycles(cycle_list):
    ## sort cycles based on how persistent they are
    sorted_list = []
    parray = np.zeros(len(cycle_list))
    for i in range(len(cycle_list)):
        parray[i] = cycle_list[i].persistence
    sortedarg = np.argsort(parray)

    for i in range(len(cycle_list)):
        sorted_list.append(cycle_list[sortedarg[-1-i]]) #decending order
    return sorted_list

# dim : represent the dimension of cycles we want to recover
# min_persistence : minimum persistency of the cycles to recover
def extract_cycles(samples, dim, min_persistence):
    sc_list = []
    for cnt in range(len(samples)):
        data = samples[cnt]
        alpha_complex = gd.AlphaComplex(points=data)
        simplex_tree = alpha_complex.create_simplex_tree()
        # Note that Gudhi changes the order of the points, thus we should extract new order based on verticies in alpha complex
        v2d = vertex2data(alpha_complex, data.shape)
        del alpha_complex
        simplicies = []
        for filtered_value in simplex_tree.get_filtration():
            simplicies.append(filtered_value)

        del filtered_value
        del simplex_tree

        # Make the filtration also in Dionysus
        f = d.Filtration(simplicies)
        m = d.homology_persistence(f)
        cycle_list = []
        for i,c in enumerate(m):
            # First condition : remove empty chains
            # Second condition : Only if the dimension is what we want
            # Third condition : only if the cycle is not 0-persistent; This is literally " Death time - Birth time"
            if c and len(f[i])-1 == dim and f[i].data - f[c[len(c)-1].index].data > min_persistence:
                vertices = set()
                for x in c:
                    for v in f[x.index]:
                        vertices.add(v)
                obj = cycle(f[c[len(c)-1].index].data, f[i].data, vertices.copy()) # we should make a shallow copy so the clear method will not change the resutls
                obj.idx2data(v2d) # change idx to data points
                cycle_list.append(obj)
                vertices.clear()

        cycle_list = sort_cycles(cycle_list)
        sc_list.append(cycle_list)
        print(cnt)
        # if cnt ==10:
        #     break

    return sc_list
