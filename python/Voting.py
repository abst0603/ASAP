import numpy as np
import scipy.spatial as spatial
from itertools import chain


class Candidate:
    def __init__(self, points, cnt, mbt, mdt):
        """
        Initialize a candidate cycle
        :param points: points
        :param cnt: how many times this feature was observed in n times of sampling
        :param mbt: mean birth time over count-times this feature was observed
        :param mdt: mean death time over count-times this feature was observed
        """
        self.points = points
        self.count = cnt
        self.mean_birth_time = mbt
        self.mean_death_time = mdt
        self.votes = None

    def cand_voting(self, radius):
        """
        find the list of unique data points with their number of iteration
        :param radius: radius for counting neighborhood vote
        """
        self.points, counts = np.unique(self.points, return_counts=True, axis=0)
        point_tree = spatial.cKDTree(self.points)
        indexes = point_tree.query_ball_tree(point_tree, r=radius)
        self.votes = np.array([np.sum(counts[idx]) for idx in indexes])


def voting(lol_cycles, threshold, radius):
    """
    vote how many times a point belong to a cycle and aggregate all these votes
    :param lol_cycles: A list of list of all cycles extracted in extract_cycles
    :param threshold: if the distance between the center of two cycles is smaller than "threshold", we assume them as
    one cycle
    :param radius: The radius of neighborhood for counting points as votes
    :return: a list of candidate cycles with all their points and votes
    """
    lol_cycles = [lol_cycle for lol_cycle in chain.from_iterable(lol_cycles)]
    visited_cycles = [False]*len(lol_cycles)

    candidates_list = []
    while not all(visited_cycles):
        # find the first cycle that is not visited yet
        for idx, cycle in enumerate(visited_cycles):
            if not cycle:
                center = np.mean(lol_cycles[idx].data, axis=0)
                break

        # find all other cycles that can be unified to the selected cycle based on centers' distance
        cnt1 = 0
        cycle_border = np.zeros((1,3))
        mbt = np.zeros(len(lol_cycles))
        mdt = np.zeros(len(lol_cycles))
        for idx in range(len(lol_cycles)):
            if visited_cycles[idx]:
                continue
            dist = np.linalg.norm(np.mean(lol_cycles[idx].data, axis=0) - center)
            if dist < threshold:
                cnt1 = cnt1 + 1
                cycle_border = np.concatenate((cycle_border, lol_cycles[idx].data), axis=0)
                mbt[idx] = lol_cycles[idx].birthtime
                mdt[idx] = lol_cycles[idx].deathtime
                visited_cycles[idx] = True
        cycle_border = np.delete(cycle_border, 0, 0)

        # voting procedure
        obj = Candidate(cycle_border, cnt1, np.sum(mbt)/cnt1, np.sum(mdt)/cnt1)
        obj.cand_voting(radius)
        candidates_list.append(obj)

    return candidates_list
