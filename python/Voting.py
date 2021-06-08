import numpy as np
import scipy.spatial as spatial

class candidate:
    def __init__(self, points, cnt, mbt, mdt):
        self.points = points
        self.count = cnt # how many times this feature was observed in n times of sampling
        self.mean_birth_time = mbt # mean birth time over count-times this feature was observed
        self.mean_death_time = mdt # mean death time over count-times this feature was observed

    def cand_voteing(self,radius):
        # find the list of unique data points with their number of iteration
        self.points, counts = np.unique(self.points, return_counts=True, axis=0)

        point_tree = spatial.cKDTree(self.points)
        indexes = point_tree.query_ball_tree(point_tree, r=radius)
        self.votes = np.zeros(self.points.shape[0])
        for i in range(len(indexes)):
            self.votes[i] = np.sum(counts[indexes[i]])

def voting(lol_cycles,threshold,radius):
    ############################
    # lol_cycles : list of list of cycles
    # threshold : if the distance between the center of two cycles is smaller than this number, we assume them as one cycle
    # keep a track of which cycles have been already observed in voting procedure
    seen_cycles = []
    for i in range(len(lol_cycles)):
        n = len(lol_cycles[i])
        seen_cycles.append([0] * n)

    candidates_list = []

    while(True):
        flag1 = 0
        for i in range(len(lol_cycles)):
            for j in range(len(lol_cycles[i])):
                if seen_cycles[i][j] == 0:
                    center = np.mean(lol_cycles[i][j].data, axis=0)
                    seen_cycles[i][j] = 1
                    flag1 = 1
                    break
            if(flag1):
                break

        if flag1==0:
            break


        cnt1 = 0
        cycle_border = np.zeros((1,3))
        mbt = np.zeros(len(lol_cycles))
        mdt = np.zeros(len(lol_cycles))
        for i in range(len(lol_cycles)):
            for j in range(len(lol_cycles[i])):
                tmp = np.mean(lol_cycles[i][j].data, axis=0)
                dist = np.linalg.norm(tmp-center)
                if(dist<threshold):
                    cnt1 = cnt1 + 1
                    cycle_border = np.concatenate((cycle_border,lol_cycles[i][j].data),axis=0)
                    mbt[i] = lol_cycles[i][j].birthtime
                    mdt[i] = lol_cycles[i][j].deathtime
                    seen_cycles[i][j] = 1
        cycle_border = np.delete(cycle_border, 0, 0)

        # voting procedure
        obj = candidate(cycle_border,cnt1,np.sum(mbt)/cnt1,np.sum(mdt)/cnt1)

        obj.cand_voteing(radius)
        candidates_list.append(obj)

    return candidates_list