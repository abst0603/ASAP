import os

import matplotlib.pyplot as plt
import numpy as np
from extract_cycles import extract_cycles
from mpl_toolkits.mplot3d import axes3d, Axes3D
from multiprocessing import Pool

import asap
from Voting import voting


def multi_run_wrapper(args):
   return asap.sampling(*args)

if __name__=="__main__":
    # Read the data
    data = np.genfromtxt('../dataset/syntheticJellyfish/synjelly_gudhi.csv', delimiter=',')
    
    # iterate subsampling 100 times
    number_of_iteration = 100
    samples = [0]*number_of_iteration
    with Pool(os.cpu_count()) as p:
        samples = p.map(multi_run_wrapper, [(data, 0.15)]*100)

    # Extract cycles in every sampling iteration
    extracted_cycles = extract_cycles(samples, 3, 0.02)
    # extracted_cycles is the list of all found cycles with dimension "dim" and minimum persistence of "min_persistence"
    # is "number_of_iteration" samples

    # Apply voting procedure and discover all candidates for superbubbles
    candidates_list = voting(extracted_cycles, 0.5, 0.3)  # radius is chosen twice the radius of sampling
    # candidates_list is the list that contains all cycles in "extracted_cycles" together with their voting value for
    # every points on the border of cycles.
    # Note that all iteration of a single cycle in different sampling iterations are represented only with one object.

    # plot and color the bigger hole with red
    thirdq = np.quantile(candidates_list[0].votes, 0.3)
    points = candidates_list[0].points[candidates_list[0].votes>thirdq, :]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], marker='o', color='r')

    # plot and color the smaller hole based on number of votes
    thirdq = np.quantile(candidates_list[1].votes,0.5)
    points = candidates_list[1].points[candidates_list[1].votes > thirdq, :]
    v = candidates_list[1].votes[candidates_list[1].votes > thirdq]
    cmhot = plt.get_cmap("hot")
    cax = ax.scatter(points[:, 0], points[:, 1], points[:, 2], v, s=50, c=np.abs(v), cmap=cmhot)
    plt.show()

    print('done')
