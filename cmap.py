from MDAnalysis import Universe 
import numpy as np 
import pandas as pd
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys

def cmap_plot(contact_map, n_frames, tick_labels. dw = 50):
    fig, ax = plt.subplots()
    im = ax.imshow(contact_map / n_frames)
    fig.colorbar(im)
    length = np.shape(contact_map)[0]
    print(length)
    ax.set_xticks(np.arange(0,length, dw))
    ax.set_yticks(np.arange(0,length, dw))
    ax.set_xticklabels(tick_labels[::dw], rotation=90)
    ax.set_yticklabels(tick_labels[::dw])
    plt.show()

def main():
    help=f"""

        Usage:
            python {sys.argv[0]} [ref.pdb] [traj.xtc]

     """
    print(help)
    # Parameters
    CUTOFF = 6.0

    ref, traj = sys.argv[1], sys.argv[2] #'../00_samples/sample3/em2.pdb', '../00_samples/sample3/npt_prod_0_skip10.xtc'

    u = Universe(ref, traj)
    u_sele = u.select_atoms('protein')
    n = len(u_sele.residues)
    contact_map = np.zeros((n,n))
    n_frames = len((u.trajectory)) 
    print(f"Cutoff = {CUTOFF}, N frames = {n_frames}")

    for frame in tqdm(u.trajectory):
        coms  = u_sele.center_of_mass(compound='residues')
        idmap = distance_matrix(coms, coms)
        icmap = np.where(idmap <= CUTOFF, 1, 0)
        contact_map += icmap 

    tick_labels = [resn+str(resi)+segid for resn,resi,segid in zip(u_sele.residues.resnames, u_sele.residues.resids, u_sele.residues.segids)]
    cmap_plot(contact_map, n_frames, tick_labels)
    pd.DataFrame(contact_map, columns = tick_labels, index = tick_labels).to_csv('contact_matrix.dat')

if __name__ == "__main__":
    main()