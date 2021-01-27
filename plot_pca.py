# This script takes in eigenvalues from plink pca
# then plot desired two PCs
# Use this code like this:
# python plot_pca.py /data100t1/home/wanying/hla_analysis/output/subset_pca/random_subset_9448_unprune_v1.eigenvec /data100t1/home/wanying/hla_analysis/output/subset_pca/random_subset_9448_unprune_v1.png

import matplotlib.pyplot as plt
import sys

try:
    argvs = sys.argv
    input_fn = argvs[1]
    output_fn = argvs[2]
except:
    print('No input file name provided, try again using following format:')
    print('  python random_list.py [input_file_name] [output_file_name]')
    quit()


# Create a function to read in and plot data
# Or add new data to existing data (ax should not be None in this case)
def plot_pca(input_fn, ax=None):
    # Plot pc 1 and 2 by default, can change this later
    pc_index_1 = 3
    pc_index_2 = 4
    # Empty list to store pc values for plotting
    pc1_lst = []
    pc2_lst = []
    with open(input_fn, 'r') as fh:
        line = fh.readline().strip()
        while line != '':
            # + 1 in .eigenve file
            # pc1 = line.split()[pc_index_1 + 1]
            # pc2 = line.split()[pc_index_2 + 1]
            pc1 = line.split()[pc_index_1 + 3]
            pc2 = line.split()[pc_index_2 + 3]
            pc1_lst.append(float(pc1))
            pc2_lst.append(float(pc2))
            line = fh.readline().strip()

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
        ax.plot(pc1_lst, pc2_lst, 'bo', markersize=2, alpha=0.3,
                markeredgecolor='k', markeredgewidth=1)
        figure_title = 'PC' + str(pc_index_1) + ' and ' + 'PC' + str(pc_index_2)
        ax.set_title(figure_title)
        ax.set_xlabel('PC' + str(pc_index_1))
        ax.set_ylabel('PC' + str(pc_index_2))
        return fig, ax
    else:  # Add data to existing ax
        ax.plot(pc1_lst, pc2_lst, 'ro', markersize=2, alpha=0.5,
                markeredgecolor='r', markeredgewidth=1, label='Subset')
        ax.legend()

# print(len(pc1_lst), pc1_lst[:10])
# print(len(pc2_lst), pc2_lst[:10])
fig, ax = plot_pca(input_fn)
# Hightlight subset data points
# highlight_fn = '/data100t1/home/wanying/hla_analysis/output/subset_pca/random_subset_9448_unprune_v1.eigenvec'
highlight_fn = '/data100t1/home/wanying/hla_analysis/output/subset_pca/new_projection_highlight_ref.sscore'
plot_pca(highlight_fn, ax)

fig.savefig(output_fn)
