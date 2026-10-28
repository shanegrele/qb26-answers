#!/usr/bin/env python

import sys

import numpy
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

#1.1: Acquiring the Data (I did this at the command line to save in my week5 directory)
#- wget https://github.com/bxlab/cmdb-quantbio/raw/refs/heads/main/assignments/lab/unsupervised_learning/extra_data/unsupervised_learning_data.tgz 
#- tar xzf unsupervised_learning_data.tgz 


def main(): #from the framework
    rapa_data = numpy.loadtxt('WT_NF_Rapa_counts.tsv')
    nf_data = numpy.loadtxt('WT_NF_Ctrl_counts.tsv')
    nf_rep_data = numpy.loadtxt('WT_NF_Ctrl_rep_counts.tsv')

    # The row and column names are stored in separate files. The sampling time and replicate can be extracted from the sample names
    genes = numpy.array([name.rstrip() for name in open('gene_names.txt')])
    sample_labels = numpy.array([name.rstrip().split('\t') for name in open('samples.txt')])
    times = numpy.array([int(time.split('_')[1]) for time, rep in sample_labels if rep == 1], numpy.int32)
    rep_times = numpy.array([int(time.split('_')[1]) for time, rep in sample_labels], numpy.int32)
    reps = numpy.array([int(rep.split('_')[1]) for time, rep in sample_labels], numpy.int32)

    # step 1.3: Print the first two PC's from step 1.2
    # I got an error message initially, and had HOPGPT suggest I place this inside def main() and outside of def plot_reps_pca() because calling main() defines plot_reps_pca() in 1.2, then calls it in here in 1.3.
    # I also used rep_times, or the value of times, to fix an issue I was having where I would get a blank plot, and Claude suggested I was getting this error because 'rep in sample_labels' is a string, and 'rep == 1' is an integer.
    plot_reps_pca(nf_rep_data, rep_times, reps, title="Replicate PCA Plot from Step 1.3", fname="step_1.3.png")
    #mTOR is a nutrient-sensing kinase that responds to feeding-based oscillations of gene expression, and mice were fed on a restricted schedule, so it makes sense that replicates from the same timepoint would be similar. This is what is shown when all gene expression from the replicates are collapsed into two components, that the colors generally are in similar positions on the plot, and are clustered together.

def plot_reps_pca(data, times, reps, title, fname): 
#step 1.2: the PCA analysis, (I asked HOPGPT first how to plot data and save those plots in python, then for help making this a function that can be called within main() when the entire script is run by typing python week5_assignment.py at the command line.)
    #create a PCA, data will be nf_rep_data when called
    #keep the first 2 components, include random_state=42
    pca = PCA(n_components=2, random_state=42)
    #transposes the matrix, flipping the axes
    pcs = pca.fit_transform(data.T)
    #create the scatterplot of PC1 and PC2, color by sampling time and size by replicate #
    fig, ax = plt.subplots()
    ax.scatter(pcs[:, 0], pcs[:, 1], c=times, cmap='Set1', s=reps * 15)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(title)
    #save the plot
    plt.savefig(fname, dpi=300, bbox_inches='tight')
   
def normalize_data(data): #from framework
    """
    Function to normalize data along the second axis of a 2D array.

    Parameters:
    -----------
    data : ndarray
           A 2D matrix containing sample data where the rows are samples and columns are features.

    Returns:
    --------
    norm_data : ndarray
                An array of the sample shape and data type as the original data with values normalized along the second axis.
    """
    norm_data = numpy.copy(data)
    norm_data -= numpy.mean(norm_data, axis=1, keepdims=True)
    norm_data /= numpy.std(norm_data, axis=1, keepdims=True)
    return norm_data


def plot_hm(*data, xlabel='', ylabel='', title='', orderby=None, showgroups=False, fname=None): #from framework
    """
    Function for plotting heatmap/matrix style data.

    Parameters:
    -----------
    data : ndarray(s)
           The array or arrays should be 2D with rows plotted along the y-axis and columns plotted along the x-axis.
           If multiple arrays are passed as positional arguments, they need to have the same number of rows.
    xlabel : str
           An x-axis label for the plotted data
    ylabel : str
           A y-axis label for the plotted data
    title : str
           A title to add above the plot
    orderby : ndarray
           A list of indices or group labels for ordering the data along the x-axis. If None, the data order is not changed
    showgroups : bool
           An indicator whether the group labels passed to `orderby` should be plotted as an additional color bar. Ignored if orderby is None.
    fname : str
           A file name to save the plot to. If None, the plot is only displayed but not saved.
    
    Returns:
    --------
    None
    """
    if orderby is None:
        order = numpy.arange(data[0].shape[0])
    else:
        order = numpy.argsort(orderby)
    if orderby is None or not showgroups:
        fig, ax = plt.subplots(1, len(data), figsize=(3, 15))
        if len(data) > 1:
            for i in range(len(data)):
                ax[i].imshow(data[i][order, :], aspect='auto')
                ax[i].set_xlabel(xlabel)
                if i > 0:
                    ax[i].set_yticks([])
            ax[0].set_ylabel(ylabel)
        else:
            ax.imshow(data[0][order, :], aspect='auto')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
    else:
        groups = orderby[order]
        fig, ax = plt.subplots(1, len(data) + 1, width_ratios=[1] + [3] * len(data), figsize=(4, 15))
        if len(data) > 1:
            for i in range(len(data)):
                ax[i + 1].imshow(data[i][order, :], aspect='auto')
                ax[i + 1].set_xlabel(xlabel)
                ax[i + 1].set_ylabel('')
                ax[i + 1].set_yticks([])
        else:
            ax[-1].imshow(data[0][order, :], aspect='auto')
            ax[-1].set_xlabel(xlabel)
            ax[-1].set_ylabel('')
            ax[-1].set_yticks([])
        ax[-1].yaxis.tick_right()
        ax2 = ax[-1].twinx()
        ax2.set_ylabel(ylabel)
        ax[0].set_ylabel('group')
        ax[0].imshow(groups.reshape(-1, 1), aspect='auto')
        bounds = numpy.r_[0, numpy.where(groups[1:] != groups[:-1])[0] + 1, groups.shape[0]] - 0.5
        ax[0].set_yticks((bounds[:-1] + bounds[1:]) / 2)
        ax[0].set_yticklabels([str(i) for i in range(bounds.shape[0] - 1)])
        ax[0].set_xticks([])
    if title:
        plt.suptitle(title)
    if fname is not None:
        plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()


def align_times(data, angle): #from framework
    """
    This function takes a 2D matrix of data and a 1D matrix of angles (radians) and shifts data along the second axis (axis=1)
    to align all rows to have an effective angle of zero.

    Parameters:
    -----------
    data : ndarray
           A 2D matrix of data where each row is a sample and each column is a measurement in a cyclical cycle
    angle : ndarray
           A 1D matrix of data of a length equal to the number of rows in data. Angles should be in radians.
           
    Returns:
    --------
    aligned : ndarray
              A new array with data shifted along the second axis to align each row to have an angle of zero
    """
    n = data.shape[1]
    offsets = n - angle * n / (numpy.pi * 2)
    offsets[numpy.where(offsets < 0)] += n
    aligned = numpy.zeros_like(data)
    for i in range(data.shape[0]):
        aligned[i, :] = (numpy.roll(data[i, :], numpy.floor(offsets[i])) * (numpy.ceil(offsets[i]) - offsets[i]) +
                         numpy.roll(data[i, :], numpy.ceil(offsets[i])) * (offsets[i] - numpy.floor(offsets[i])))
    return aligned

main()
