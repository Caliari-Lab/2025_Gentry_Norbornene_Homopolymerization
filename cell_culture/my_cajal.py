# Functions for sampling points from a 2D segmented image
import os
import warnings
from typing import List, Iterator, Tuple
import numpy as np
import numpy.typing as npt
from skimage import measure
import tifffile
from PIL import Image
from scipy.spatial.distance import pdist
import itertools as it
from pathos.pools import ProcessPool
import glob
#from .utilities import write_csv_block


import itertools as it
from typing import Iterator, Iterable, Optional, TypeVar, Generic, Union, Callable
from dataclasses import dataclass
import csv

T = TypeVar("T")

@dataclass
class Err(Generic[T]):
    code: T


def _filter_to_cells(segmask: npt.NDArray[np.int_], background: int) -> list[int]:
    """
    Return a list of identifiers for cells in the interior of the image.
    """
    cell_ids = set(np.unique(segmask))
    remove_cells = set()
    remove_cells.add(background)
    remove_cells.update(np.unique(segmask[0, :]))
    remove_cells.update(np.unique(segmask[-1, :]))
    remove_cells.update(np.unique(segmask[:, 0]))
    remove_cells.update(np.unique(segmask[:, -1]))
    return list(cell_ids.difference(remove_cells))


def cell_boundaries(
    imarray: npt.NDArray[np.int_],
    n_sample: int,
    background: int = 0,
    discard_cells_with_holes: bool = False,
    only_longest: bool = False,
    area_cutoff: int = 0,
) -> List[Tuple[int, npt.NDArray[np.float64]]]:
    """
    Sample n coordinates from the boundary of each cell in a segmented image,
    skipping cells that touch the border of the image

    :param imarray: 2D segmented image where the pixels belonging to\
          different cells have different values
    :param n_sample: number of pixel coordinates to sample from boundary of each cell
    :param background: value of background pixels, this will not be saved as a boundary
    :param discard_cells_with_holes: \
          if discard_cells_with_holes is true, we discard any cells \
          with more than one boundary (e.g., an annulus) with a \
          warning. Else, the behavior is determined by only_longest.
    :param only_longest: if discard_cells_with_holes is true, \
          only_longest is irrelevant. Otherwise, this determines whether \
          we sample points from only the longest boundary (presumably \
          the exterior) or from all boundaries, exterior and interior.
    :return:
       list of float numpy arrays of shape (n_sample, 2) \
       containing points sampled from the contours.
    """

    cell_id_list = _filter_to_cells(imarray, background)
    outlist: List[Tuple[int, npt.NDArray[np.float64]]] = []
    for cell in cell_id_list:
        cell_imarray = (imarray == cell) * 1
        if sum(sum(cell_imarray)) > area_cutoff:
            boundary_pts_list = measure.find_contours(
                cell_imarray, 0.5, fully_connected="high"
            )
            if discard_cells_with_holes and len(boundary_pts_list) > 1:
                warnings.warn("More than one boundary for cell " + str(cell))
                continue
            boundary_pts: npt.NDArray[np.float64]
            if only_longest:
                boundary_pts_list.sort(key=lambda ell: ell.shape[0])
                boundary_pts = boundary_pts_list[0]
            else:
                boundary_pts = np.concatenate(boundary_pts_list)
            if boundary_pts.shape[0] < n_sample:
                warnings.warn("Fewer than " + str(n_sample) + " pixels around boundary of cell " + str(cell))
            indices = np.linspace(0, boundary_pts.shape[0] - 1, n_sample)
            outlist.append((cell, boundary_pts[indices.astype("uint32")]))
    return list(outlist)


def _compute_intracell_all(
    infolder: str,
    file_names: list[str],
    n_sample: int,
    pool: ProcessPool,
    background: int,
    discard_cells_with_holes: bool,
    only_longest: bool,
    area_cutoff: int = 0,
) -> Iterator[Tuple[str, npt.NDArray[np.float64]]]:
#    file_names = [
#        file_name
#        for file_name in glob.glob(infolder + '/*_masks.tif')[:2]#os.listdir(infolder) #TRITC 556,600]_001
#        #if os.path.splitext(file_name)[1] in [".tif", ".tiff", ".TIF", ".TIFF"]
#    ]
    #cell_names = [os.path.splitext(file_name)[0] for file_name in file_names]
    cell_names = ['file_'+str(i) for i in range(len(file_names))] #reduces storage space and can be accessed with meta file
    # compute_cell_boundaries: Callable[[str], List[Tuple[int, npt.NDArray[np.float64]]]]
    def compute_cell_boundaries(file_name: str):
        return cell_boundaries(
            #tifffile.imread(os.path.join(infolder, file_name)),  # type: ignore
            measure.label(np.array(Image.open(os.path.join(infolder, file_name)))[:,:,0])-1, #modified to accept png file of masks, but tiff file is better
            n_sample,
            background,
            discard_cells_with_holes, 
            only_longest,
            area_cutoff #remove erroneous masks
        )

    cell_names_repeat: Iterator[Iterator[str]]
    cell_names_repeat = map(it.repeat, cell_names)
    cell_bdary_lists: Iterator[
        Tuple[Iterator[str], Iterator[Tuple[int, npt.NDArray[np.float64]]]]
    ]

    cell_bdary_lists = zip(
        cell_names_repeat, pool.imap(compute_cell_boundaries, file_names, chunksize=100)
    )
    cell_bdary_list_iters: Iterator[
        Iterator[Tuple[str, Tuple[int, npt.NDArray[np.float64]]]]
    ]
    cell_bdary_list_iters = map(lambda tup: zip(tup[0], tup[1]), cell_bdary_lists)
    cell_bdary_list_flattened: Iterator[Tuple[str, Tuple[int, npt.NDArray[np.float64]]]]
    cell_bdary_list_flattened = it.chain.from_iterable(cell_bdary_list_iters)

    def restructure_and_get_pdist(
        tup: tuple[str, tuple[int, npt.NDArray[np.float64]]]
    ) -> tuple[str, npt.NDArray[np.float64]]:
        name = tup[0] + "_" + str(tup[1][0])
        pd = pdist(tup[1][1])
        return name, pd

    return pool.imap(
        restructure_and_get_pdist, cell_bdary_list_flattened, chunksize=1000
    )


def compute_icdm_all(
    infolder: str,
    out_csv: str,
    file_names: list[str],
    n_sample: int,
    num_processes: int = 8,
    background: int = 0,
    discard_cells_with_holes: bool = False,
    only_longest: bool = False,
    area_cutoff: int = 0,
) -> None:
    """
    Read in each segmented image in a folder (assumed to be .tif), \
    save n pixel coordinates sampled from the boundary
    of each cell in the segmented image, \
    skipping cells that touch the border of the image.

    :param infolder: path to folder containing .tif files.
    :param out_csv: path to csv file to save cell boundaries.
    :param n_sample: number of pixel coordinates \
           to sample from boundary of each cell
    :param discard_cells_with_holes: \
        if discard_cells_with_holes is true, we discard any cells \
        with more than one boundary (e.g., an annulus) with a \
        warning. Else, the behavior is determined by only_longest.
    :param background: value which characterizes the color of the background pixels, \
         this will not be saved as a boundary
    :param only_longest: if discard_cells_with_holes is true, \
         only_longest is irrelevant. Otherwise, this determines whether \
         we sample points from only the longest boundary (presumably \
         the exterior) or from all boundaries, exterior and interior.

    :param num_processes: How many threads to run while sampling.
    :return: None (writes to file)
    """

    pool = ProcessPool(nodes=num_processes)
    name_dist_mat_pairs = _compute_intracell_all(
        infolder, 
        file_names, 
        n_sample, pool, background, discard_cells_with_holes, only_longest, area_cutoff
    )
    batch_size: int = 1000
    write_csv_block(out_csv, n_sample, name_dist_mat_pairs, batch_size)
    pool.close()
    pool.join()
    pool.clear()
    return None



def write_csv_block(
    out_csv: str,
    sidelength: int,
    dist_mats: Iterator[tuple[str, Union[Err[T], npt.NDArray[np.float64]]]],
    batch_size: int,
) -> list[tuple[str, Err[T]]]:
    """
    :param sidelength: The side length of all matrices in dist_mats.
    :param dist_mats: an iterator over pairs (name, arr), where arr is an
    vector-form array (rank 1) or an error code.
    """

    failed_cells: list[tuple[str, Err[T]]] = []
    with open(out_csv, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=",")
        firstline = ["cell_id"] + [
            "d_%d_%d" % (i, j) for i, j in it.combinations(range(sidelength), 2)
        ]
        csvwriter.writerow(firstline)

        while next_batch := list(it.islice(dist_mats, batch_size)):

            good_cells: list[list[Union[str, float]]] = []
            for name, cell in next_batch:

                if isinstance(cell, Err):
                    failed_cells.append((name, cell))
                else:
                    good_cells.append([name] + cell.round(decimals=3).tolist())
            csvwriter.writerows(good_cells)
    return failed_cells




from scipy.spatial.distance import squareform
from scipy.sparse import coo_array
import itertools as it
from typing import Iterator, Iterable, Optional, TypeVar, Generic, Union, Callable

from scipy.sparse.csgraph import dijkstra


def avg_shape(
    obj_names: list[str],
    gw_dist_dict: dict[tuple[str, str], float],
    iodms: dict[str, npt.NDArray[np.float64]],
    gw_coupling_mat_dict: dict[tuple[str, str], coo_array],
):
    """
    Compute capped and uncapped average distance matrices. \
    In both cases the distance matrix is rescaled so that the minimal distance between two points \
    is 1. The "capped" distance matrix has a max distance of 2.

    :param obj_names: Keys for the gw_dist_dict and iodms.
    :param gw_dist_dict: Dictionary mapping ordered pairs (cellA_name, cellB_name) \
    to Gromov-Wasserstein distances.
    :param iodms: (intra-object distance matrices) - \
    Maps object names to intra-object distance matrices. Matrices are assumed to be given \
    in vector form rather than squareform.
    :param gw_coupling_mat_dict: Dictionary mapping ordered pairs (cellA_name, cellB_name) to \
    Gromov-Wasserstein coupling matrices from cellA to cellB.
    """
    num_objects = len(obj_names)
    medoid = identify_medoid(obj_names, gw_dist_dict)
    medoid_matrix = iodms[medoid]
    # Rescale to unit step size.
    ss = step_size(medoid_matrix)
    assert ss > 0
    medoid_matrix = medoid_matrix / step_size(medoid_matrix)
    dmat_accumulator_uncapped = np.copy(medoid_matrix)
    dmat_accumulator_capped = cap(medoid_matrix, 2.0)
    others = (obj for obj in obj_names if obj != medoid)
    for obj_name in others:
        iodm = iodms[obj_name]
        # Rescale to unit step size.
        iodm = iodm / step_size(iodm)
        reoriented_iodm = squareform(
            orient(
                medoid,
                obj_name,
                squareform(iodm, force="tomatrix"),
                gw_coupling_mat_dict,
            ),
            force="tovector",
        )
        # reoriented_iodm is not a distance matrix - it is a "pseudodistance matrix".
        # If X and Y are sets and Y is a metric space, and f : X -> Y, then \
        # d_X(x0, x1) := d_Y(f(x0),f(x1)) is a pseudometric on X.
        dmat_accumulator_uncapped += reoriented_iodm
        dmat_accumulator_capped += cap(reoriented_iodm, 2.0)
    # dmat_avg_uncapped can have any positive values, but none are zero,
    # because medoid_matrix is not zero anywhere.
    # dmat_avg_capped has values between 0 and 2, exclusive.
    return (
        dmat_accumulator_capped / num_objects,
        dmat_accumulator_uncapped / num_objects,
    )


def step_size(icdm: npt.NDArray[np.float64]) -> float:
    """
    Heuristic to estimate the step size a neuron was sampled at.
    :param icdm: Vectorform distance matrix.
    """
    return float(np.min(icdm[icdm>0]))

def avg_shape_spt(
    cell_names: list[str],
    gw_dist_dict: dict[tuple[str, str], float],
    icdms: dict[str, npt.NDArray[np.float64]],
    gw_coupling_mat_dict: dict[tuple[str, str], coo_array],
    k: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Given a set of cells together with their intracell distance matrices and
    the (precomputed) pairwise GW coupling matrices between cells, construct a
    morphological "average" of cells in the cluster. This function:

    * aligns all cells in the cluster with each other using the coupling matrices
    * takes a "local average" of all intracell distance matrices, forming a
      distance matrix which models the average local connectivity structure of the neurons
    * draw a neighborh
    * draws a minimum spanning tree through the intracell distance graph,
      allowing us to visualize this average morphology

    :param cell_names: The cluster you want to take the average of,
        expressed as a list of names of cells in the cluster. These should be
        names that occur in the keys for the other dictionary arguments.
    :param gw_dist_dict: Dictionary mapping ordered pairs (cellA_name, cellB_name)
        to Gromov-Wasserstein distances between them, as returned by
        cajal.utilities.dist_mat_of_dict.
    :param icdms: (intra-cell distance matrices) -
        Maps cell names to intra-cell distance matrices. Matrices are assumed to be given
        in vector form rather than squareform. Intracell distances are
        computed by any of the sampling functions in sample_swc, sample_seg, etc.
        and are read from file by cell_iterator_csv.
    :param gw_coupling_mat_dict: Dictionary mapping ordered pairs (cellA_name, cellB_name) to
        Gromov-Wasserstein coupling matrices from cellA to cellB, with
        cellA_name < cellB_name lexicographically
    :param k: how many neighbors in the nearest-neighbors graph in step 3
    :returns: A pair (adjacency_matrix, confidence) where adjacency_matrix
        is a Numpy matrix of shape (n, n)  (where n is the number of points in each sampled cell)
        and confidence is an array of shape (n)  adjacency_matrix has values between 0 and 2.
        When "confidence" at a node in the average graph is high, the node is not
        very close to its nearest neighbor.  We can think of this as saying that
        this node in the averaged graph is a kind of poorly amalgamated blend of
        different features in different graphs.  Conversely, when confidence is
        low, and the node is close to its nearest neighbor, we interpret this as
        meaning that this node and its nearest neighbor appear together in many
        of the graphs being averaged, so this is potentially a good
        representation of some edge that really appears in many of the graphs.
    """
    dmat_avg_capped, dmat_avg_uncapped = avg_shape(
        cell_names, gw_dist_dict, icdms, gw_coupling_mat_dict
    )
    dmat_avg_uncapped = squareform(dmat_avg_uncapped)
    # So that 0s along diagonal don't get caught in min
    np.fill_diagonal(dmat_avg_uncapped, np.max(dmat_avg_uncapped))
    # When confidence at a node in the average graph is high, the node is not
    # very close to its nearest neighbor.  We can think of this as saying that
    # this node in the averaged graph is a kind of poorly amalgamated blend of
    # different features in different graphs.  Conversely, when confidence is
    # low, and the node is close to its nearest neighbor, we interpret this as
    # meaning that this node and its nearest neighbor appear together in many
    # of the graphs being averaged, so this is potentially a good
    # representation of some edge that really appears in many of the graphs.
    confidence = np.min(dmat_avg_uncapped, axis=0)
    d = squareform(dmat_avg_capped)
    G = knn_graph(d, k)
    d = np.multiply(d, G)
    # Get shortest path tree

    spt = dijkstra(d, directed=False, indices=0, return_predecessors=True)
    # Get graph representation by only keeping distances on edges from spt
    mask = np.array([True] * (d.shape[0] * d.shape[1])).reshape(d.shape)
    for i in range(1, len(spt[1])):
        if spt[1][i] == -9999:
            print("Disconnected", i)
            continue
        mask[i, spt[1][i]] = False
        mask[spt[1][i], i] = False
    retmat = squareform(dmat_avg_capped)
    retmat[mask] = 0
    return retmat, confidence


def identify_medoid(
    cell_names: list, gw_dist_dict: dict[tuple[str, str], float]
) -> str:
    """
    Identify the medoid cell in cell_names.
    """
    return cell_names[
        np.argmin(
            dist_mat_of_dict(gw_dist_dict, cell_names, as_squareform=True).sum(axis=0)
        )
    ]


def cap(a: npt.NDArray[np.float64], c: float) -> npt.NDArray[np.float64]:
    """
    Return a copy of `a` where values above `c` in `a` are replaced with `c`.
    """
    a1 = np.copy(a)
    a1[a1 >= c] = c
    return a1


def dist_mat_of_dict(
    gw_dist_dictionary: dict[tuple[str, str], float],
    cell_names: Iterable[str],
    as_squareform: bool = True,
) -> npt.NDArray[np.float64]:
    """
    Given a distance dictionary and a list of cell names, return a square distance \
    matrix containing the pairwise GW distances between all cells in `cell_names`, and \
    in the same order.\

    It is assumed that the keys in `gw_dist_dictionary` are in alphabetical order.
    """
    if cell_names is None:
        names = set()
        for key in gw_dist_dictionary:
            names.add(key[0])
            names.add(key[1])
        cell_names = sorted(names)
    dist_list: list[float] = []
    for first_cell, second_cell in it.combinations(cell_names, 2):
        first_cell, second_cell = sorted([first_cell, second_cell])
        dist_list.append(gw_dist_dictionary[(first_cell, second_cell)])
    arr = np.array(dist_list, dtype=np.float64)
    if as_squareform:
        return squareform(arr, force="tomatrix")
    return arr


def orient(
    medoid: str,
    obj_name: str,
    iodm: npt.NDArray[np.float64],
    gw_coupling_mat_dict: dict[tuple[str, str], coo_array],
) -> npt.NDArray[np.float64]:
    """
    :param medoid: String naming the medoid object, its key in iodm
    :param obj_name: String naming the object to be compared to
    :param iodm: intra-object distance matrix given in square form
    :param gw_coupling_mat_dict: maps pairs (objA_name, objB_name) to scipy COO matrices
    :return: "oriented" squareform distance matrix
    """
    if obj_name < medoid:
        gw_coupling_mat = gw_coupling_mat_dict[(obj_name, medoid)]
    else:
        gw_coupling_mat = gw_coupling_mat_dict[(medoid, obj_name)].transpose()

    i_reorder = np.argmax(gw_coupling_mat.todense(), axis=0)
    return iodm[i_reorder][:, i_reorder]


def knn_graph(dmat: npt.NDArray[np.float64], nn: int) -> npt.NDArray[np.int_]:
    """
    :param dmat: squareform distance matrix
    :param nn: (nearest neighbors) - in the returned graph, nodes v and w will be \
    connected if v is one of the `nn` nearest neighbors of w, or conversely.
    :return: A (1,0)-valued adjacency matrix for a nearest neighbors graph, same shape as dmat.
    """
    a = np.argpartition(dmat, nn + 1, axis=0)
    sidelength = dmat.shape[0]
    graph = np.zeros((sidelength, sidelength), dtype=np.int_)
    for i in range(graph.shape[1]):
        graph[a[0 : (nn + 1), i], i] = 1
    graph = np.maximum(graph, graph.T)
    np.fill_diagonal(graph, 0)
    return graph
