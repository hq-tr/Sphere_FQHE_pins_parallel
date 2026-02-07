# Diagonalization of two- and one-body potentials on the sphere

This routine contains exact diagonalization (ED) routines for a Hamiltonian of the form

$$\hat H = \hat V^{2bdy} + \lambda \hat V^{1bdy}$$

where $\hat V^{2bdy}$ and $\hat V^{1bdy}$ denotes the two- and one-body potentials, respectively.

This is a work in progress. Detailed documentation will be made available soon.

## Quick explanations
• Given a fixed number of electrons and a fixed number of orbitals, a basis for many-body Hilbert space is constructed. The vectors in this basis are indiced by $|\psi_1\rangle$, $|\psi_2\rangle$,...,$|\psi_d\rangle$ where $d$ is the dimension of the Hilbert space

• For each pair of indices $(i,j)$ the matrix element is calculated $h_{ij} =\langle\psi_i|\hat H|\psi_j\rangle$.

• The collection of matrix data `(i,j,h_ij)` are stored in two separate files, `rows-cols.txt` stores `i,j` in each line, and `vals.txt` stores the corresponding `h_ij` on each line.

• Using `rows-cols.txt`, a zero sparse matrix is initialized. After that, the values in `vals.txt` are updated into this matrix.

• The main reason for this method is that updating the sparse pattern of a sparse matrix is computationally expensive. It is easier to pre-allocate memory for the matrix of a given sparse pattern, given that it is known. For this reason, the matrix elements are computed and recorded first.

• Saving the matrix to files also enables us to rerun the same systems with different parameters (e.g. different values of $\lambda$ without constructing the matrix again.)

## To-do
• The output file `rows-cols.txt` may contains repeated incidents of values `(row,col)` if the data are created in different subzones. Combining the data could reduce the amount of disk space needed to store the matrix files.


## Known bugs
• The output matrix files are corrupted for 2 electrons. (For 3 electrons or more, it is fine and the results are correct.)


## Change log
See `changelog.txt`.