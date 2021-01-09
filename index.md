## CellAssign Code Walkthrough

This page walks through the code for [cellassign](https://github.com/Irrationone/cellassign/), a model that assigns single-cell RNA-seq data to known cell types. Information about known marker cell types is provided as input to the model in the form of a binary marker gene by cell-type matrix. 

### cellassign.R

cellassign takes in the following parameters:
- `exprs_obj`: Either a matrix representing gene expression counts or a \code{SummarizedExperiment}
- `marker_gene_info`: Information relating marker genes to cell types
- `s`: Numeric vector of cell size factors
- `min_delta`: The minimum log fold change a marker gene must be over-expressed by in its cell type
- `X`: Numeric matrix of external covariates
- `B`: Number of bases to use for RBF dispersion function
- `shrinkage`: Boolean - should the delta parameters have hierarchical shrinkage?
- `n_batches`: Number of data subsample batches to use in inference
- `dirichlet_concentration`: Dirichlet concentration parameter for cell type abundances
- `rel_tol_adam`: The change in Q function value (in pct) below which each optimization round is considered converged
- `rel_tol_em`: The change in log marginal likelihood value (in pct) below which the EM algorithm is considered converged
- `max_iter_adam`: Maximum number of ADAM iterations to perform in each M-step
- `max_iter_em`: Maximum number of EM iterations to perform
- `learning rate`: Learning rate of ADAM optimization
- `verbose`: Boolean - should running info be printed?
- `sce_assay`: The assay from teh input \code{SummarizedExperiment} to use: this assay should always represent raw counts
- `return_SCE`: Boolean - should a SingleCellExperiment be returned with the cell type annotations added?
- `num_runs`: Number of EM optimizations to perform (the one with the maximum log-marginal likelihood value will be used as the final).

Initialize the following variables for inference.
Variable      | What it is | What it represents
| ----------- | ----------- | ----------- |
| rho | matrix from `marker_gene_info` | Binary gene by cell type matrix (where a 1 indicates that a gene is a marker for a cell type, and 0 otherwise)
| Y | expression matrix from `exprs_obj` | Expression matrix of gene expression counts
| N      | # of rows in `Y`       | Number of cells
| X   | `X` matrix       | Cleaned covariate matrix
| G   | # cols in `Y`        | Number of genes
| C   | # of cols in `rho`        | Number of cell types
| P   | # of cols in `X`        | Number of covariates

Then compute size factors for each cell using `scran::computeSumFactors(t(Y))` and assign to `s`. Call `inference_tensorflow`.

### inference_tensorflow.R

inference_tensorflow takes in the following parameters:
- `Y`: initialized to shape *(null, G)*
- `rho`: initialized to shape *(null, C)*
- `s`: initialized to shape *(null)*
- `X`: initialized to shape *(null, P)*
- `G`,`C`,`N`,`P`, `B`, `shrinkage`, `verbose`, `n_batches`, `rel_tol_adam`, `rel_tol_em`, `max_iter_adam`, `max_iter_em`, `learning_rate`, `random_seed`, `min_delta`, `dirichlet_concentration`, `threads`


```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/wukathy/cellassign/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
