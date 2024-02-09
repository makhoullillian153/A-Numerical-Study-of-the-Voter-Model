---
editor_options: 
  markdown: 
    wrap: 72
---

This repository was cloned from the github account I used in college (lem224).

The following files include code that was used to generate and visualize results for the final report (Final Report.pdf).

## Voter-model-simulations.R

### `Neighbor`

Uniformly selects neighbor of an individual using classic movements.

**Parameters:**

-   nbhd: Matrix representing the different individuals

-   i: Row number of individual

-   j: Column number of individual

-   boundaries: If TRUE, we assume boundaries, otherwise edges and
    corners are

-   linked to the other side of the matrix

-   Default is set to TRUE

**Returns**: row and column number of a randomly selected neighbor

### `VM`

Classic voter mode, used in VM-theory.R file. Assumes boundaries.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `VM_func`

Models the time to consensus of `VM` as a function over a number of
columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns

**Returns:** Average of, variance of, second moment of consensus of each
dimension (3 vectors total)

### `CVM_marg`

Confidence voter model, marginal. Assumes boundaries

**Parameters:**

-   row: number of rows

-   col: number of columns

-   s: number of observations

**Returns:** Time to consensus of each observation

### `CVM_extr`

Confidence voter model, extremal. Assumes boundaries

**Parameters:**

-   row: number of rows

-   col: number of columns

-   s: number of observations

**Returns:** Time to consensus of each observation

### `CVM_func`

Models the time to consensus of the confidence voter models as a
function over a number of columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns. Each will be used twice: once in
    marginal, and once in extremal case.

**Returns:** Average of, variance of, second moment of consensus of each
dimension, for both marginal and extremal cases (6 vectors total).

### `VM_noBound`

Runs the same as `VM` but assume **no boundaries**.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

**Returns:** Time to consensus of each observation

### `VM_func_noBound`

Models the time to consensus of `VM_noBound` as a function over a number
of columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns

**Returns:** Average of, variance of, second moment of consensus of each
dimension (3 vectors total)

### `CVM_marg_V2`

Confidence voter model, marginal (Version 2). Runs similarly to
`CVM_marg`, but state/confidence only change if neighbor is confident.
We can also adjust the probability that an individual is convinced by a
neighbor (conditioned on that neighbor being confident). Assumes
boundaries.

**Parameters:**

-   row: Number of rows

-   col: Number of columns

-   s: Number of observations

-   change: Probability that an individual will be convinced by a
    confident neighbor (constant for each individual). Default is 1.

**Returns:** Time to consensus of each observation

### `CVM_extr_V2`

Confidence voter model, extremal (Version 2). Runs similarly to
`CVM_extr`, but state/confidence only change if neighbor is confident.
We can also adjust the probability that an individual is convinced by a
neighbor (conditioned on that neighbor being confident). Assumes
boundaries.

**Parameters:**

-   row: Number of rows

-   col: Number of columns

-   s: Number of observations

-   change: Probability that an individual will be convinced by a
    confident neighbor (constant for each individual). Default is 1.

**Returns:** Time to consensus of each observation

### `CVM_func_V2`

Models the time to consensus of the confidence voter models (Version 2)
as a function over a number of columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns. Each will be used twice: once in
    marginal, and once in extremal case.

-   change: Probability that an individual will be convinced by a
    confident neighbor (constant for each individual). Default is 1

**Returns:** Average of, variance of, second moment of consensus of each
dimension, for both marginal and extremal cases (6 vectors total).

### `VM_anim_special_1`

Animating the classic voter model, special case 1. For fun but a good
exercise with `gganimate`.

### `VM_anim_special_2`

Animating the classic voter model, special case 2. For fun but a good
exercise with `gganimate`.

### `VM_win`

Runs the same as `VM_noBound`, so assumes **no boundaries**.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Vector of which value wins consensus in each simulation

### `VM_percents`

Runs the same as `VM_noBound`, so assumes **no boundaries**.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Percent of individuals exist with state '1' at each time

### `VM_diag`

Individuals may now interact with neighbors "diagonal" to them. Uses an
internal neighbor function, which includes both diagonal movements and
classic movements. Assumes **no boundaries**.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `VM_diag_func`

Models the time to consensus of `VM_diag` as a function over a number of
columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns

**Returns:** Average of, variance of, second moment of consensus of each
dimension (3 vectors total)

### `CVM_marg_V3`

Confidence voter model, marginal (Version 3). Each individual will have
their own probability of being influenced by a confident voter. Assumes
boundaries, and only classic movements.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `CVM_extr_V3`

Confidence voter model, extremal (Version 3). Each individual will have
their own probability of being influenced by a confident voter. Assumes
boundaries, and only classic movements.

**Parameters:**

-   row: Number of rows

-   col: Number of coumns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `CVM_func_V3`

Models the time to consensus of the confidence voter models (Version 3)
as a function over a number of columns.

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns. Each will be used twice: once in
    marginal, and once in extremal case.

**Returns:** Average of, variance of, second moment of consensus of each
dimension, for both marginal and extremal cases (6 vectors total).

### `couplingVM`

Uses internal neighbor methods (no boundaries). Runs model with only
classic movements, and model with both classic and diagonal movements
side by side. For fun, but good coupling exercise.

### `myModified_percent`

My Modified Voter model. Here, the probability of being convinced is
dependent on the proportion of agreeing voters. That is, if we have 3
individuals labeled "0" and 7 voters labeled "1", a "0" voter has a
probability of 3/10 of being convinced, and a "1" voter has a
probability of 7/10 of being convinced. We only run 1 observation.
Assumes boundaries.

**Parameters:**

-   row: number of rows

-   col: Number of columns

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Percent of population labeled "1" throughout observation

### `myModified_time`

My Modified Voter model, used in VM-theory.R file. Here, the probability
of being convinced is dependent on the proportion of agreeing voters.
That is, if we have 3 individuals labeled "0" and 7 voters labeled "1",
a "0" voter has a probability of 3/10 of being convinced, and a "1"
voter has a probability of 7/10 of being convinced. Assumes boundaries.

**Parameters:**

-   row: number of rows

-   col: Number of columns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `myModified_func`

Models the time to consensus of `myModified_time` as a function over a
number of columns

**Parameters:**

-   fixed.row: Number of rows that will remain fixed through each
    simulation

-   col.range: Range of columns

**Returns:** Average of, variance of, second moment of consensus of each
dimension (3 vectors total)

### `myModified_complete`

Looks at my modified model, but on a complete graph, rather than a
square lattice graph.

**Parameters:**

-   row: number of rows

-   col: Number of columns

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `VM_complete`

Looks at the classic voter model, but on a complete graph, rather than a
square lattice graph.

**Parameters:**

-   N: Number of individuals

-   s: Number of observations

-   pOne: Probability that a spot in the matrix is initialized with '1'.
    Default is 0.5.

**Returns:** Time to consensus of each observation

### `VM_complete_func`

Models the time to consensus of VM_complete as a function over the
initial density of agreeing voters.

**Parameters:**

-   N: Number of individuals

-   s: Number of observations

**Returns:** Average time to consensus of each initial density

## VoterModel

An R package including the functions used to generate results for my
report. To install, use
`devtools::install_github("makhoullillian153/A-Numerical-Study-of-the-Voter-Model/VoterModel")`.

Functions included:

-   `CVM_extr`

-   `CVM_extr_V2`

-   `CVM_extr_V3`

-   `CVM_func`

-   `CVM_func_V2`

-   `CVM_func_V3`

-   `CVM_marg`

-   `CVM_marg_V2`

-   

-   `CVM_marg_V3`

-   `VM`

-   `VM_complete`

-   `VM_complete_func`

-   `VM_func`

-   `myModified_percent`

-   `myModified_time`

-   `neighbor`

## VM-theory.R

Brief outline of script:

-   Calculation of the time to absorption of a 3x2 classic voter model

-   Calculation of the time to absorption of a 2x2 modified voter model

-   Calculation of the time to absorption of a 3x2 modified voter model

See `myModified_time` in Voter-model-simulations.R for more detail on
the modified model.

## Matrix-generator.R

Functions:

-   `states`: gives all possible states of an nxp model

-   `transitionMatrix.Classic`: gives the transition matrix of an nxp
    classic model

-   `transitionMatrix.Modified`: gives the transition matrix of an nxp
    modified model (see `myModified_time` in "Voter model simulations.R"
    for details)

-   `calculate`: returns the expected time to absorption of a given
    transition matrix

## Graph-theory.R

Functions:

-   `transitionMatrix.RandomWalk.lattice`: Generates a transition matrix
    of a random walk on a lattice graph

-   `transitionMatrix.RandomWalk.complete`: Generates a transition
    matrix of a random walk on a complete graph

## VM-SimResults.R

The scratchwork behind the project.

## Modified-Model-Results.R

Code and scratchwork of results of modified model.
