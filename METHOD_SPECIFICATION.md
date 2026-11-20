# Method

Weighted FWL regression on treatment and additive covariates after internal cyclic absorption of cluster (or supplied panel unit) and time FE. Target is a TWFE coefficient, not automatically ATT.

At fixed sample, assignment and design the scalar is linear in Y, with fixed
least-squares and aggregation operators. Tests check additivity and scaling;
designs/contrasts must still be rebuilt after deletion or reassignment.
TWFE uses weighted cyclic FE absorption without giant FE dummy matrices.
Other methods currently use bounded transparent SVD designs (12 million estimated
entries), failing explicitly above the cap. No shared runtime package dependency.
