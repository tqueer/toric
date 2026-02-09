# Classification of colored path graphs with 4 vertices.

Colored Gaussian graphical models are statistical models realized by colored graphs. The only constraints is that colored used by vertices is different than the ones used by edges. When each vertex and each edge has distinct color, the model coincides with the usual Gaussian graphical model. The space of concentration matrices associated with the graph is a linear space of symmetric matrices determined by the missing edges and colors in the graph. Its inverse is the space of covariance matrices. We are interested in the polynomials vanishing on this space with a focus on when these polynomials form a toric ideal, possibly after a linear change of variables. We provide a complete description of the 44 models coming from coloring of the path on 4 nodes. Surprisingly, most of them are not toric.

This repository contains the following files:

1. P4Table.pdf: A summary table of all 44 colored path graphs on 4 vertices, including the dimension of the corresponding model and whether the vanishing ideal is toric (before after a linear change of coordinates).
2. P4Graphs.m2: Macaulay 2 code (following the frame work of https://doi.org/10.1137/21M1466943) generates the vanishing ideal of all 44 colored path graphs and outputs their generators in LaTeX format.
3. Macaulay_to_Sage.ipynb: Python code that translates the output of P4Graphs.m2 into a format suitable for SageMath computations.
4. ToricCasesP4.ipynb: Code incorporating functions from [villjulian/isToric](https://github.com/villjulian/isToric) to determine which of the 44 models are toric after a linear change of variables. The results are summarized in P4Table.pdf.

References:
1. Maraj, Aida, and Arpan Pal. "Symmetry Lie algebras of varieties with applications to algebraic statistics." arXiv preprint arXiv:2309.10741 (2023).
2. Kahle, Thomas, and Julian Vill. "Efficiently deciding if an ideal is toric after a linear coordinate change." SIAM Journal on Applied Algebra and Geometry 9, no. 4 (2025): 727-740.
3. Coons, Jane I., Aida Maraj, Pratik Misra, and Miruna-Ştefana Sorea. "Symmetrically colored Gaussian graphical models with toric vanishing ideals." SIAM Journal on Applied Algebra and Geometry 7, no. 1 (2023): 133-158.

This project was conducted by tahda queer and advised by Dr. Aida Maraj.
Special thanks to Dr. Maximilian Wiesmann and Joan Pascual Ribes for helpful discussions.

The repository is licensed under MIT, except for ToricCasesP4.ipynb, which is licensed under GPL-3.0 (see LICENSE_ToricCasesP4).
