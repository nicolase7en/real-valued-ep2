# -*- coding: utf-8 -*-

# This file is part of `real-valued-ep2`, a library to locate an EP2 for
# real-valued parameters.
# Copyright (C) 2023  Nicolas Even, Benoit Nennig

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sympy as sp
from scipy import linalg
from real_valued_ep2 import evaluate_qep_matrices, find_exceptional_point, plot_eigenvalues, qep_to_gep

# Parameters values
m1 = 1
m2 = 1
k1 = 0.9
k2 = 0.9
k12 = 0.1
c2 = 0

# Create the symbolic variables (use `_` suffix)
m1_, m2_ = sp.symbols(r'm_1 m_2', real=True, positive=True)
k1_, k2_, k12_ = sp.symbols(r'k_1 k_2 k_{12}', real=True, positive=True)
c2_ = sp.symbols(r'c_2', real=True, positive=True)

mu_ = c2_
nu_ = k2_

# Match the variables with their values
replacements_ = [(m1_, m1),
                 (m2_, m2),
                 (k1_, k1),
                 (k2_, k2),
                 (k12_, k12),
                 (c2_, c2)]

# Create the mass, stiffness and damping symbolic matrices
M_ = sp.Matrix([[m1_,  0],
                [0, m2_]])
K_ = sp.Matrix([[k1_ + k12_,     -k12_],
                [-k12_, k2_ + k12_]])
C_ = sp.Matrix([[0,  0],
                [0, c2_]])

# Find the real-valued EPs and the associated eigenvalue
sols = find_exceptional_point(mu=mu_, nu=nu_,
                              M=M_, C=C_, K=K_,
                              replacements=replacements_)

# Index to chose among the EPs found
idx = 0
sol = sols[idx]
M, C, K = evaluate_qep_matrices(sol, mu=mu_, nu=nu_,
                                M=M_, C=C_, K=K_,
                                replacements=replacements_)

A, B = qep_to_gep(M, C, K)
w, _ = linalg.eig(A, B)
plot_eigenvalues(sol, w)
