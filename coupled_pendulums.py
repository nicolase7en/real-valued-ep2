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
l1 = 0.9534116593508425
moi1 = 0.05627128550779616
m2 = 1
l2 = 0.6749721473410476
moi2 = 1.0001497454962438
w1 = 3.7821204567028786
w2 = 2.9638804625382984
k12 = 0.1750615764438783
c1 = 1e-03
c2 = 1e-02
c12 = 1e-04

# Create the symbolic variables (use `_` suffix)
m1_, l1_, moi1_, m2_, l2_, moi2_ = sp.symbols(r'm_1 l_1 I_1 m_2 l_2 I_2',
                                              real=True, positive=True)
g_, w1_, w2_, k12_ = sp.symbols(r'g \Omega_1 \Omega_2 k_{12}', real=True,
                                positive=True)
c1_, c2_, c12_ = sp.symbols(r'c_1 c_2 c_{12}', real=True, positive=True)

mu_ = c2_
nu_ = l2_

# Match the variables with their values
replacements_ = [(m1_, m1),
                 (l1_, l1),
                 (moi1_, moi1),
                 (m2_, m2),
                 (l2_, l2),
                 (moi2_, moi2),
                 (g_, 9.80665),
                 (w1_, w1),
                 (w2_, w2),
                 (k12_, k12),
                 (c1_, c1),
                 (c2_, c2),
                 (c12_, c12)]

# Create the mass, stiffness and damping symbolic matrices
M_ = sp.Matrix([[m1_*l1_**2 + moi1_, 0],
                [0, m2_*l2_**2 + moi2_]])
K_ = sp.Matrix([[m1_*g_*l1_ + w1_**2*moi1_ + k12_, -k12_],
                [-k12_, m2_*g_*l2_ + w2_**2*moi2_ + k12_]])
C_ = sp.Matrix([[c1_ + c12_, -c12_],
                [-c12_, c2_ + c12_]])

# Find the real-valued EPs and the associated eigenvalue
sols = find_exceptional_point(mu=mu_, nu=nu_,
                              M=M_, C=C_, K=K_,
                              replacements=replacements_)

# Index to choose among the EPs found
idx = 0
sol = sols[idx]

M, C, K = evaluate_qep_matrices(sol, mu=mu_, nu=nu_,
                                M=M_, C=C_, K=K_,
                                replacements=replacements_)
A, B = qep_to_gep(M, C, K)
w, _ = linalg.eig(A, B)
plot_eigenvalues(sol, w)
