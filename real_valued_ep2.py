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


import matplotlib.pyplot as plt
import numpy as np
import pypolsys
from scipy import integrate
from scipy import linalg
import sympy as sp


def find_exceptional_point(mu, nu, M, C, K, replacements, nu0=5):
    """
    Find EP for real-valued parameters mu and nu.

    Parameters
    ----------
    mu : symbol
        The first parameter.
    nu : symbol
        The second parameter.
    M : matrix
        The symbolic mass matrix.
    C : matrix
        The symbolic damping matrix.
    K : matrix
        The symbolic stiffness matrix.
    replacements : list
        A list referencing the variables and their values.
    nu0 : float, optional
        The initial guess of the second parameter. The default is 5.

    Returns
    -------
    sols : list
        The list of solutions.

    """
    def fun(p, x):
        """Evaluate the tracking ODE."""
        temp = linalg.solve(J(x, p), -dFdp(x, p))
        return temp.reshape((temp.shape[0],))

    def real(p, x):
        """Define the event 'being real valued' test by the ODE solver."""
        return x[1].imag

    real.terminal = False

    lamda = sp.symbols(r'$\lambda$')

    # Define the solution x and parameter p for numerical continuation
    x = (lamda, mu)
    p = nu
    p0 = nu0

    # Delete mu_ and nu_ from replacements
    replacements = replacements.copy()
    replacements_dict = dict(replacements)
    replacements_dict.pop(mu)
    replacements_dict.pop(nu)
    replacements = tuple(replacements_dict.items())

    # Get the parameterized equations
    det = (lamda**2*M + lamda*C + K).det()
    F = sp.Matrix([det, sp.diff(det, lamda)])

    # Compute the derivative of F with respect to p
    dFdp = sp.lambdify([x, p],
                       sp.diff(F, p).subs(replacements), 'numpy')
    # Compute the Jacobian of F
    J = sp.lambdify([x, p], F.jacobian(x).subs(replacements), 'numpy')

    # Find initial solutions using pypolsys
    poly = det.subs(replacements)
    poly = poly.subs(p, p0)
    pol = pypolsys.utils.fromSympy([sp.poly(poly, x),
                                    sp.poly(sp.diff(poly, lamda), x)])
    # Pass it to POLSYS_PLP
    pypolsys.polsys.init_poly(*pol)
    # Create homogeneous partition
    part = pypolsys.utils.make_h_part(2)
    # Pass it to POLSYS_PLP
    pypolsys.polsys.init_partition(*part)
    # Solve
    pypolsys.polsys.solve(1e-12, 1e-12, 0)
    # Get the roots
    r = pypolsys.polsys.myroots

    # We are not interested in the "spurious" roots
    idx = np.where(abs(r[2, :]) > 1e-06)
    sols = []

    # We track these roots until a real parameter solution is found
    for i in idx[0]:
        w0 = r[0, i]
        a_sym0 = r[1, i]
        u0 = (w0, a_sym0)

        p_span = (p0, 0)
        p_eval = np.linspace(p_span[0], p_span[1], 101)
        sol = integrate.solve_ivp(fun, t_span=p_span, y0=u0,
                                  method='RK45', events=real,
                                  t_eval=p_eval, rtol=1e-12,
                                  atol=1e-12)

        if sol.t_events[0].shape != (0,):
            for idx in range(len(sol.t_events[0])):
                if sol.y_events[0][idx, 1].real > 0:
                    sols.append([sol.y_events[0][idx, 0],
                                sol.t_events[0][idx],
                                sol.y_events[0][idx, 1]])

    return sols


def qep_to_gep(M, C, K):
    """
    Convert the matrices from a quadratic eigenvalue problem to the
    corresponding matrices of a generalized eigenvalue problem.

    Parameters
    ----------
    M : ndarray of shape (2, 2)
        Mass matrix of a quadratic eigenvalue problem.
    C : ndarray of shape (2, 2)
        Damping matrix of a quadratic eigenvalue problem.
    K : ndarray of shape (2, 2)
        Stiffness matrix of a quadratic eigenvalue problem.

    Returns
    -------
    A : ndarray of shape (4, 4)
        Left-hand side matrix of a generalized eigenvalue problem.
    B : ndarray of shape (4, 4)
        Right-hand side matrix of a generalized eigenvalue problem.

    """
    I = np.eye(2)
    Z = np.zeros((2, 2))

    A = np.block([[Z, I],
                  [-K, -C]])

    B = np.block([[I, Z],
                  [Z, M]])

    return A, B


def evaluate_qep_matrices(sol, mu, nu, M, C, K, replacements):
    # Change the values of the replacements to corresponds to t
    replacements_ep = replacements.copy()
    replacements_dict = dict(replacements_ep)
    replacements_dict.pop(mu)
    replacements_dict.pop(nu)
    replacements_ep = list(replacements_dict.items())
    replacements_ep.append((mu, sol[2].real))
    replacements_ep.append((nu, sol[1].real))

    # Evaluate the matrices at the EP
    M = np.array(M.subs(replacements_ep)).astype(np.float64)
    C = np.array(C.subs(replacements_ep)).astype(np.float64)
    K = np.array(K.subs(replacements_ep)).astype(np.float64)

    return M, C, K


def plot_eigenvalues(sol, w):
    fig, ax = plt.subplots()
    ax.scatter(w.imag, w.real, label='Eigenvalues found by scipy.linalg.eig')
    ax.scatter(sol[0].imag, sol[0].real, marker='x',
               label='Exceptional point')
    ax.set(xlabel=r'$\operatorname{Im} \lambda$ ($\mathrm{rad}\,\mathrm{s}^{-1}$)',
           ylabel=r'$\operatorname{Re} \lambda$ ($\mathrm{s}^{-1}$)')
    ax.axis('equal')
    ax.legend()
