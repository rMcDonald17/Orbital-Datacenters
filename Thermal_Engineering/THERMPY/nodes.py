"""
nodes.py -- Layer 2 lumped-parameter node network (Part C).

Ten nodes: eight diffusion, two boundary. Thermal-Desktop/SINDA conventions --
linear conductors in W/K, radiation conductors carry sigma so they are W/K^4,
and the pumped loop uses a directional advective term that is NOT a symmetric
conductance.

    C_i dT_i/dt = Q_i(t) + sum_j G_ij (T_j - T_i) + sum_j R_ij (T_j^4 - T_i^4)

PLACEHOLDERS, flagged in Network below and in the doc:
  - h_cp_A     cold-plate conductance: needs the spreading-resistance sub-model
  - UA_rad     loop-to-manifold conductance
  - C_* several capacitances are order-of-magnitude estimates
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

SIGMA = 5.670374419e-8
T_SPACE = 3.0

NODE_NAMES = ['die+plate', 'loop cold', 'loop hot', 'manifold+evap',
              'radiator fin', 'bus', 'array', 'battery']
N_DIFF = 8
SPACE = 8

# Earth is deliberately NOT a node. Writing a surface as "emit sigma*eps*T^4 to
# everything, absorb eps*q_IR*F from Earth" is algebraically identical to "emit
# to space over (1-F), exchange with Earth over F", because q_IR = sigma*T_e^4
# defines Earth's effective radiating temperature. The two agree to 3e-6 W/m2
# across the full range including deep cold, and the sign handles itself: net
# rejection crosses zero at -99.3 C, below which the flux form simply returns a
# negative number. Earth therefore enters only through env_flux().


@dataclass
class Network:
    """Node capacitances, conductors and loop parameters."""
    # --- capacitances, J/K
    C_die: float = 1.8e4
    C_loop_cold: float = 1.2e4
    C_loop_hot: float = 1.2e4
    C_manifold: float = 6.0e3
    C_fin: float = 2.76e5            # 67 m2 panel x 4.58 kg/m2 x 900 J/kgK
    C_bus: float = 9.0e4
    C_array: float = 4.0e4
    C_batt: float = 5.0e4

    # --- linear conductors, W/K
    h_cp_A: float = 2.0e3            # PLACEHOLDER -- spreading resistance model
    UA_rad: float = 8.0e3            # PLACEHOLDER
    G_fin: float = 1.0e4             # manifold -> fin; heat pipe + fin, eta folded in
    G_bus_batt: float = 15.0
    G_bus_manifold: float = 5.0

    # --- loop
    mdot: float = 0.90               # kg/s ammonia
    cp: float = 4700.0               # J/kgK
    eta_fin: float = 0.753           # from the fin sub-model

    # --- radiating areas, m2 (radiating surface, not physical panel)
    A_fin: float = 134.0
    A_array: float = 232.0
    A_bus: float = 12.0
    eps_fin: float = 0.85
    eps_array: float = 0.85
    eps_bus_mli: float = 0.02        # effective emittance through MLI

    # --- view factors to other spacecraft surfaces (0 until geometry defined)
    F_array_fin: float = 0.0

    @property
    def mdot_cp(self):
        return self.mdot * self.cp

    @property
    def C(self):
        return np.array([self.C_die, self.C_loop_cold, self.C_loop_hot,
                         self.C_manifold, self.C_fin, self.C_bus,
                         self.C_array, self.C_batt])

    def hx_effectiveness(self):
        """
        Effectiveness-NTU for the loop-to-manifold exchange.

        A two-node loop cannot resolve the temperature profile along the
        manifold, so writing UA*(T_outlet - T_manifold) understates the driving
        dT by the full loop rise. Using effectiveness with the INLET temperature
        is the correct two-node form and stays right if mdot or UA change.
        """
        ntu = self.UA_rad / self.mdot_cp
        return 1.0 - np.exp(-ntu)


def derivatives(t, T, net: Network, loads, env_flux):
    """
    Right-hand side of the node ODEs. T is the 8 diffusion node temperatures.

    loads(t)    -> dict of internal dissipation per node index, W
    env_flux(t) -> dict of absorbed environmental flux per node index, W
    """
    d = np.zeros(N_DIFF)
    Q = np.zeros(N_DIFF)
    for i, q in loads(t).items():
        Q[i] += q
    for i, q in env_flux(t).items():
        Q[i] += q

    die, cold, hot, man, fin, bus, arr, bat = range(8)

    # --- pumped loop. Advective terms are DIRECTIONAL: heat moves with the
    #     fluid, not down the gradient, so this is not a symmetric conductance.
    Q[hot] += net.mdot_cp * (T[cold] - T[hot])
    Q[cold] += net.mdot_cp * (T[hot] - T[cold])

    # cold plate: payload -> hot leg
    q_cp = net.h_cp_A * (T[die] - T[hot])
    Q[die] -= q_cp
    Q[hot] += q_cp

    # loop -> manifold, effectiveness-NTU on the INLET temperature.
    # The exchange happens as fluid transits hot leg -> cold leg, so the heat
    # leaves the COLD node. Placing it on the hot node forces T_hot == T_cold
    # at steady state and the loop rise vanishes.
    q_hx = net.hx_effectiveness() * net.mdot_cp * (T[hot] - T[man])
    Q[cold] -= q_hx
    Q[man] += q_hx

    # manifold -> fin (heat pipe + fin, eta_net folded into G_fin)
    q_fin = net.G_fin * (T[man] - T[fin])
    Q[man] -= q_fin
    Q[fin] += q_fin

    # structure
    for a, b, G in [(bus, bat, net.G_bus_batt), (bus, man, net.G_bus_manifold)]:
        q = G * (T[a] - T[b]); Q[a] -= q; Q[b] += q

    # --- radiation to space
    #  A_fin is already radiating area; the array is physical, so it gets 2x
    #  for front and back. Radiating from one face only puts it at ~91 C.
    for i, A, eps in [(fin, net.A_fin, net.eps_fin),
                      (arr, 2 * net.A_array, net.eps_array),
                      (bus, net.A_bus, net.eps_bus_mli)]:
        Q[i] -= eps * SIGMA * A * (T[i] ** 4 - T_SPACE ** 4)

    # --- array <-> fin coupling. Zero until vehicle geometry is defined.
    if net.F_array_fin > 0:
        q = net.eps_array * SIGMA * net.A_array * net.F_array_fin * \
            (T[arr] ** 4 - T[fin] ** 4)
        Q[arr] -= q
        Q[fin] += q

    return Q / net.C


def integrate(net: Network, T0, loads, env_flux, t_end, method='Radau'):
    """
    Integrate to t_end. Stiff by construction -- the loop responds in ~3 s and
    the radiator in ~330 s, a 117:1 ratio before the die node is considered --
    so an explicit method is not an option.
    """
    return solve_ivp(derivatives, [0, t_end], np.asarray(T0, float),
                     args=(net, loads, env_flux), method=method,
                     rtol=1e-6, atol=1e-3, dense_output=True)


def steady_state(net: Network, loads, env_flux, T0=None, t_settle=20000.0):
    """Layer 1 check: integrate to equilibrium and return node temperatures."""
    T0 = T0 if T0 is not None else np.full(N_DIFF, 300.0)
    s = integrate(net, T0, loads, env_flux, t_settle)
    return s.y[:, -1]
