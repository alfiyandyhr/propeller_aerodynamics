"""
Written by: alfiyandyhariansyah@gmail.com
"""
from MCEVS.Analyses.Aerodynamics.BEMT.Solver import BEMTSolverOMGroup
from scipy.interpolate import interp1d
import openmdao.api as om
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class BEMTSolver(object):
    """
    docstring for BEMTSolver
    """
    def __init__(self, propDict: dict, sectionDict: dict, airfoilDict: dict):
        super(BEMTSolver, self).__init__()
        self.propDict = propDict
        self.sectionDict = sectionDict
        self.airfoilDict = airfoilDict

    def run(self, flystate, n_iter_max=500, verbose=False):

        rho_air, RPM, v_inf_list = define_flow_conditions(flystate)
        omega = RPM * 2 * np.pi / 60  # rad/s
        B = self.propDict['num_blades']

        # Interpolate CL and CD values
        CL_func = interp1d(self.airfoilDict['alpha_data'],
                           self.airfoilDict['CL_data'],
                           kind='cubic', bounds_error=False,  # Allows extrapolation
                           fill_value=1.5856)  # CL_max
        CD_func = interp1d(self.airfoilDict['alpha_data'],
                           self.airfoilDict['CD_data'],
                           kind='cubic', bounds_error=False,   # Allows extrapolation
                           fill_value=0.04603)  # CD_max

        # Discretizing the blade radius
        d_radius = self.sectionDict['radius_list'][1] - self.sectionDict['radius_list'][0]

        # Initialize thrust, torque, and power
        T = np.zeros(len(v_inf_list))
        Q = np.zeros(len(v_inf_list))
        P = np.zeros(len(v_inf_list))

        # Initialize non-dimensional parameters
        CT = np.zeros(len(v_inf_list))
        CQ = np.zeros(len(v_inf_list))
        CP = np.zeros(len(v_inf_list))
        J = np.zeros(len(v_inf_list))
        eff = np.zeros(len(v_inf_list))

        # Do for every v_inf
        for i, v_inf in enumerate(v_inf_list):
            
            for j, radius in enumerate(self.sectionDict['radius_list']):
                
                pitch = self.sectionDict['pitch_list'][j] + self.propDict['global_twist']

                # Initial inflow factors
                a = 0.1  # axial
                ao = 0.01  # angular
                
                # Convergence stopping
                stop = False
                n_iter = 0

                while not stop:

                    v_o = (1 - ao) * omega * radius  # local angular speed (with inflow)
                    v_ax = (1 + a) * v_inf  # local axial speed (with inflow)
                    v_local = np.sqrt(v_o**2 + v_ax**2)

                    phi = np.atan(v_ax / v_o)
                    AoA = pitch - phi

                    if AoA < self.airfoilDict['alpha_data'][0]:
                        CL_int = -0.5609
                        CD_int = 0.09348
                    else:
                        CL_int = CL_func(AoA)
                        CD_int = CD_func(AoA)

                    # Inifinitesimal thrust and torque
                    dT = 0.5 * rho_air * v_local**2 * (CL_int * np.cos(phi) - CD_int * np.sin(phi)) * B * self.sectionDict['chord_list'][j]
                    dQ = 0.5 * rho_air * v_local**2 * (CL_int * np.sin(phi) + CD_int * np.cos(phi)) * B * self.sectionDict['chord_list'][j] * radius

                    # New inflow factors
                    anew = dT / (rho_air * v_local**2 * (1 + a) * 4 * np.pi * radius)
                    aonew = dQ / (rho_air * v_local * (1 + a) * omega * 4 * np.pi * radius**3)

                    # Middle inflow factor values
                    a_middle = (a + anew) / 2
                    ao_middle = (ao + aonew) / 2

                    # Convergence criteria
                    n_iter += 1
                    if abs(a_middle - a) < 1E-5 and abs(ao_middle - ao) < 1E-5 or n_iter > n_iter_max:
                        stop = True

                    # Verbosity
                    if verbose:
                        print(f'iter= {n_iter}; a= {a}; ao= {ao}; a_resid= {abs(a_middle - a)}; ao_resid= {abs(ao_middle - ao)}')

                    # Update inflow factors and n_iter
                    a = a_middle
                    ao = ao_middle

                # Accumulate total thrust and torque
                T[i] += dT * d_radius
                Q[i] += dQ * d_radius

            # Power is torque times omega
            P[i] = Q[i] * omega

            # Calculate dimensionless parameters
            n = omega / (2 * np.pi)  # rotation per second
            CT[i] = T[i] / (rho_air * n**2 * self.propDict['diameter']**4)  # coefficient of thrust
            CQ[i] = Q[i] / (rho_air * n**2 * self.propDict['diameter']**5)  # coefficient of torque
            CP[i] = P[i] / (rho_air * n**3 * self.propDict['diameter']**5)  # coefficient of power
            J[i] = v_inf / (n * self.propDict['diameter'])  # advance ratio
            eff[i] = T[i] * v_inf / P[i]  # efficiency of propeller

            # Windmilling condition
            if CT[i] < 0 or CP[i] < 0:
                eff[i] = 0.0

        return T, Q, P, CT, CQ, CP, J, eff, v_inf_list


def run_BEMTSolver_MCEVS(flystate, propDict, sectionDict, airfoilDict):
        
    rho_air, RPM, v_inf_list = define_flow_conditions(flystate)

    # Initialize thrust, torque, and power
    T = np.zeros(len(v_inf_list))
    Q = np.zeros(len(v_inf_list))
    P = np.zeros(len(v_inf_list))

    # Initialize non-dimensional parameters
    CT = np.zeros(len(v_inf_list))
    CQ = np.zeros(len(v_inf_list))
    CP = np.zeros(len(v_inf_list))
    J = np.zeros(len(v_inf_list))
    eff = np.zeros(len(v_inf_list))

    for i, v_inf in enumerate(v_inf_list):

        # --- OpenMDAO problem --- #
        prob = om.Problem(reports=False)
        indeps = prob.model.add_subsystem('indeps', om.IndepVarComp(), promotes=['*'])

        # --- Design parameters --- #
        indeps.add_output('v_inf', v_inf, units='m/s')
        indeps.add_output('rpm', RPM, units='rpm')
        indeps.add_output('diameter', propDict['diameter'], units='m')
        indeps.add_output('blade_radius', propDict['diameter'] / 2, units='m')
        indeps.add_output('hub_radius', propDict['hub_radius'], units='m')
        indeps.add_output('global_twist', propDict['global_twist'], units='rad')

        for j in range(len(sectionDict['radius_list'])):
            if j == 0:
                width = sectionDict['radius_list'][j] - propDict['hub_radius']
            else:
                width = sectionDict['radius_list'][j] - sectionDict['radius_list'][j - 1]
            indeps.add_output(f'Section{j+1}|radius', sectionDict['radius_list'][j], units='m')
            indeps.add_output(f'Section{j+1}|chord', sectionDict['chord_list'][j], units='m')
            indeps.add_output(f'Section{j+1}|pitch', sectionDict['pitch_list'][j], units='rad')
            indeps.add_output(f'Section{j+1}|width', width, units='m')

        prob.model.add_subsystem('BEMT_Solver',
                                 BEMTSolverOMGroup(nblades=propDict['num_blades'],
                                                   airfoil_list=sectionDict['airfoil_list'],
                                                   rho=rho_air,
                                                   trim_rpm=False),
                                 promotes_inputs=['*'],
                                 promotes_outputs=['*'])

        prob.setup(check=False)
        prob.run_model()
        # prob.check_partials(compact_print=True)

        # Results book-keeping
        T[i] = prob.get_val('Thrust', 'N')[0]
        Q[i] = prob.get_val('Torque', 'N*m')[0]
        P[i] = prob.get_val('Power', 'W')[0]
        CT[i] = prob.get_val('CT')[0]
        CQ[i] = prob.get_val('CQ')[0]
        CP[i] = prob.get_val('CP')[0]
        J[i] = prob.get_val('J')[0]
        eff[i] = prob.get_val('eta')[0]

    return T, Q, P, CT, CQ, CP, J, eff


def define_flow_conditions(flystate: 'str'):

    allowed_states = ['takeoff', 'climb', 'cruise']
    if flystate not in allowed_states:
        raise ValueError(f"flystate must be one of {allowed_states}, got '{flystate}'")
    print(f"Flow condition set to: {flystate}")

    if flystate == 'takeoff':
        rho_air = 1.225
        RPM = 2700
        v_inf_list = np.arange(1, 61, 1)  # actual speed= 37.04 m/s

    elif flystate == 'climb':
        rho_air = 1.024
        RPM = 2650
        v_inf_list = np.arange(1, 61, 1)  # actual speed= 45.27 m/s

    elif flystate == 'cruise':
        rho_air = 0.849
        RPM = 2680
        v_inf_list = np.arange(60, 121, 1)  # actual speed= 95.17 m/s

    return rho_air, RPM, v_inf_list


if __name__ == '__main__':

    # Choose a flight state
    flystate = 'cruise'  # ['takeoff', 'climb', 'cruise']

    prop_data = pd.read_excel('aerodata.xlsx', sheet_name=0, skiprows=1, nrows=12, usecols='A:E')
    aero_data = pd.read_excel('aerodata.xlsx', sheet_name=1, skiprows=list(range(9)) + [10])

    section = prop_data['SECTION (m)']
    pitch = prop_data['PITCH (rad)']
    chord = prop_data['CHORD LENGTH (m)']
    thtoch = prop_data['Thickness to chord ratio']

    alpha_data = aero_data['alpha']
    CL_data = aero_data['CL']
    CD_data = aero_data['CD']

    propDict = {'global_twist': 15.6 * np.pi / 180,
                'num_blades': 2,
                'diameter': 1.525,
                'hub_radius': 0.10}
    sectionDict = {'radius_list': section.to_numpy(),
                   'pitch_list': pitch.to_numpy(),
                   'chord_list': chord.to_numpy(),
                   'airfoil_list': len(section) * ['Jabiru']}
    airfoilDict = {'alpha_data': alpha_data.to_numpy() * np.pi / 180,
                   'CL_data': CL_data.to_numpy(),
                   'CD_data': CD_data.to_numpy()}

    solver1 = BEMTSolver(propDict, sectionDict, airfoilDict)
    T1, Q1, P1, CT1, CQ1, CP1, J1, eff1, v_inf_list = solver1.run(flystate, verbose=False)

    T2, Q2, P2, CT2, CQ2, CP2, J2, eff2 = run_BEMTSolver_MCEVS(flystate, propDict, sectionDict, airfoilDict)

    # Plot results
    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Plot eff on the left y-axis
    color_eff = 'tab:blue'
    ax1.set_xlabel(r'Advance Ratio ($J$)')
    ax1.set_ylabel(r'Efficiency ($\eta$)', color=color_eff)
    ax1.plot(J1, eff1, 'o', color=color_eff, label=r'$\eta$ from course')
    ax1.plot(J2, eff2, '-', color=color_eff, label=r'$\eta$ from MCEVS')
    ax1.tick_params(axis='y', labelcolor=color_eff)

    # Create a second y-axis for CT and CQ (right side)
    ax2 = ax1.twinx()
    color_ct = 'tab:red'
    color_cq = 'tab:green'
    ax2.set_ylabel(r'Thrust coeff ($C_\text{T}$) or Torque coeff ($C_\text{Q}$)', color='black')
    ax2.plot(J1, CT1, 'o', color=color_ct, label=r'$C_\text{T}$ from course')
    ax2.plot(J1, CQ1, 'o', color=color_cq, label=r'$C_\text{Q}$ from course')
    ax2.plot(J2, CT2, '-', color=color_ct, label=r'$C_\text{T}$ from MCEVS')
    ax2.plot(J2, CQ2, '-', color=color_cq, label=r'$C_\text{Q}$ from MCEVS')
    ax2.tick_params(axis='y', labelcolor='black')

    # Legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='center left')

    plt.title(r'BEMT Results: $C_\text{T}$, $C_\text{Q}$, $\eta$ vs $J$')
    plt.tight_layout()
    plt.show()

    # README !!!
    # The reason for the discrepancy:
    # MCEVS BEMT implements Prandtl tip/root-loss factor
    # While this course does not do that
    # That is why MCEVS underpredicts the efficiency
