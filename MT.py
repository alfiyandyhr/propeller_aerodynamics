"""
Written by: alfiyandyhariansyah@gmail.com
"""
from MCEVS.Analyses.Power.HoverClimb.Constant_Speed import PowerHoverClimbConstantSpeedMT
from MCEVS.Analyses.Power.HoverDescent.Constant_Speed import PowerHoverDescentConstantSpeed
import openmdao.api as om
import numpy as np


def calc_power_using_MT(state: str, rho_air, g, v_inf, r_rotor, FM, thrust_kg):

    thrust_N = thrust_kg * g
    S_disk = np.pi * r_rotor**2

    if state == 'HoverClimb':
        P_rotor = thrust_N / FM * ((v_inf / 2) + np.sqrt((v_inf / 2)**2 + thrust_N / (2 * rho_air * S_disk)))

    elif state == 'HoverDescent':
        if v_inf >= 2.0 * np.sqrt(thrust_N / (2 * rho_air * S_disk)):
            P_rotor = thrust_N / FM * ((v_inf / 2) - np.sqrt((v_inf / 2)**2 - thrust_N / (2 * rho_air * S_disk)))
        else:
            P_rotor = 1 / FM * np.sqrt((thrust_N**3) / (2 * rho_air * S_disk))
    else:
        raise ValueError('State must be either "HoverClimb" or "HoverDescent"')

    return P_rotor


def calc_power_using_MCEVS_MT(state: str, rho_air, g, v_inf, r_rotor, FM, thrust_kg):

    if state == 'HoverClimb':

        model = om.Group()
        model.add_subsystem('MT',
                            PowerHoverClimbConstantSpeedMT(N_rotor=1,
                                                           hover_FM=FM,
                                                           rho_air=rho_air,
                                                           g=g))

        prob = om.Problem(model, reports=False)
        prob.setup()

        prob.set_val('MT.Weight|takeoff', thrust_kg)
        prob.set_val('MT.LiftRotor|radius', r_rotor)
        prob.set_val('MT.Mission|hover_climb_speed', v_inf)

        prob.run_model()
        
        return prob.get_val('MT.Power|HoverClimbConstantSpeed')[0]

    elif state == 'HoverDescent':

        model = om.Group()
        model.add_subsystem('MT',
                            PowerHoverDescentConstantSpeed(N_rotor=1,
                                                           hover_FM=FM,
                                                           rho_air=rho_air,
                                                           g=g))

        prob = om.Problem(model, reports=False)
        prob.setup()

        prob.set_val('MT.Weight|takeoff', thrust_kg)
        prob.set_val('MT.LiftRotor|radius', r_rotor)
        prob.set_val('MT.Mission|hover_descent_speed', v_inf)

        prob.run_model()
        
        return prob.get_val('MT.Power|HoverDescentConstantSpeed')[0]

    else:
        raise ValueError('State must be either "HoverClimb" or "HoverDescent"')


if __name__ == '__main__':

    # Choose a flow condition
    rho_air = 1.225  # kg/m**3
    v_inf = 5.0  # m/s
    g = 9.81  # m/s**2
    
    # Assumption
    thrust_kg = 375.0  # kg
    FM = 0.75  # figure of merit
    r_rotor = 1.5  # m

    P_MT_climb = calc_power_using_MT('HoverClimb', rho_air, g, v_inf, r_rotor, FM, thrust_kg)
    P_MCEVS_MT_climb = calc_power_using_MCEVS_MT('HoverClimb', rho_air, g, v_inf, r_rotor, FM, thrust_kg)
    P_MT_descent = calc_power_using_MT('HoverDescent', rho_air, g, v_inf, r_rotor, FM, thrust_kg)
    P_MCEVS_MT_descent = calc_power_using_MCEVS_MT('HoverDescent', rho_air, g, v_inf, r_rotor, FM, thrust_kg)

    print(f'P_MT at hover climb= {P_MT_climb / 1000} kW')
    print(f'P_MCEVS_MT at hover climb= {P_MCEVS_MT_climb / 1000} kW')

    print(f'P_MT at hover descent= {P_MT_descent / 1000} kW')
    print(f'P_MCEVS_MT at hover descent= {P_MCEVS_MT_descent / 1000} kW')
