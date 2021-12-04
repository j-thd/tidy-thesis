# File to contain all the convection proposed in literature

from basic.chamber import velocity_from_mass_flow
from thermo.prop import FluidProperties
import numpy as np

def Nu_DB(args):
    """Quick and dirty placeholder for Nusselt number from Dittus Boelter

    Args:
        args ([type]): [description]

    Returns:
        [type]: [description]
    """
    Nu = 0.023* args['Re']**0.8*args['Pr']**0.4
    return Nu

def Nusselt_Dittus_Boelter(T_inlet, T_outlet, p, D_hydraulic, L_channel, m_dot, A, fp: FluidProperties, heating=True, supressExceptions=False):
    """Calculate the Nusselt number based on the emperical relation D.10 from the TRP reader, which is in turn based on the work of Dittus and Boelter in 1930

    Args:
        T_inlet (K): Channel inlet temperature
        T_outlet (K): Channel outlet temperature (also, nozzle inlet)
        p (Pa): Channel pressure (assumed to be inlet pressure)
        D_hydraulic (m): Hydraulic diameter (for non-circular channels)
        L_channel (m): Channel length
        m_dot (kg/s): Mass flow (used to calculate velocity at the correct density)
        A (m^2): Area through which fluid flows
        fp (FluidProperties): Object to access properties of the fluid
        heating (bool, optional): Is True if the channel wall is heated. Defaults to True.
        supressExceptions (bool, optional): Surpresses the exceptions thrown if the channel parameters lie outside the range of validity. Defaults to False.

    Raises:
        ValueError: L_channel/D_hydraulic < 60 and is outside range of validity 
        ValueError: Bulk Reynolds number is outside range of validity.
        ValueError: Bulk Prandtl number is outside of range of validity.

    Returns:
        Nu (-) : Nusselt number
    """
    if (L_channel/D_hydraulic < 60):
        msg_error = "L_channel/D_hydraulic < 60 and is outside range of validity "
        if supressExceptions:
            print(msg_error)
        else:
            raise ValueError(msg_error)
    
    # Relation D.10 from TRP reader
    # The relation must be evaluated at the bulk temperature. This means that either the temperature differential must be set to a certain temperature
    # such that the inlet and outlet temperature are in agreement with the heat flow calculated, or this function is part of an iteration that continues until the temperatures
    # match the desired conditions
    # Implicitly this also means that any pressure drop due to other reason is not evaluated, but I expect to have no big effect on the Re and Pr anyway
    T_bulk = (T_inlet+T_outlet)/2 # [K]
    rho_bulk = fp.get_density(T=T_bulk, p=p) # [kg/m^3] Density
    u_bulk = velocity_from_mass_flow(A=A, m_dot=m_dot, rho=rho_bulk) # [m/s] Velocity at bulk temperature
    Re_bulk = fp.get_Reynolds_from_velocity(T=T_bulk, p=p, L_ref=D_hydraulic, u=u_bulk) # The reference length is the tube diameter, so hydraulic is used instead for this channel
    Pr_bulk = fp.get_Prandtl(T=T_bulk,p=p)

    # Check if flow paramaters are in range of validity
    if Re_bulk < 2500 or Re_bulk > 120000:
        msg_error = "Bulk Reynolds number is outside range of validity. Re_D = {:10.8f} ".format(Re_bulk)
        if supressExceptions:
            print(msg_error)
        else:
            raise ValueError(msg_error)
        

    if Pr_bulk < 0.7 or Pr_bulk > 120:
        msg_error = "Bulk Prandtl number is outside of range of validity. Pr = {:6.3f}".format(Pr_bulk)
        if supressExceptions:
            print(msg_error)
        else:
            raise ValueError(msg_error)

    # Determine the exponent of the Prandtl number, based on heating or cooling
    if heating is True:
        n = 0.4
    else:
        n = 0.3

    return 0.023 * Re_bulk**0.8 * Pr_bulk**n

def Stanton_from_Nusselt_and_velocity(Nu, T_ref, p_ref, L_ref, m_dot, A, fp: FluidProperties):
    """Returns the Stanton number based on the relations between the Nusselt number, the Prandtl number and the Reynolds numbers.

    WARNING: REFERENCE LENGTH L_ref MUST BE EQUAL TO THE SAME LENGTH THE NUSSELT NUMBER WAS DETERMINED WITH
    WARNING: REFERENCE THERMODYNAMIC STATE (T_ref, p_ref) MUST BE EQUAL TO THOSE WITH WHICH NUSSELT NUMBER WAS DETERMINED

    Args:
        Nu (-): Nusselt number
        T_ref (K): Reference temperature (WARNING: MUST BE EQUAL TO THE REFERENCE STATE NUSSELT NUMBER WAS DETERMINED WITH)
        p_ref (Pa): Reference pressure (WARNING: MUST BE EQUAL TO THE REFERENCE STATE NUSSELT NUMBER WAS DETERMINED WITH)
        m_dot (kg/s): Mass flow (used for calculating velocity at correct density)
        A (m^2): Channel area
        L_ref (m): Reference length (for calculation of Reynold's number) WARNING: MUST BE EQUAL TO THE REFERENCE LENGTH THE NUSSELT NUMBER WAS DETERMINED WITH
        fp (FluidProperties): Object to access fluid properties

    Returns:
        St (-): Stanton number.
    """
    
    Re = fp.get_Reynolds_from_mass_flow(T=T_ref,p=p_ref, L_ref=L_ref,m_dot=m_dot,A=A)
    Pr = fp.get_Prandtl(T=T_ref,p=p_ref)

    return Nu/(Re*Pr)

def Stanton_from_Nusselt_func_and_velocity(Nu_func, m_dot, A, T_ref, p_ref, L_ref, fp: FluidProperties, kwargs: dict = None):
    """Returns the Stanton number based on the relations between the Nusselt number, the Prandtl number and the Reynolds numbers.
    This time with a Nusselt function built-in to it

    WARNING: REFERENCE LENGTH L_ref MUST BE EQUAL TO THE SAME LENGTH THE NUSSELT NUMBER WAS DETERMINED WITH
    WARNING: REFERENCE THERMODYNAMIC STATE (T_ref, p_ref) MUST BE EQUAL TO THOSE WITH WHICH NUSSELT NUMBER WAS DETERMINED

    Args:
        Nu (-): Nusselt number
        T_ref (K): Reference temperature (WARNING: MUST BE EQUAL TO THE REFERENCE STATE NUSSELT NUMBER WAS DETERMINED WITH)
        p_ref (Pa): Reference tressure (WARNING: MUST BE EQUAL TO THE REFERENCE STATE NUSSELT NUMBER WAS DETERMINED WITH)
        m_dot (kg/s): Mass flow (used for calculating velocity at correct density)
        A (m^2): Channel area
        L_ref (m): Reference length (for calculation of Reynold's number) WARNING: MUST BE EQUAL TO THE REFERENCE LENGTH THE NUSSELT NUMBER WAS DETERMINED WITH
        fp (FluidProperties): Object to access fluid properties
        kwargs (dict) : Additional keywords arguments that must be passed into the Nusselt functions

    Returns:
        St (-): Stanton number.
    """
    Re = fp.get_Reynolds_from_mass_flow(T=T_ref,p=p_ref, L_ref=L_ref,m_dot=m_dot,A=A)
    Pr = fp.get_Prandtl(T=T_ref,p=p_ref)

    # List of arguments that the Nusselt function could use. Just passing them all on
    arguments = {   'T_ref': T_ref,
                    'p_ref': p_ref,
                    'Re': Re,
                    'Pr': Pr,
                    'm_dot': m_dot,
                    'A': A,
                    'L_ref': L_ref,
                    'fp': fp,
                    'kwargs': kwargs } 

    return Nu_func(args=arguments)/(Re*Pr)

def heat_transfer_coefficient_from_Nu(Nu, kappa, L_ref):
    """ Heat transfer coefficient, from Nusselt number

    Args:
        Nu (-): Nusselt number
        kappa (W/(m*K)): Conduction
        L_ref (m): Reference length (such as hydraulic diameter)

    Returns:
        h_conv (W/(m^2*K)) - Heat transfer coefficient
    """

    return Nu * kappa / L_ref 

def Nu_laminar_developed_constant_wall_temp_square(args):
    """Nusselt number for fully developed laminar flow in a square channel
    with constant wall temperature
    """
    return 3.66*np.ones_like(args['Re'])

def Nu_laminiar_developed_constant_heat_flux_square(args):
    """Nusselt number for fully developed laminiar flow in a square channel
    with constant heat flux
    """
    return 2.98

def Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water(args):
    """Nusselt number for a Re<100 in a channel with a constant wall temp and square shape

    Args:
        args (): Dictionary with flow similarity parameters

    """
    return 1058*args['Bo']**0.7*3.614 # [W/(m^2*K)] Heat transfer cofficient according to 12.164 from Carey

