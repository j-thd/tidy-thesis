# File to test functions of oneD.py

from basic.chamber import area_ratio_contraction
import unittest
import numpy as np

from thermo.prop import FluidProperties
import thermo.convection
import thermo.two_phase as tp
import models.one_D as oneD


class TestPrepareSinglePhaseLiquid(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_inlet = 300 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-3 # [kg/s]

        self.prep = oneD.prepare_single_phase_liquid(T_inlet=T_inlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the LAST element in the array
        exp_T = 393.36 # [K] Saturation temperature
        exp_rho = 942.94 # [kg/m^3] Density
        exp_h = 504.70e3 # [J/kg] Enthalpy
        exp_mu = 0.00023162 # [Pa*s] Viscosity
        exp_kappa = 0.68321 # [W/(m*K)] Thermal conductivity
        exp_Pr_gas = 1.438755460253802 # [-] Prandtl

        self.assertAlmostEqual(exp_T, self.prep['T'][-1], delta=exp_T*1e-5)
        self.assertAlmostEqual(exp_rho, self.prep['rho'][-1], delta=exp_rho*1e-5)
        self.assertAlmostEqual(exp_h, self.prep['h'][-1], delta=exp_h*1e-5)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][-1], delta=exp_mu*1e-4)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][-1], delta=exp_kappa*2.5e-3)
        self.assertAlmostEqual(exp_Pr_gas, self.prep['Pr'][-1], delta=exp_Pr_gas*2.5e-3)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 391.95 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*1e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from T_inlet to T_sat
        #  so the third value should 2/9 along the between saturation and inlet values
        exp_T = (393.36-300)*(2/9) + 300 # [K] 320.746666 

        self.assertAlmostEqual(exp_T, self.prep['T'][2], delta=exp_T*1e-5)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=320.746666&TLow=320.746666&TInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # NOTE: Slight difference in reference value as webbook rounds temperature to 320.75 K
        exp_rho = 989.15 # [kg/m^3]
        exp_h = 199.46e3 # [J/kg]
        exp_mu = 0.00056965 # [Pa*s]
        exp_kappa = 0.64072 # [W/(m*K)]
        exp_Pr = 3.7167902125

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-5)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*1e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*5-3)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*5e-3)
        
    def testDeltaT(self):
        # Check if delta T is evenly spaced as expected
        exp_dT = (393.36-300)/9 # [K] temperature step
        self.assertAlmostEqual(exp_dT, self.prep['dT'], delta=1e-6*exp_dT)



class TestCalcSinglePhase(unittest.TestCase):
    # Test single-phase calculations with liquid example as input
    # Results based on manual calculation in excelsheet Verification_one_d.xlsx
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_inlet = 300 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 3 # [-] 3 steps for convenience
        m_dot = 0.1e-3 # [kg/s]
        # Get the prepared thermodynamic values
        self.prep = oneD.prepare_single_phase_liquid(T_inlet=T_inlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        # set the Nusselt function
        Nu_func_liquid = thermo.convection.Nu_DB # [-] Nusselt number for liquid phase
        # Set remaining inputs
        T_wall = 500 # [K] Wall temperature\

        #Geometry
        w_h = 100e-6 # [m] Width and height
        A_channel = w_h**2 # [m^2]
        wetted_perimeter = 4*w_h # [m]
        D_hydr = 4*A_channel/wetted_perimeter # [m]

        self.res = oneD.calc_channel_single_phase(\
            T = self.prep['T'],
            Q_dot= self.prep['Q_dot'],
            rho = self.prep['rho'],
            Pr = self.prep['Pr'],
            kappa = self.prep['kappa'],
            mu = self.prep['mu'],
            p_ref=p_inlet,
            m_dot=m_dot,
            T_wall=T_wall,
            D_hydr=D_hydr,
            wetted_perimeter=wetted_perimeter,
            A_channel=A_channel,
            Nu_func=Nu_func_liquid,
            fp=fp
            )
        return super().setUp()

    def testLength(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_L0 = 0 # [m]
        exp_L1 = 2.709706792e-3 # [m]
        exp_L2 = 5.840860059e-3 # [m]

        self.assertEqual(exp_L0, self.res['L'][0])
        self.assertAlmostEqual(exp_L1, self.res['L'][1], delta=2e-3*exp_L1)
        self.assertAlmostEqual(exp_L2, self.res['L'][2], delta=2e-3*exp_L2)

    def testReynolds(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Re1 = 2.60E+03 # [-]
        exp_Re2 = 4.32E+03 # [-]

        self.assertAlmostEqual(exp_Re1, self.res['Re'][1], delta=1e-3*exp_Re1)
        self.assertAlmostEqual(exp_Re2, self.res['Re'][2], delta=1e-3*exp_Re2)

    def testNusselt(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Nu1 = 1.7665E+01 # [-]
        exp_Nu2 = 2.1533E+01 # [-]

        self.assertAlmostEqual(exp_Nu1, self.res['Nu'][1], delta=5e-3*exp_Nu1)
        self.assertAlmostEqual(exp_Nu2, self.res['Nu'][2], delta=5e-3*exp_Nu2)

    def testVelocity(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_u1 = 1.024842E+01 # [-]
        exp_u2 = 1.060513E+01 # [-]

        self.assertAlmostEqual(exp_u1, self.res['u'][1], delta=1e-5*exp_u1)
        self.assertAlmostEqual(exp_u2, self.res['u'][2], delta=1e-5*exp_u2)

class TestPrepareSinglePhaseGas(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        T_outlet = 450 # [K]
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-3 # [kg/s]

        self.prep = oneD.prepare_single_phase_gas(T_outlet=T_outlet,steps=steps, p_ref = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_T = 393.36 # [K] Saturation temperature
        exp_rho = 1.1291 # [kg/m^3] Density
        exp_h = 2706.2e3 # [J/kg] Enthalpy
        exp_mu = 1.2963e-05 # [Pa*s] Viscosity
        exp_kappa = 0.027493 # [W/(m*K)] Thermal conductivity
        exp_Pr = 1.0270253009857053 # [-] Prandtl

        self.assertAlmostEqual(exp_T, self.prep['T'][0], delta=exp_T*1e-5)
        self.assertAlmostEqual(exp_rho, self.prep['rho'][0], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][0], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][0], delta=exp_mu*2.5e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][0], delta=exp_kappa*3e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][0], delta=exp_Pr*3e-2)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 117.8 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*5e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from T_inlet to T_sat
        #  so the third value should 2/9 along the between saturation and inlet values
        exp_T = (450-393.36)*(2/9) + 393.36 # [K] 405.95

        self.assertAlmostEqual(exp_T, self.prep['T'][2], delta=exp_T*1e-5)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=450&TLow=393.36&TInc=6.293333333333332&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # NOTE: Slight difference in reference value as webbook rounds temperature to 320.75 K
        exp_rho = 1.0901 # [kg/m^3]
        exp_h = 2733.2e3 # [J/kg]
        exp_mu = 1.3454e-05 # [Pa*s]
        exp_kappa = 0.028312 # [W/(m*K)]
        exp_Pr = 1.004250430912687 # [-]

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*2e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*2.5e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*2.5e-2)
        
    def testDeltaT(self):
        # Check if delta T is evenly spaced as expected
        exp_dT = (450-393.36)/9 # [K] temperature step
        self.assertAlmostEqual(exp_dT, self.prep['dT'], delta=1e-5*exp_dT)


class TestPrepareHomogeneousTransition(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 10 # [-] Number of data point
        m_dot = 1e-4 # [kg/s]

        self.prep = oneD.prepare_homogenous_transition(steps=steps, p = p_inlet, m_dot=m_dot, fp=fp)
        return super().setUp()
    
    def testSaturationValues(self):
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_T_sat = 393.36 # [K] Saturation temperature
        exp_rho_gas = 1.1291 # [kg/m^3] Density
        exp_h_gas = 2706.2e3 # [J/kg] Enthalpy
        exp_mu_gas = 1.2963e-05 # [Pa*s] Viscosity
        exp_kappa_gas = 0.027493 # [W/(m*K)] Thermal conductivity
        exp_Pr_gas = 1.0270253009857053 # [-] Prandtl

        self.assertAlmostEqual(exp_T_sat, self.prep['T_sat'], delta=exp_T_sat*1e-5)
        self.assertAlmostEqual(exp_rho_gas, self.prep['rho'][-1], delta=exp_rho_gas*1e-4)
        self.assertAlmostEqual(exp_h_gas, self.prep['h'][-1], delta=exp_h_gas*1e-4)
        self.assertAlmostEqual(exp_mu_gas, self.prep['mu'][-1], delta=exp_mu_gas*2.5e-3)
        self.assertAlmostEqual(exp_kappa_gas, self.prep['kappa'][-1], delta=exp_kappa_gas*3e-2)
        self.assertAlmostEqual(exp_Pr_gas, self.prep['Pr'][-1], delta=exp_Pr_gas*3e-2)

        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Saturation values which must be true for the FIRST element in the array
        exp_rho_liquid = 942.94 # [kg/m^3] Density
        exp_h_liquid = 504.70e3 # [J/kg] Enthalpy
        exp_mu_liquid = 0.00023162 # [Pa*s] Viscosity
        exp_kappa_liquid = 0.68321 # [W/(m*K)] Thermal conductivity
        exp_Pr_liquid = 1.438755460253802 # [-] Prandtl

        self.assertAlmostEqual(exp_rho_liquid, self.prep['rho'][0], delta=exp_rho_liquid*1e-5)
        self.assertAlmostEqual(exp_h_liquid, self.prep['h'][0], delta=exp_h_liquid*1e-5)
        self.assertAlmostEqual(exp_mu_liquid, self.prep['mu'][0], delta=exp_mu_liquid*1e-3)
        self.assertAlmostEqual(exp_kappa_liquid, self.prep['kappa'][0], delta=exp_kappa_liquid*2e-3)
        self.assertAlmostEqual(exp_Pr_liquid, self.prep['Pr'][0], delta=exp_Pr_liquid*2e-3)

    def testTotalQDot(self):
        # The total power put in must of course equal the specific enthalpy change times mass flow
        # Reference: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=IsoBar&Digits=5&P=2&THigh=500&TLow=300&TInc=100&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        exp_total_Q_dot = 220.15 # [W] 
        res_total_Q_dot = np.sum(self.prep['Q_dot'])

        self.assertAlmostEqual(exp_total_Q_dot, res_total_Q_dot, delta=exp_total_Q_dot*5e-4)

    def testThirdValue(self):
        # Just to check if values are proper at intermediate sections
        # There are ten values in the array, so 9 steps from x=0 to x=1
        #  The third value depends both on void fraction alpha and vapour quality
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_x = 2/9 # [-] Vapour quality
        exp_alpha = 0.995826503 # [-] Void fraction

        self.assertAlmostEqual(exp_x, self.prep['x'][2], delta=exp_x*1e-5)
        self.assertAlmostEqual(exp_alpha, self.prep['alpha'][2], delta=exp_alpha*1e-7)
        # This expected value can be used as reference to check other values at this point
        # Ref: https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C7732185&Type=SatT&Digits=5&PLow=2&PHigh=2&PInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
        # Calculating in steps of 2/9 from the data points in the reference
        
        
        exp_rho = 5.059744672 # [kg/m^3]
        exp_h = 993.9222222e3 # [J/kg]
        exp_mu = 1.3876E-05 # [Pa*s]
        exp_kappa =0.030229633 # [W/(m*K)]
        exp_Pr = 1.028743655# [-]

        self.assertAlmostEqual(exp_rho, self.prep['rho'][2], delta=exp_rho*1e-4)
        self.assertAlmostEqual(exp_h, self.prep['h'][2], delta=exp_h*1e-4)
        self.assertAlmostEqual(exp_mu, self.prep['mu'][2], delta=exp_mu*2.5e-3)
        self.assertAlmostEqual(exp_kappa, self.prep['kappa'][2], delta=exp_kappa*3e-2)
        self.assertAlmostEqual(exp_Pr, self.prep['Pr'][2], delta=exp_Pr*3e-2)

class TestCalcHomogeneousTransition(unittest.TestCase):
    def setUp(self):
        # Inputs
        fp = FluidProperties('water')
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 3 # [-] Number of data point
        m_dot = 0.1e-3 # [kg/s]
        T_wall = 500 # [K]

        self.prep = oneD.prepare_homogenous_transition(steps=steps, p = p_inlet, m_dot=m_dot, fp=fp)
        
        # Nusselt functions
        Nu_func_two_phase = tp.Nu_Kandlikar_NBD_dryout # [-] Function to calculate Nusselt number (two-phase)
        Nu_func_le = thermo.convection.Nu_DB # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
        # NOTE: Can be set to one, so only the last value should change from 0 to 1. As it allows for better testing
        Nu_func_dryout = thermo.two_phase.Nu_DB_two_phase #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it

        #Geometry
        w_h = 100e-6 # [m] Width and height
        A_channel = w_h**2 # [m^2]
        wetted_perimeter = 4*w_h # [m]
        D_hydr = 4*A_channel/wetted_perimeter # [m]

        x_tp = self.prep['x']
        alpha_tp = self.prep['alpha']
        T_sat = self.prep['T_sat']
        rho_tp_l = self.prep['rho_l']
        rho_tp_g = self.prep['rho_g']
        rho_tp = self.prep['rho']
        mu_tp_l = self.prep['mu_l']
        mu_tp = self.prep['mu']
        Pr_tp_l = self.prep['Pr_l']
        Pr_tp = self.prep['Pr']
        kappa_tp_l = self.prep['kappa_l']
        kappa_tp = self.prep['kappa']
        Q_dot_tp = self.prep['Q_dot']

        self.res = oneD.calc_homogenous_transition(
            p_sat=p_inlet,
            x=x_tp,
            alpha=alpha_tp,
            T_sat=T_sat,
            rho_l=rho_tp_l,
            rho_g=rho_tp_g,
            rho=rho_tp,
            m_dot=m_dot,
            mu_l=mu_tp_l,
            mu=mu_tp,
            Pr_l=Pr_tp_l,
            Pr=Pr_tp,
            kappa_l=kappa_tp_l,
            kappa=kappa_tp,
            Q_dot=Q_dot_tp,
            T_wall=T_wall,
            D_hydr=D_hydr,
            wetted_perimeter=wetted_perimeter,
            A_channel=A_channel,
            Nu_func_tp=Nu_func_two_phase,
            Nu_func_le=Nu_func_le,
            Nu_func_dryout=Nu_func_dryout,
            fp=fp
        )

        return super().setUp()

    def testVelocity(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_u0 = 10.6051 # [m/s]
        exp_u1 = 4433.61 # [m/s]
        exp_u2 = 8856.61 # [m/s]

        self.assertAlmostEqual(exp_u0, self.res['u'][0], delta=1e-4*exp_u0)
        self.assertAlmostEqual(exp_u1, self.res['u'][1], delta=1e-4*exp_u1)
        self.assertAlmostEqual(exp_u2, self.res['u'][2], delta=1e-4*exp_u2)

    def testReynolds(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Re0 = 4.3174E+03 # [m/s]
        exp_Re1 = 7.5617E+04 # [m/s]
        exp_Re2 = 7.7143E+04 # [m/s]

        self.assertAlmostEqual(exp_Re0, self.res['Re'][0], delta=1e-3*exp_Re0)
        self.assertAlmostEqual(exp_Re1, self.res['Re'][1], delta=2.5e-3*exp_Re1)
        self.assertAlmostEqual(exp_Re2, self.res['Re'][2], delta=2.5e-3*exp_Re2)

    def testNusselt(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_Nu0 = 2.6043E+02 # [m/s]
        exp_Nu1 = 1.6577E+02 # [m/s]
        # Dry-out function is set to 1 to only change the last Nusselt number, to allow for better testing of Kandlikar function
        exp_Nu2 = 7.6009E+00 # [m/s]

        self.assertAlmostEqual(exp_Nu0, self.res['Nu'][0], delta=1e-3*exp_Nu0)
        self.assertAlmostEqual(exp_Nu1, self.res['Nu'][1], delta=2.5e-3*exp_Nu1)
        self.assertAlmostEqual(exp_Nu2, self.res['Nu'][2], delta=2e-2*exp_Nu2)

    def testNusseltLiquidEntire(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        # The Nusselt number for the case as if the flow was entirely liquid
        exp_Nu_le = 2.1533E+01 # [-] In this case Dittus Boelter, as seen in setup.

        self.assertAlmostEqual(exp_Nu_le, self.res['Nu_le'], delta=1e-3*exp_Nu_le)

    def testLength(self):
         # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_L0 = 0 # [m]
        exp_L1 = 2.2785E-03 # [m]
        exp_L2 = 5.1971E-02 # [m]

        self.assertEqual(exp_L0, self.res['L'][0])
        self.assertAlmostEqual(exp_L1, self.res['L'][1], delta=1e-3*exp_L1)
        self.assertAlmostEqual(exp_L2, self.res['L'][2], delta=2e-2*exp_L2)
        # So, 2 percent difference in length calculations

class TestRectangularMultiChannelHomogenousCalculation(unittest.TestCase):
    def setUp(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        fp = FluidProperties('water')
        p_inlet = 2e5  # [Pa] Inlet pressure
        steps = 3 # [-] Number of data point
        m_dot = 0.1e-3 # [kg/s]
        T_wall = 500 # [K]
        T_inlet = 300 # [K]
        T_outlet = 450
        channel_amount = 1 

        
        # Nusselt functions
        Nusselt_relations = {
            'Nu_func_gas': thermo.convection.Nu_DB, # [-] Function to calculate Nusselt number (gas phase)
            'Nu_func_liquid': thermo.convection.Nu_DB ,  # [-] Function to caculate Nusselt number (liquid phase)
            'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
            'Nu_func_le': thermo.convection.Nu_DB, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
            'Nu_func_dryout': tp.Nu_DB_two_phase, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
        }
        
        # Pressure drop functions. Testing 3 seperate functions to check integraotin, modularity, as well as separate functions
        pressure_drop_relations = {
            'l': oneD.calc_single_phase_frictional_pressure_drop_low_Reynolds,
            'tp': oneD.calc_two_phase_frictional_pressure_drop_low_Reynolds,
            'g': oneD.calc_single_phase_frictional_pressure_drop_high_Reynolds,
            'contraction': oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015,
        }
       

        #Geometry
        w_channel = 100e-6 # [m] Width and height
        h_channel = 100e-6 # [m] Channel height
        area_ratio_contraction = 0 # [-] Arbitrary input for contraction pressure drop
        steps = 3 

        self.prepared_values = oneD.full_homogenous_preparation(\
            T_inlet=T_inlet,
            T_outlet=T_outlet,
            m_dot=m_dot,
            p_ref=p_inlet,
            steps_l=steps,
            steps_tp=steps,
            steps_g=steps,
            fp=fp)

        self.res = oneD.rectangular_multi_channel_homogenous_calculation(\
            channel_amount=channel_amount,
            prepared_values=self.prepared_values,
            Nusselt_relations=Nusselt_relations,
            pressure_drop_relations=pressure_drop_relations,
            w_channel=w_channel,
            h_channel=h_channel,
            m_dot=m_dot,
            T_wall=T_wall,
            p_inlet=p_inlet,
            fp=fp,
            area_ratio_contraction=area_ratio_contraction)
        return super().setUp()

    def testPressureDrop(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_dP_tp = 1.9093e7 # [Pa] 
        res_dP_tp = self.res['dP_tp'] # [Pa]
        exp_dP_l = 5.884E+04 # [Pa]
        res_dP_l = self.res['dP_l']
        exp_dP_g = 9.143E+07 # [Pa] 
        res_dP_g = self.res['dP_g'] # [Pa]
        exp_dP_contraction = 7.6266E+04 # [Pa]
        res_dP_contraction = self.res['dP_contraction'] # [Pa]
        exp_dP_total = 1.106544E+08 # [Pa]
        res_dP_total = self.res['dP_total']
        # Yes, it is negative, but shouldn't matter, as long as calculation does what it is supposed to
        exp_p_chamber = -1.103781593E+08 # 2e5- res_dP_tp - res_dP_l - res_dP_g #-1.8893E+07 # [Pa]
        res_p_chamber = self.res['p_chamber'] # [Pa]

        self.assertAlmostEqual(exp_dP_l, res_dP_l, delta=exp_dP_l*1.5e-3)
        self.assertAlmostEqual(exp_dP_tp, res_dP_tp, delta=exp_dP_tp*1.5e-2)
        self.assertAlmostEqual(exp_dP_g, res_dP_g, delta=exp_dP_g*1.5e-2)
        self.assertAlmostEqual(exp_dP_contraction, res_dP_contraction, delta=exp_dP_contraction*1e-4)
        self.assertAlmostEqual(exp_dP_total, res_dP_total, delta=exp_dP_total*5e-3)
        self.assertAlmostEqual(exp_p_chamber, res_p_chamber, delta=-exp_p_chamber*6e-3)

    def testLength(self):
        # Results based on manual calculation in excelsheet Verification_one_d.xlsx
        exp_L_total = 0.067157307 # [m]
        res_L_total = self.res['L_total'] # [m]
        # Check total two-phased-based length # [m]
        res_L_tp_total = self.res['res_tp']['L'][-1] # [m]
        exp_L_tp_total = 0.051970928 # [m]
        # Also check total gas-based length
        res_L_g_total = self.res['res_g']['L'][-1] # [m]
        exp_L_g_total = 0.009345518 # [m]

        
        self.assertAlmostEqual(exp_L_total,res_L_total,delta=exp_L_total*1.5e-2)
        self.assertAlmostEqual(exp_L_tp_total,res_L_tp_total,delta=exp_L_tp_total*1.75e-2)
        self.assertAlmostEqual(exp_L_g_total,res_L_g_total,delta=exp_L_g_total*5e-3)

    def testGasDensity(self):
        exp_rho_0 = 1.1291 # [kg/m^3]
        exp_rho_1 = 1.0458 # [kg/m^3]
        exp_rho_2 = 0.97558 # [kg/m^3]

        self.assertAlmostEqual(exp_rho_0, self.res['p_g']['rho'][0], delta=exp_rho_0*1e-4)
        self.assertAlmostEqual(exp_rho_1, self.res['p_g']['rho'][1], delta=exp_rho_1*1e-4)
        self.assertAlmostEqual(exp_rho_2, self.res['p_g']['rho'][2], delta=exp_rho_2*1e-4)

    def testEnthalpy(self):
        exp_h_0 = 2706200 # [J/kg]
        exp_h_1 = 2766100 # [J/kg]
        exp_h_2 = 2824000 # [J/kg]

        self.assertAlmostEqual(exp_h_0, self.res['p_g']['h'][0], delta=exp_h_0*1e-4)
        self.assertAlmostEqual(exp_h_1, self.res['p_g']['h'][1], delta=exp_h_1*1e-4)
        self.assertAlmostEqual(exp_h_2, self.res['p_g']['h'][2], delta=exp_h_2*1e-4)

    def testQDot(self):
        # Test if total power consumption and partial consumption is correct
        exp_Q_dot_total = 11.78 # [W] Total power need in gas section\
        res_Q_dot_total = np.sum(self.res['p_g']['Q_dot'])

        self.assertAlmostEqual(exp_Q_dot_total, res_Q_dot_total, delta=exp_Q_dot_total*1e-3)

    def testKappa(self):
        exp_kappa_0 = 0.027493 # [W/(m*K)] 
        exp_kappa_1 = 0.029434 # [W/(m*K)]
        exp_kappa_2 = 0.031677 # [W/(m*K)]

        self.assertAlmostEqual(exp_kappa_0, self.res['p_g']['kappa'][0], delta=exp_kappa_0*3e-2)
        self.assertAlmostEqual(exp_kappa_1, self.res['p_g']['kappa'][1], delta=exp_kappa_1*3e-2)
        self.assertAlmostEqual(exp_kappa_2, self.res['p_g']['kappa'][2], delta=exp_kappa_2*3e-2)

    def testViscosity(self):
        exp_mu_0 = 1.2963E-05 # [Pa*s] 
        exp_mu_1 = 1.4074E-05 # [Pa*s]
        exp_mu_2 = 1.5207E-05 # [Pa*s]

        self.assertAlmostEqual(exp_mu_0, self.res['p_g']['mu'][0], delta=exp_mu_0*5e-3)
        self.assertAlmostEqual(exp_mu_1, self.res['p_g']['mu'][1], delta=exp_mu_1*5e-3)
        self.assertAlmostEqual(exp_mu_2, self.res['p_g']['mu'][2], delta=exp_mu_2*5e-3)

    def testPrandtl(self):
        exp_Pr_0 = 1.0270 # [Pa*s] 
        exp_Pr_1 = 9.8916E-01 # [Pa*s]
        exp_Pr_2 = 9.7338E-01 # [Pa*s]

        self.assertAlmostEqual(exp_Pr_0, self.res['p_g']['Pr'][0], delta=exp_Pr_0*3e-2)
        self.assertAlmostEqual(exp_Pr_1, self.res['p_g']['Pr'][1], delta=exp_Pr_1*3e-2)
        self.assertAlmostEqual(exp_Pr_2, self.res['p_g']['Pr'][2], delta=exp_Pr_2*3e-2)

    def testVelocity(self):
        # Test velocities in various phases
        exp_u_g = [\
            8856.61,
            9562.06,
            10250.31]
        res_u_g = self.res['res_g']['u'] # [m/s]
        for exp, res in zip(exp_u_g, res_u_g):
            self.assertAlmostEqual(exp,res,delta=0.5)

        exp_u_tp = [\
            10.6051,
            4433.61,
            8856.61]
        res_u_tp = self.res['res_tp']['u'] # [m/s]
        for exp, res in zip(exp_u_tp, res_u_tp):
            self.assertAlmostEqual(exp,res,delta=0.5)

        exp_u_l = [\
            10.0341,
            10.2484,
            10.6051]
        res_u_l = self.res['res_l']['u'] # [m/s]
        for exp, res in zip(exp_u_l, res_u_l):
            self.assertAlmostEqual(exp,res,delta=0.001)

    def testReynolds(self):
        # Gas phase check
        exp_Re_1 = 7.1053E+04 # [Pa*s]
        exp_Re_2 = 6.5759E+04 # [Pa*s]

        self.assertAlmostEqual(exp_Re_1, self.res['res_g']['Re'][1], delta=exp_Re_1*1e-3)
        self.assertAlmostEqual(exp_Re_2, self.res['res_g']['Re'][2], delta=exp_Re_2*1e-3)

    def testNusselt(self):
        # Gas phase check
        exp_Nu_1 = 1.7422E+02 # [Pa*s]
        exp_Nu_2 = 1.6271E+02 # [Pa*s]

        self.assertAlmostEqual(exp_Nu_1, self.res['res_g']['Nu'][1], delta=exp_Nu_1*3e-2)
        self.assertAlmostEqual(exp_Nu_2, self.res['res_g']['Nu'][2], delta=exp_Nu_2*3e-2)

class TestCalcSinglePhaseContractionPressureDropKawahara2015(unittest.TestCase):
    def testSimple(self):
        args = {'area_ratio_contraction': 1,
                'Re_Dh_downstream': 1,
                'total_dynamic_pressure': 1}
        
        exp_dP = 0 # [Pa] Area ratio of 1 gives no pressure drop of course
        res_dP = oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015(args=args)
        self.assertAlmostEqual(exp_dP, res_dP, delta=1e-30)

        args = {'area_ratio_contraction': 0.5,
                'Re_Dh_downstream': 100,
                'total_dynamic_pressure': 1}

        exp_dP = 4.510669484996588 # [Pa] Area ratio of 1 gives no pressure drop of course
        res_dP = oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015(args=args)
        self.assertAlmostEqual(exp_dP, res_dP, delta=exp_dP*1e-10)

        args = {'area_ratio_contraction': 0.5,
                'Re_Dh_downstream': 100,
                'total_dynamic_pressure': 100}

        exp_dP = 4.510669484996588e2 # [Pa] Area ratio of 1 gives no pressure drop of course
        res_dP = oneD.calc_single_phase_contraction_pressure_drop_Kawahara2015(args=args)
        self.assertAlmostEqual(exp_dP, res_dP, delta=exp_dP*1e-10)


class TestCalcWallEffectsSimple(unittest.TestCase):
    def setUp(self):
        h_conv = 100 # [W/(m^2*K)]
        w_channel_spacing = 50e-6
        kappa_wall = 100 # [W/(m*K)]
        h_channel = 100e-6 # [m]
        w_channel = 200e-6 # [m]
        T_wall_top = 1000 # [K]
        T_wall_bottom = 600 # [K]
        T_fluid = 500 # [K]
        emissivity_chip_bottom = 0.5 # [-]

        self.we = oneD.calc_wall_effect_parameters(\
            h_conv=h_conv,
            w_channel_spacing=w_channel_spacing,
            kappa_wall=kappa_wall,
            h_channel=h_channel,
            w_channel=w_channel,
            T_wall_top=T_wall_top,
            T_wall_bottom=T_wall_bottom,
            T_fluid=T_fluid)

        # Test bottom plane heat transfer

        wall_args = {
            'kappa_wall': kappa_wall,
            'T_wall_bottom': T_wall_bottom,
            'h_channel': h_channel,
            'w_channel': w_channel,
            'w_channel_spacing': w_channel_spacing,
            'emissivity_chip_bottom': emissivity_chip_bottom,
        }

        delta_L = 1e-4*np.array([1,5]) # [m] Lengths of channel sections
        self.Q_bottom_plane_heat_transfer = oneD.calc_bottom_plane_heat_balance(\
            h_conv=h_conv,
            T_fluid=T_fluid,
            we=self.we,
            wall_args=wall_args,
            delta_L=delta_L)

    def testMSimple(self):
        # Test the wall parameter
        exp_m = 200 # [1/m]
        res_m = self.we['m']

        self.assertAlmostEqual(exp_m, res_m, delta=exp_m*1e-9)

    def testC1_C2(self):
        # Test the coefficients of the boundary value problem
        exp_C1 = -20003.66656 # [K]
        exp_C2 = 500 # [K]
        res_C1 = self.we['C1'] # [K]
        res_C2 = self.we['C2'] # [K]

        self.assertAlmostEqual(exp_C1, res_C1, delta=abs(exp_C1*1e-9))
        self.assertAlmostEqual(exp_C2, res_C2, delta=exp_C2*1e-9)

    def testT_wall_side_effective(self):
        exp_T_wall_side_effective = 799.990000400 # [K]
        res_T_wall_side_effective = self.we['T_wall_side_effective']
        self.assertAlmostEqual(exp_T_wall_side_effective, res_T_wall_side_effective, delta=exp_T_wall_side_effective*1e-9)

    def testT_wall_effective(self):
        exp_T_wall_effective = 799.996666800 # [K]
        res_T_wall_effective = self.we['T_wall_effective']
        self.assertAlmostEqual(exp_T_wall_effective, res_T_wall_effective, delta=exp_T_wall_effective*1e-9)

    def test_grad_theta_L(self):
        exp_grad_theta_L = -3999533.352 # [K/m]
        res_grad_theta_L = self.we['grad_theta_L'] 
        self.assertAlmostEqual(exp_grad_theta_L, res_grad_theta_L, delta=abs(exp_grad_theta_L*1e-9))
    
    def test_bottom_plane_heat_balance(self):
        exp_Q_heat_balance = 11.99684890 # [W]
        res_Q_heat_balance = self.Q_bottom_plane_heat_transfer # [W]

        self.assertAlmostEqual(exp_Q_heat_balance, res_Q_heat_balance, delta=abs(exp_Q_heat_balance*1e-9))