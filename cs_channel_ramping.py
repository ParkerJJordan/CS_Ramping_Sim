import pandas as pd
import numpy as np
from simple_pid import PID
from scipy.integrate import odeint

#Initial Conditions
time = np.linspace(0, 1, 100)
precoat_supply_lvl = 50 #%
proveup_supply_lvl = 70 #%
ix_supply_lvl = 63 #%
pix_supply_lvl = 60 #%
evap_supply_lvl = 30 #%

num_running_pu = 3
evap_supply_ds = 0.341
evap_prod_ds = 0.805

ics_dict = {
    'pclvl': precoat_supply_lvl,  
    'pulvl': proveup_supply_lvl,
    'ixlvl': ix_supply_lvl,
    'pixlvl' : pix_supply_lvl,
    'evaplvl' : evap_supply_lvl,

    'int_pu_on': num_running_pu,
    'evapsupds' : evap_supply_ds,
    'evapprodds' : evap_prod_ds
    }

ics = [proveup_supply_lvl, ix_supply_lvl, pix_supply_lvl, evap_supply_lvl]

class CornSyurpRampingSim():
    def __init__(self, ics, time, foresight=False):
        self.ics = ics
        self.time = time
        self.puflowpv = [270]
        self.ixflowpv = [320]
        self.pixflowpv = [300]
        self.evapflowpv = [315]

        self.puflowsp = []
        self.ixflowsp = []
        self.pixflowsp = []
        self.evapflowsp = []

        self.pid = PID()
        self.pid.sample_time = 1.0
        self.rampingmodel = self.simulate()


    def tanklevels(self):
        def PUsuplvl(pulvl=None, ixlvl=None, pixlvl=None, evaplvl=None):
            Vtank = 8500

            V_in = 300
            V_out, V_outsp = self.getproveupflow(pulvl=pulvl, ixlvl=ixlvl)
            dldt = (V_in - V_out)/Vtank
            return dldt
        
        def IXsuplvl(pulvl=None, ixlvl=None, pixlvl=None, evaplvl=None):
            Vtank = 8500

            V_in, V_insp = self.getproveupflow(pulvl=pulvl, ixlvl=ixlvl)
            V_out, V_outsp = self.getixflow(ixlvl=ixlvl)
            dldt = (V_in - V_out)/Vtank
            return dldt

        def PIXsuplvl(pulvl=None, ixlvl=None, pixlvl=None, evaplvl=None):
            Vtank = 8500

            V_in, V_insp = self.getixflow(ixlvl=ixlvl)
            V_out, V_outsp = self.getpixflow(pixlvl=pixlvl)
            dldt = (V_in - V_out)/Vtank
            return dldt
    
        def EVAPsuplvl(pulvl=None, ixlvl=None, pixlvl=None, evaplvl=None):
            Vtank = 8500

            V_in, V_insp = self.getpixflow(pixlvl=pixlvl)
            V_out, V_outsp = self.getevapflow(evaplvl, auto=True)
            dldt = (V_in - V_out)/Vtank
            return dldt

        def dLdt(L, t):
            pulvl, ixlvl, pixlvl, evaplvl = L
            # Create f = (proveuplvl/dt, ixlvl/dt, pixlvl/dt, evaplvl/dt)
            f = [PUsuplvl(pulvl=pulvl, ixlvl=ixlvl, pixlvl=pixlvl, evaplvl=evaplvl), 
                 IXsuplvl(pulvl=pulvl, ixlvl=ixlvl, pixlvl=pixlvl, evaplvl=evaplvl), 
                 PIXsuplvl(pulvl=pulvl, ixlvl=ixlvl, pixlvl=pixlvl, evaplvl=evaplvl), 
                 EVAPsuplvl(pulvl=pulvl, ixlvl=ixlvl, pixlvl=pixlvl, evaplvl=evaplvl)]
            return f
        
        solution = odeint(dLdt, self.ics, self.time)

        return solution
    

# Takes in the time and tank levels and calculates flow setpoints for the precoarts
    def getprecoatflow(self):
        return
    
# Takes in the time and tank levels and calculates flow setpoints for the proveups
    def getproveupflow(self, pulvl, ixlvl):
        sp_high_lim = 370
        sp_low_lim =270
        lvl_high_lim_ix = 75
        lvl_low_lim_ix = 50
        lvl_high_lim_pu = 100
        lvl_low_lim_pu = 0
        Kc, Ki, Kd = [1.0, 1.0, 0.0]
        self.pid.tunings = (Kc, Ki, Kd)
        self.pid.output_limits = (sp_low_lim, sp_high_lim)

        spfromix = ((sp_high_lim - sp_low_lim)/(lvl_high_lim_ix - lvl_low_lim_ix))*(lvl_high_lim_ix - ixlvl) + sp_low_lim
        spfrompu = ((sp_high_lim - sp_low_lim)/(lvl_high_lim_pu - lvl_low_lim_pu))*(pulvl - lvl_low_lim_pu) + sp_low_lim
        flowsp = min(spfromix, spfrompu)
        self.pid.setpoint = flowsp
        flowpv = self.pid(self.puflowpv[-1])
        self.puflowpv.append(flowpv)
        self.puflowsp.append(flowsp)
        return flowpv, flowsp

# Takes in the time and tank levels and calculates flow setpoints for the ion exhcange units
    def getixflow(self, ixlvl):
        sp_high_lim = 360
        sp_low_lim = 225
        lvl_high_lim = 68   #Tunable
        lvl_low_lim = 40    #Tunable
        Kc, Ki, Kd = [1.0, 1.0, 0.0]
        self.pid.tunings = (Kc, Ki, Kd)
        self.pid.output_limits = (sp_low_lim, sp_high_lim)

        flowsp = ((sp_high_lim - sp_low_lim)/(lvl_high_lim - lvl_low_lim))*(ixlvl - lvl_low_lim) + sp_low_lim
        self.pid.setpoint = flowsp
        flowpv = self.pid(self.ixflowpv[-1])
        self.ixflowpv.append(flowpv)
        self.ixflowsp.append(flowsp)
        return flowpv, flowsp

# Takes in the time and tank levels and calculates flow setpoints for the polishers
    def getpixflow(self, pixlvl):
        ratio_val = 5.5     #Tunable
        scaler = 50         #Tunable
        Kc, Ki, Kd = [1.0, 1.0, 0.0]
        self.pid.tunings = (Kc, Ki, Kd)

        flowsp = (pixlvl)*ratio_val + scaler
        self.pid.setpoint = flowsp
        flowpv = self.pid(self.pixflowpv[-1])
        self.pixflowpv.append(flowpv)
        self.pixflowsp.append(flowsp)
        return flowpv, flowsp

# Takes in the time and tank levels and calculates flow setpoints for the evaporator
    def getevapflow(self, evaplvl, auto=False):
        evap_auto_sp = 315  #Tunable
        sp_high_lim = 340
        sp_low_lim = 290
        lvl_high_lim = 65   #Tunable
        lvl_low_lim = 30    #Tunable
        Kc, Ki, Kd = [1.0, 1.0, 0.0]
        self.pid.tunings = (Kc, Ki, Kd)

        if auto:
            flowsp = evap_auto_sp
        else:
            flowsp = ((sp_high_lim - sp_low_lim)/(lvl_high_lim - lvl_low_lim))*(evaplvl - lvl_low_lim) + sp_low_lim
        self.pid.setpoint = flowsp
        flowpv = self.pid(self.evapflowpv[-1])
        self.evapflowpv.append(flowpv)
        self.evapflowsp.append(flowsp)
        return flowpv, flowsp

    def simulate(self):
        # for t in self.time:
        #     result = self.tanklevels(t, self.ics)

        test = self.tanklevels()
        return test

sim = CornSyurpRampingSim(ics, time)
print('Corn Syurp Ramping:')
print(sim.rampingmodel)
print()
print('Process Values:')
print('PU', sim.puflowpv)
print('IX', sim.ixflowpv)
print('PIX', sim.pixflowpv)
print('EVAP', sim.evapflowpv)
print()
print('Setpoints')
print('PU', sim.puflowsp)
print('IX', sim.ixflowsp)
print('PIX', sim.pixflowsp)
print('EVAP', sim.evapflowsp)