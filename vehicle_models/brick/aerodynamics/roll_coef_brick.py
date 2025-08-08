
def Cl_brick(Clp, Clr, p_b_rps, r_b_rps, b_m, true_airspeed_mps):
    '''Roll moment coefficient for dampened tumbling brick simulation'''

    Cl = Clp*p_b_rps*b_m/(2*true_airspeed_mps) + Clr*r_b_rps*b_m/(2*true_airspeed_mps)
    
    return Cl