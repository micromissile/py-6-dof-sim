
def Cn_brick(Cnp, Cnr, p_b_rps, r_b_rps, b_m, true_airspeed_mps):

    Cn = Cnp*p_b_rps*b_m/(2*true_airspeed_mps) + Cnr*r_b_rps*b_m/(2*true_airspeed_mps)
    
    return Cn