
def Cm_brick(Cmq, q_b_rps, c_m, true_airspeed_mps):

    Cm = Cmq*q_b_rps*c_m/(2*true_airspeed_mps)
    
    return Cm