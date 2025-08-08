import math

def sphere_drag(Mach):

    if Mach <= 0.722:
        CD = 0.45*Mach**2 + 0.424   
    else:   
        CD = 2.1*math.exp(-1.16*(Mach+0.35)) - 8.9*math.exp(-2.2*(Mach+0.35)) + 0.92

    return CD