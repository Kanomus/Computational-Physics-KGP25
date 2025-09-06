# General integrator for finite bounds 
def integ(eqn : callable, bounds : tuple[float, float], method : str = "simp3", step : float = 0.001, return_list : bool = False):
    base_method = globals()[method.lower()]

    out = 0
    x = bounds[0]

    if return_list:
        out_vals = []
        idx = 0
        while x < bounds[1]:
            out_vals.append(out)
            idx += 1
            out += base_method(eqn, x, step)
            x += step
        out_vals.append(out)
        return out_vals
    
    else:    
        while x < bounds[1]:
            out += base_method(eqn, x, step)
            x += step
        return out

# Methods for integration

def rect(eqn : callable, x : float, step : float):
    return eqn(x)*step

def trapez(eqn : callable, x : float, step : float):
    return (eqn(x)+eqn(x+step))*step/2

def simp3(eqn : callable, x : float, step : float):
    return step*(eqn(x) + 4*eqn(x+(step/2)) + eqn(x + (step)))/6