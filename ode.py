from linear import solve_system as solve_sys
import matrix

# Methods to solve ordinary differential equations

# All arguments should follow the order (x, y, yd, ydd, ...)

# function to convert given ode to appropriate system for use in these functions

def make_system(eqn : callable, order : int) -> callable :
    # this system is to be returned
    def system(x : float, y : list) -> list :
        yd = [0.0]*order
        yd[:-1] = y[1:]
        yd[-1] = eqn(x, *y)
        return yd
    return system
   
# general ode solver

def ode(order : int, method : str, eqn : callable, indep_var_range : tuple, init_values : list, h : float = 0.1, tolerance : float|None = None) -> list[list]:
    base_method = globals()[method.lower()]
    system = make_system(eqn, order)
    vars = list(init_values)
    vars_list = [[] for i in range(len(vars))]

    for i in range(len(vars)):
        vars_list[i].append(vars[i])

    while vars[0] > indep_var_range[0]:
        vars, h = base_method(system, vars[0], vars[1:], h, direction = False, tolerance = tolerance)
        for i in range(len(vars)):
            vars_list[i].insert(0, vars[i])

    vars = list(init_values)

    while vars[0] < indep_var_range[1]:
        vars, h = base_method(system, vars[0], vars[1:], h, direction = True, tolerance = tolerance)
        for i in range(len(vars)):
            vars_list[i].append(vars[i])
    
    return vars_list

# shooting method

def bound_ode(order : int, method : str, eqn : callable, indep_var_range : tuple, points : list[list], h : float = 0.01, tolerance : float|None = None, **kwargs) -> list[list]:
    
    init_point = list(points[0])
    unknown_indices = []
    for i in range(len(points[0])):
        if points[0][i] == None: unknown_indices.append(i)
    
    if len(unknown_indices) == 0: return ode(order, method, eqn, indep_var_range, init_point, h, tolerance)
        
    req_indices = [[] for _ in range(len(points) - 1)]
    for i in range(1, len(points)):
        for j in range(len(points[i])):
            if points[i][j] != None: req_indices[i-1].append(j)
    
    num_req_indices = 0
    for point in req_indices: num_req_indices += len(point)-1
    
    # take initial jacobian (identity) and initial residual (all zeroes), which will generate first guess to be all zeroes
    # within newton loop:
        # take guess using current jacobian and check result
        # populate residual from the current guess
        # perturb all degrees of freedom in a for loop, and measure error for each perturbation (fill jacobian by measuring errors compared to the previous residual)
        # use jacobian to generate next best guess
    
    if num_req_indices < len(unknown_indices): raise ValueError(f"Too few boundary values given, needed {len(unknown_indices)}, got {num_req_indices}")
    if num_req_indices > len(unknown_indices): raise ValueError(f"Too many boundary values to fit, needed {len(unknown_indices)}, got {num_req_indices}")

    # Initial Jacobian
    jacobian = matrix.identity(num_req_indices)
    residual = matrix.zeroes(num_req_indices, 1)
    guess = matrix.zeroes(num_req_indices, 1)

    while True:
        # Take guess using jacobian (using linear equation solver)
        guess -= solve_sys(jacobian, residual)
        guessed_init_point = [guess[unknown_indices.index(i),0] if i in unknown_indices else init_point[i] for i in range(len(init_point))]
        vals = ode(order, method, eqn, indep_var_range, init_values = guessed_init_point, h = h, tolerance = tolerance)

        # Populate residual
        row = 0
        for i in range(len(req_indices)):
            indices = req_indices[i]
            pos = min(range(len(vals[0])), key = lambda x: abs(vals[0][x] - points[i+1][indices[0]]))
            for j in indices[1:]:
                residual[row, 0] = points[1+i][j] - vals[j][pos]
                row += 1
        
        # check if all required values are within tolerance
        if all(abs(residual[i,0]) < tolerance for i in range(num_req_indices)):
            print(guessed_init_point)
            return vals

        # Perturb all degrees of freedom
        perturbed_residual = matrix.zeroes(num_req_indices, 1)
        for i in range(len(unknown_indices)):
            perturbed_guess = guess.copy()
            perturbed_guess[i,0] += h
            perturbed_init_point = [perturbed_guess[unknown_indices.index(i),0] if i in unknown_indices else init_point[i] for i in range(len(init_point))]
            perturbed_vals = ode(order, method, eqn, indep_var_range, init_values = perturbed_init_point, h = h, tolerance = tolerance)

            row = 0
            for k in range(len(req_indices)):
                indices = req_indices[k]
                pos = min(range(len(vals[0])), key = lambda x: abs(vals[0][x] - points[k+1][indices[0]]))
                for j in indices[1:]:
                    perturbed_residual[row, 0] = points[1+k][j] - perturbed_vals[j][pos]
                    row += 1
                
                jacobian[k,i] = (perturbed_residual[k,0] - residual[k,0])/h

# General nth order ODEs, fixed step

def euler(func : callable, indep_var : float, dep_var : list, h : float, direction : bool = True, **kwargs) -> float :
    sign = 1 if direction else -1
    
    dydx = func(indep_var, dep_var)
    
    return [indep_var + sign*h] + [dep_var[i] + sign*h*dydx[i] for i in range(len(dydx))], h

def heun(func : callable, indep_var : float, dep_var : list, h : float, direction : bool = True, **kwargs) -> list :
    sign = 1 if direction else -1

    k1 = func(indep_var, dep_var)
    k2 = func(indep_var + (sign*h), [dep_var[i] + (sign*h)*k1[i] for i in range(len(k1))])

    return [indep_var + (sign*h)] + [dep_var[i] + (sign*h/2)*(k1[i] + k2[i]) for i in range(len(dep_var))], h

def rk2(func : callable, indep_var : float, dep_var : list, h : float, direction : bool = True, **kwargs) -> list :
    sign = 1 if direction else -1

    k1 = func(indep_var, dep_var)
    k2 = func(indep_var + (sign*h), [dep_var[i] + (sign*h)*k1[i] for i in range(len(k1))])

    return [indep_var + (sign*h)] + [dep_var[i] + (sign*h/2)*(k1[i] + k2[i]) for i in range(len(dep_var))], h

def rk4(func : callable, indep_var : float, dep_var : list, h : float, direction : bool = True, **kwargs) -> list :
    sign = 1 if direction else -1

    k1 = func(
        indep_var, 
        dep_var
        )
    k2 = func(
        indep_var + (sign*h/2),
        [dep_var[i] + (sign*h/2)*k1[i] for i in range(len(k1))]
        )
    k3 = func(
        indep_var + (sign*h/2),
        [dep_var[i] + (sign*h/2)*k2[i] for i in range(len(k2))]
        )
    k4 = func(
        indep_var + (sign*h),
        [dep_var[i] + (sign*h)*k3[i] for i in range(len(k3))]
        )

    return [indep_var + (sign*h)] + [dep_var[i] + (sign*h/3)*((k1[i]/2) + k2[i] + k3[i] + (k4[i]/2)) for i in range(len(dep_var))], h

# Adaptive step methods

def euler2(func : callable, indep_var : float, dep_var : list, h : float, tolerance : float = 0.001, direction : bool = True, **kwargs) -> list :
    sign = 1 if direction else -1
    h_new = h

    outvalues = []
    while True:
        step = sign*h_new

        k1 = func(indep_var, dep_var)
        y_euler = [dep_var[i] + step*k1[i] for i in range(len(k1))]
        k2 = func(indep_var + step, y_euler)
        y_rk2 = [dep_var[i] + (step/2)*(k1[i] + k2[i]) for i in range(len(dep_var))]

        error = max(abs(y_rk2[i] - y_euler[i]) for i in range(len(y_euler)))

        if error < tolerance :
            outvalues = [indep_var + step] + y_rk2
        
        if error == 0: scale = 3
        else: scale = 0.9*((tolerance/error)**0.5)

        h_new *= max(0.1, min(3, scale))

        if error < tolerance :
            return outvalues, h_new
