
import numpy
import sys
import random
import dadi
import imp

cur_temp = 0
number_of_steps = 4
T_total = 20000
dimention = 2
step_of_split = 3
old_value = 0
best_value = 0

T_weights = [1, 1.001, 0, 0]

def get_start_str(dim):
    result = "import dadi\ndef generated_model((n1,n2), pts):\n"
    if (dim == 1):
        result += "\txx = dadi.Numerics.default_grid(pts)\n"
    if (dim == 2):
        result += "\txx = yy = dadi.Numerics.default_grid(pts)\n"
    result += "\tphi = dadi.PhiManip.phi_1D(xx)\n\tm12, m21 = 0, 0\n"
    return result

def get_start_str_par(dim):
    result = "import dadi\ndef generated_model(params, (n1,n2), pts):\n"
    if (dim == 1):
        result += "\tpar1 = params[0]\n\txx = yy = dadi.Numerics.default_grid(pts)\n"
    if (dim == 2):
        result += "\tpar1, par2 = params\n\txx = yy = dadi.Numerics.default_grid(pts)\n"
    result += "\tphi = dadi.PhiManip.phi_1D(xx)\n\tm12, m21 = 0, 0\n"
    return result


def get_change_size_1d_str_par(T):
    return "\tphi = dadi.Integration.one_pop(phi, xx, " + str(T) + ", par1)\n"

def get_change_size_1d_str(T, nu):
    return "\tphi = dadi.Integration.one_pop(phi, xx, " + str(T) + ", " + str(nu) + ")\n"

def get_change_size_2d_str_par(T, migration_flag):
    if (migration_flag):
        m = random.random()
    else:
        m = 0
    strg = "\tT = " + str(T) + "\n"
    if (migration_flag):
        strg += "\tm12 = " + str(m) + "\n"
        strg += "\tm21 = " + str(m) + "\n"
    else:
        strg += "\tm12 = 0\n"
        strg += "\tm21 = 0\n"
    strg += "\tphi = dadi.Integration.two_pops(phi, xx,  T, par1, par2, m12, m21)\n"
    return strg

def get_change_size_2d_str(T, nu1, nu2, migration_flag):
    if (migration_flag):
        m = random.random()
    else:
        m = 0
    strg = "\tT = " + str(T) + "\n"
    strg += "\tnu1 = " + str(nu1) + "\n"
    strg += "\tnu2 = " + str(nu2) + "\n"
    if (migration_flag):
        strg += "\tm12 = " + str(m) + "\n"
        strg += "\tm21 = " + str(m) + "\n"
    else:
        strg += "\tm12 = 0\n"
        strg += "\tm21 = 0\n"
    strg += "\tphi = dadi.Integration.two_pops(phi, xx,  T, nu1, nu2, m12, m21)\n"
    return strg


def get_split_str():
    return "\tphi = dadi.PhiManip.phi_1D_to_2D(xx, phi)\n"

def get_finish_str(dim):
    return "\tsfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))\n\treturn sfs\n"

def generate_model_code(changes_id, before_par, after_par, migr_fl, filename, extra_data):
    global T_total
    global T_weights

    if (not len(extra_data) == 0):
        pos_var = extra_data[0]
    else:
        pos_var = -1

    par = -1
    sum_T = 0
    for i in xrange(number_of_steps):
        sum_T += T_weights[i]

    with open(filename, 'w') as output:
        if (not len(extra_data) == 0):
            output.write(get_start_str_par(1))
        else:
            output.write(get_start_str(2))
        for i in xrange(step_of_split):
            T = T_weights[i] * T_total / sum_T
            if (not i == pos_var):
                if (changes_id[i] == 0 and not T == 0):
                    output.write(get_change_size_1d_str(T, 1.0 / before_par[i]))
                if (changes_id[i] == 2 and not T == 0):
                    output.write(get_change_size_1d_str(T, before_par[i]))
                if (changes_id[i] == 1 and not T == 0):
                    output.write(get_change_size_1d_str(T, 1))
            else:
                output.write(get_change_size_1d_str_par(T))

        output.write(get_split_str())

        for j in xrange(len(after_par)):
            i = step_of_split + j - 1
            T = T_weights[i] * T_total / sum_T
            if (not i == pos_var):
                if (changes_id[i] == 0 and not T == 0):
                    output.write(get_change_size_2d_str(T, 1.0 / after_par[j][0], 1.0 / after_par[j][1], migr_fl[i]))
                if (changes_id[i] == 2 and not T == 0):
                    output.write(get_change_size_2d_str(T, after_par[j][0], after_par[j][1], migr_fl[i]))
                if (changes_id[i] == 1 and not T == 0):
                    output.write(get_change_size_2d_str(T, 1, 1, migr_fl[i]))
            else:
                output.write(get_change_size_2d_str_par(T, migr_fl[i]))

        output.write(get_finish_str(2))


def get_dimention(pos):
    if (pos < step_of_split):
        return 1
    return 2

def get_probs(pos, changes_id):
    incr = 0
    decr = 0
    for i in xrange(pos):
        if changes_id[i] == 0:
            decr += 1
        if changes_id[i] == 2:
            incr += 1
    if (incr == 0 or decr == 0): 
        return (0.5, 0.5)
    return (float(decr) / (decr+incr), float(incr) / (decr+incr))



def energy_func(ch, bf_par, af_par, mig_f):
    global popt_gl
    #print ch, popt_gl
    generate_model_code(ch, bf_par, af_par, mig_f, "demographic_models.py", [])

    import demographic_models

    imp.reload(dadi)
    imp.reload(demographic_models)


    dd = dadi.Misc.make_data_dict('fin_result_dadi.txt')
    data = dadi.Spectrum.from_data_dict(dd, ['Tan', 'Nam'], [5,7])
    ns = [5, 7]
    #data = dadi.Spectrum.from_file('YRI_CEU.fs')
    #ns = data.sample_sizes

    pts_l = [20, 40, 60]

    func = demographic_models.generated_model

#    upper_bound = [3, 3, 3, 3, 3, 3]
#    lower_bound = [1e-1, 1e-1, 1e-1, 0, 0, 0]

#    p0 = popt_gl
    func_ex = dadi.Numerics.make_extrap_log_func(func)

#    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
#                                  lower_bound=lower_bound)
#    popt_gl = dadi.Inference.optimize_log(data, func_ex, pts_l, 
#                                   lower_bound=lower_bound,
#                                   upper_bound=upper_bound,
#                                   verbose=len(p0), maxiter=100500)

    model = func_ex(ns, pts_l)
    return dadi.Inference.ll_multinom(model, data)

def get_param(ch, bf_par, af_par, mig_f, file_name, pos):
    global popt_gl
    #print ch, popt_gl
    print file_name
    generate_model_code(ch, bf_par, af_par, mig_f, file_name, [pos])

    import tmp
    imp.reload(dadi)
    imp.reload(tmp)

    dd = dadi.Misc.make_data_dict('fin_result_dadi.txt')
    data = dadi.Spectrum.from_data_dict(dd, ['Tan', 'Nam'], [5,7])
    ns = [5, 7]
    #data = dadi.Spectrum.from_file('YRI_CEU.fs')
    #ns = data.sample_sizes

    pts_l = [20, 40, 60]

    func = tmp.generated_model

#    upper_bound = [3, 3, 3, 3, 3, 3]
#    lower_bound = [1e-1, 1e-1, 1e-1, 0, 0, 0]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params([random.randint(1, 100500)], fold=1, upper_bound=[100500],
                                  lower_bound=[1])
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=[1],
                                   upper_bound=[100500],
                                   verbose=len(p0), maxiter=3)

    return popt


def check_next_state(changes_id, befor_split_params, after_split_params, migration_flags):
    global old_value
    global best_value
    global T_weights

    new_changes_id = list(changes_id)
    new_befor_split_params = list(befor_split_params)
    new_after_split_params = list(after_split_params)
    new_migration_flags = list(migration_flags)

    position = random.randint(0, number_of_steps - 1) #not split
    while (position == step_of_split):
        position = random.randint(0, number_of_steps - 1)

    dim = get_dimention(position)
    change = numpy.random.choice([0,1], 1, p = get_probs(position, changes_id))

    print "Change at " + str(position) + "pos from " + str(changes_id[position]) + " to " + str(new_changes_id[position])

#    if (change == 0):
#        new_changes_id[position] = 0
#    else:
#        new_changes_id[position] = 2

    if (position < step_of_split):
        new_befor_split_params[position] = get_param(new_changes_id, new_befor_split_params, new_after_split_params, new_migration_flags, "tmp.py", position) [0]
    else:
        new_after_split_params[position - step_of_split - 1] = (random.randint(1, 100), random.randint(1, 100))

    if (dim == 2):
        new_migration_flags[position] = random.choice([0,1])

    old_T_value = T_weights[position]
    T_weights[position] = random.randint(0, T_total)
    
    new_value = energy_func(new_changes_id, new_befor_split_params, new_after_split_params, new_migration_flags)
    delta = abs(old_value) - abs(new_value)
    if (delta > 0 or (numpy.exp(-float(delta) / cur_temp) > random.random())):
        print "Change at " + str(position) + "pos from " + str(changes_id[position]) + " to " + str(new_changes_id[position])
        old_value = new_value
        if (abs(new_value) < abs(best_value)):
            best_value = new_value
            print "write"
            generate_model_code(new_changes_id, new_befor_split_params, new_after_split_params, new_migration_flags, "best_model.py", [])

        print new_value
        return new_changes_id, new_befor_split_params, new_after_split_params, new_migration_flags

    T_weights[position] = old_T_value
    return changes_id, befor_split_params, after_split_params, migration_flags



def decrease_temp(temp, iteration):
    return (0.999 ** iteration) * float(temp)

def SA(initial_temp, final_temp, max_iter, changes_id, befor_split_params, after_split_params, migration_flags):
    global cur_temp
    cur_temp = initial_temp
    cur_iter = 0
    while(cur_iter < max_iter and cur_temp > final_temp):
        changes_id, befor_split_params, after_split_params, migration_flags = check_next_state(changes_id, befor_split_params, after_split_params, migration_flags)
        print changes_id
        #print changes_id
        cur_iter += 1
        cur_temp = decrease_temp(cur_temp, cur_iter)



migration_flags = [False] * number_of_steps
befor_split_params = [72381, 23710.8, 1] 
after_split_params = [] #[(100,100)] * (number_of_steps - step_of_split)
changes_id = [2,2,0,1] # 0 - decrease, 1 - nothing, 2 - increase

old_value = energy_func(changes_id, befor_split_params, after_split_params, migration_flags)
best_value = old_value
generate_model_code(changes_id, befor_split_params, after_split_params, migration_flags, "best_model.py", [])

SA(1000, 10, 100, changes_id, befor_split_params, after_split_params, migration_flags)

