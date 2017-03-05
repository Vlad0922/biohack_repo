import Queue
import multiprocessing
import random

import dadi
import numpy

queue = Queue.Queue()

def main(draw_img=True):
    data = dadi.Spectrum.from_file('data.fs')
    ns = data.sample_sizes
    # These are the grid point settings will use for extrapolation.
    pts_l = [40,50,60]


    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(gen_model)

    upper_bound = [20000, 500000, 20000, 500000, 20000, 500000, ]
    lower_bound = [   10,   1000,    10,   1000,    10,   1000, ]

    # This is our initial guess for the parameters, which is somewhat arbitrary.
    p0 = [random.randint(a, b) for a, b in zip(lower_bound, upper_bound)]
    # while not validate_params(p0):
    #     p0 = [random.randint(a, b) for a, b in zip(lower_bound, upper_bound)]
    #
    # p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
    #                               lower_bound=lower_bound)
    # print('Beginning optimization ************************************************')
    # popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
    #                                    lower_bound=lower_bound,
    #                                    upper_bound=upper_bound,
    #                                    maxiter=3)
    # The verbose argument controls how often progress of the optimizer should be
    # printed. It's useful to keep track of optimization process.
    #print('Finshed optimization **************************************************')

    #print('Best-fit parameters: [{0}]'.format(', '.join('{:.1f}'.format(x) for x in popt)))
    popt = [18940.2, 138227.5, 17237.8, 226833.5, 17202.8, 3179.2]
    # Calculate the best-fit model AFS.
    model = func_ex(popt, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)
    #print('Maximum log composite likelihood: {0}'.format(ll_model))
    # The optimal value of theta given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    #print('Optimal value of theta: {0}'.format(theta))

    my_likehoodv = my_likehood(model, data)
    #print("my likehood={}".format(my_likehoodv))

    with open("model_4_log.txt", "a") as log:
        log.write("{0}; {1}; {2}\n".format(ll_model, my_likehoodv, popt))

    if draw_img:
        # Plot a comparison of the resulting fs with the data.
        import pylab
        pylab.figure(1)
        dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                            pop_ids =('Tan','Nam'))

        pylab.show()

    return ll_model, popt

def my_likehood(model, data):
    model2 = dadi.Inference.optimally_scaled_sfs(model, data)
    return numpy.sum((data - model2) ** 2)


def gen_model(params, ns, pts):
    T1, p1, T2, p2, T3, p3 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, p1)

    phi = dadi.Integration.one_pop(phi, xx, T2, p2)
    phi = dadi.Integration.one_pop(phi, xx, T3, p3)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs




def validate_params(params):
    return params[0] > params[2] > params[4] and params[1] < params[3] > params[5] and params[4] - params[2] < 15

main()