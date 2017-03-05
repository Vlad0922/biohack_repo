import dadi
import numpy




def main():
    data = dadi.Spectrum.from_file('data.fs')
    ns = data.sample_sizes
    # These are the grid point settings will use for extrapolation.
    pts_l = [40,50,60]


    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(gen_model)

    upper_bound = [999, 248, 1400, 900]
    lower_bound = [697, 73, 900, 600]

    # This is our initial guess for the parameters, which is somewhat arbitrary.
    p0 = [100, 100, 1000, 700]
    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                                  lower_bound=lower_bound)
    print('Beginning optimization ************************************************')
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound,
                                       verbose=len(p0), maxiter=3)
    # The verbose argument controls how often progress of the optimizer should be
    # printed. It's useful to keep track of optimization process.
    print('Finshed optimization **************************************************')

    print('Best-fit parameters: {0}'.format(popt))

    # Calculate the best-fit model AFS.
    model = func_ex(popt, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)
    print('Maximum log composite likelihood: {0}'.format(ll_model))
    # The optimal value of theta given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    print('Optimal value of theta: {0}'.format(theta))

    my_likehoodv = my_likehood(model, data)
    print("my likehood={}".format(my_likehoodv))

    # Plot a comparison of the resulting fs with the data.
    import pylab
    pylab.figure(1)
    dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                        pop_ids =('Tan','Nam'))

    pylab.show()


def my_likehood(model, data):
    model2 = dadi.Inference.optimally_scaled_sfs(model, data)
    return numpy.sum((data - model2) ** 2)


def gen_model(params, ns, pts):
    nuB, nuF, TB, TF = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


main()
