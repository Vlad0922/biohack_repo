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
    func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics2D.split_mig)

    upper_bound = [500000, 500000, 20000, 0.001]
    lower_bound = [100, 100, 20, 0.001]

    # This is our initial guess for the parameters, which is somewhat arbitrary.
    p0 = [random.randint(a, b) for a, b in zip(lower_bound, upper_bound)[:-1]]
    while not validate_params(p0):
        p0 = [random.randint(a, b) for a, b in zip(lower_bound, upper_bound)[:-1]]
    p0.append(a + random.random()*(b-a))

    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                                  lower_bound=lower_bound)
    print('Beginning optimization ************************************************')
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound,
                                       verbose=len(p0), maxiter=3)
    # The verbose argument controls how often progress of the optimizer should be
    # printed. It's useful to keep track of optimization process.
    #print('Finshed optimization **************************************************')

    #print('Best-fit parameters: [{0}]'.format(', '.join('{:.1f}'.format(x) for x in popt)))

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


    with open("model_2_log.txt", "a") as log:
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



def validate_params(params):
    return True

def get_best(i, cnt=10):
    try:
        out_f = open("out_{}.txt".format(i), 'a')
        best = 100000000
        while True:
            queue.get_nowait()
            score, params = main(draw_img=False)
            if validate_params(params) and abs(score) < best:
                best = abs(score)
                best_params = params
                ss = '[{}]'.format(', '.join('{:.1f}'.format(x) for x in params))
                print('!!!!!Best-fit parameters score={0} params: {1}'.format(score, ss))
                out_f.write("{0}; {1}\n".format(score, best_params))
                out_f.flush()
            queue.task_done()
    except Queue.Empty:
        out_f.close()
        print("Queue is empty")
#
# for a in range(80):
#     queue.put(a)
#
# procs = []
# for i in range(8):
#     procs.append(multiprocessing.Process(target=get_best, args=(i,)))
#
# for proc in procs:
#     proc.start()
#
# for proc in procs:
#     proc.join()

main()