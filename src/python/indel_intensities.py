from common import *

def plot(in_fn, out_fn, xlabel="Intensity"):
    data = read_table(in_fn)
   
    approx_density = (data['count'] / 0.01) / np.sum(data['count'])

    intensities = data['intensity'] / 100.0

    plt.plot(intensities, approx_density, 'b-')
    plt.fill_between(intensities, approx_density, color="grey", alpha=0.8)
    plt.xticks(np.arange(0.5, np.max(intensities)))
    plt.xlabel(xlabel)
    plt.ylabel("Probability density")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.savefig(out_fn)
