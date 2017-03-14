# -- external imports --
import numpy as np
from onezone import imf as _imf # to do IMF sampling
import time
from fastloop import loop_over_rnum

#
#
#
np.random.seed(12345)
_imf.np.random.seed(12345)

# make IMF function a global for now
# and set to assumed defaults --- make this 
# easily modifiable in the future
_IMF = _imf.salpeter()
_IMF.M_min =   1.0
_IMF.M_max = 100.0

def _check_scalar_input(x):

    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None]
        scalar_input = True

    return x, scalar_input


def _check_scalar_output(x, scalar_input):

    if scalar_input:
        return np.squeeze(x)
    else:
        return x


def exponential_disk(r, z, a, b):
    """
    An un-normalized Miyamoto-Nagai density
    function. The full density function can be computed
    by multiplying by (b^2 * M)/(4*pi), where M is the 
    total disk mass. This is not needed for establishing
    the IC's since it gets normalized out in computing the probability
    distribution anyway.
    """

    # allow for computing at many r/z values with only one z/r value
    if np.size(r) != np.size(z):
        if not (np.size(r) == 1 or np.size(z) == 1):
            return "r and z must be same size OR one (or both) must be scalar"
        if np.size(z) > np.size(r):
            r = np.ones(np.size(z)) * r
        else:
            z = np.ones(np.size(r)) * z

    zs_bs = np.sqrt(z*z + b*b)
    numerator   = a*r*r + (a + 3.0 * zs_bs)*(a + zs_bs)
    denominator = (r*r +  (a + zs_bs)**2)**(2.5) * (zs_bs)**3

    return numerator / denominator

def _find_bin(x, array):
    """
    Wrapper on argmin to ensure selected bin
    is less than desired value. I.e:

        array[index] <= x < array[index+1]

    where index is the found and returned bin index.
    array must be sorted in increasing order.
    """

    if x < array[0] or x > array[-1]:
        print "search value out of bounds on array"
        print x, array[0], array[-1]

        return -1

    index = np.abs(array - x).argmin()

    if x < array[index] and index > 0:
        index = index - 1
    if x < array[index] and index > 0:
        index = index - 1
    if x < array[index] and index > 0:
        print i+2, array[i+2], x
        print i, array[i], x
        print "Failure finding bin"

    return index


class particleIC(object):

    def __init__(self, M_star, a, b, rmax, zmax,
                       Z = 0.0043, rmin = 0.1, zmin = 0.1,
                       IMF = _IMF, npoints = 1000):

        #
        # Set properties
        #
        self._M_star  = M_star
        self._a       = a
        self._b       = b
        self._rmin    = rmin
        self._zmin    = zmin
        self._rmax    = rmax
        self._zmax    = zmax
        self._IMF     = IMF
        self._npoints = npoints

        self._metallicity = Z

        return

    def generate(self):
        """
        Given initial condition properties, generate a particle
        distribution. Particles are random sampled from an IMF
        and deposited over an exponential disk profile with 
        specified scale radius, height, and cutoff radius.
        """

        if not self._is_dist_tabulated():
            self._tabulate_probability_distribution()
        

        # sample the IMF - this sets the number of stars we need
        self.M = self._IMF.sample(M = self._M_star)
        self.number_of_particles = np.size(self.M)

        # initialize position arrays (r,theta,z)
        self.r = np.zeros(self.number_of_particles)
        self.z = np.zeros(self.number_of_particles)
        self.theta = np.zeros(self.number_of_particles)

        # get random numbers for radial position and set r
        delta_r = self._tabulated_prob_r[-1] - self._tabulated_prob_r[0]
        rnum   = (np.random.rand(self.number_of_particles) * delta_r) + self._tabulated_prob_r[0]
        for i in np.arange(self.number_of_particles):
            bin = _find_bin(rnum[i], self._tabulated_prob_r)
            self.r[i] = 10.0**(bin*self._tabulated_dr + self._r_o)

        # get random numbers for vertical position and set z
        delta_z = self._tabulated_prob_z[-1] - self._tabulated_prob_z[0]
        rnum   = np.random.rand(self.number_of_particles) * delta_z + self._tabulated_prob_z[0]
        for i in np.arange(self.number_of_particles):
            bin = _find_bin(rnum[i], self._tabulated_prob_z)
            self.z[i] = 10.0**(bin*self._tabulated_dz + self._z_o)

        # now randomly set z's to positive or negative
        rnum = np.random.rand(self.number_of_particles)
        ones = np.ones(self.number_of_particles)
        ones[rnum < 0.5] = -1
        self.z = self.z * ones

        # and set angle (easy)
        rnum = np.random.rand(self.number_of_particles)
        self.theta = 2.0 * np.pi * rnum

        self.metallicity = np.ones(self.number_of_particles)*self._metallicity        
        # done
        start = time.time()
        self.write_IC()
        end = time.time()
        print "write out took ", end - start
       
        return

    def _tabulate_probability_distribution(self):

        dr = np.log10(self._rmax / self._rmin) / (1.0*(self._npoints-1))
        dz = np.log10(self._zmax / self._zmin) / (1.0*(self._npoints-1))

        r_o = np.log10(self._rmin)
        z_o = np.log10(self._zmin)

        r   = 10.0**(r_o + np.arange(0,self._npoints)*dr)
        z   = 10.0**(z_o + np.arange(0,self._npoints)*dz)

        # sample exponential disk at z = 0
        r_dist = np.cumsum( exponential_disk(r, 0.0, self._a, self._b))
        r_dist = r_dist / (r_dist[-1])

        # sample exponential disk at r = 0
        z_dist = np.cumsum( exponential_disk(0.0, z, self._a, self._b))
        z_dist = z_dist / (z_dist[-1])

        # save tabulated properties
        self._r_o = r_o
        self._z_o = z_o
        self._tabulated_dr = dr
        self._tabulated_dz = dz
        self._tabulated_r  = r
        self._tabulated_z  = z
        self._tabulated_prob_r = r_dist
        self._tabulated_prob_z = z_dist

        return

    def _is_dist_tabulated(self):
        
        if hasattr(self, '_tabulated_r'):
            return (np.size(self._tabulated_r) == self._npoints)
        else:
            return False

    def write_IC(self, outfile = './particle_IC.in'):

        with open(outfile, 'w') as f:
            header = "# M Z x y z\n"
            fmt = "%.3f %3.3E %5.5E %5.5E %5.5E\n"

            f.write(header)

            x = self.x
            y = self.y

            for i in np.arange(self.number_of_particles):
             
                f.write(fmt%(self.M[i], self.metallicity[i],\
                             x[i], y[i], self.z[i]))


        print "wrote IC's for %i particles to "%(self.number_of_particles) + outfile

        return

    @property
    def x(self):
        return self.r * np.cos(self.theta)

    @property
    def y(self):
        return self.r * np.sin(self.theta)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # perform a test of the IC's and plot
    SFR = 1.0E-4
    dt  = 25.0E6

 
    M_star = SFR*dt
    a      = 250.0 # parsec
    b      = 125.0 # parsec

    rmax   = a*2.0
    zmax   = b*2.0

    start = time.time()
    IC = particleIC(M_star, a, b, rmax, zmax)
    IC.generate()
    end = time.time()
    print "generation took ", end - start

    # now plot these
    fig, ax = plt.subplots(1,3)

    ax[0].scatter(IC.x, IC.z, color='black', s = IC.M)
    ax[1].scatter(IC.y, IC.z, color='black', s = IC.M)
    ax[2].scatter(IC.x, IC.y, color='black', s = IC.M)

    ax[0].set_xlim(-rmax,rmax)
    ax[0].set_ylim(-rmax,rmax)
#    ax[0].set_ylim(-b,b)
    ax[1].set_xlim(-rmax,rmax)
    ax[1].set_ylim(-rmax,rmax)
#    ax[1].set_ylim(-b,b)
    ax[2].set_xlim(-rmax,rmax)
    ax[2].set_ylim(-rmax,rmax)

    ax[0].set_xlabel("x (pc)")
    ax[0].set_ylabel("z (pc)")
    ax[1].set_xlabel("y (pc)")
    ax[1].set_ylabel("z (pc)")
    ax[2].set_xlabel("x (pc)")
    ax[2].set_ylabel("y (pc)")

    fig.set_size_inches(18,6)
    plt.tight_layout()

    plt.savefig("particle_IC_test.png")
    plt.close(fig)

    # now plot the surface density profile
    rmin = 0.0
    rmax = rmax
    dr   = 50.0

    rbin = np.arange(rmin*rmin, (rmax)**2 + dr**2, dr**2)
    rbin = np.sqrt(rbin)
    centers = 0.5 * ( rbin[1:] + rbin[:-1])
    Mtot = np.zeros(np.size(rbin)-1)
    SD   = np.zeros(np.size(rbin)-1)
 
    for i in np.arange(1, np.size(rbin)-1):
        Mtot[i] = np.sum(  IC.M[ (IC.r < rbin[i]) * (IC.r >= rbin[i-1])])
        SD[i]   = Mtot[i] / (np.pi * ( rbin[i]**2 - rbin[i-1]**2))

    fig, ax = plt.subplots(1,2)

    x = rbin[:-1]
    ax[0].plot(x, Mtot,            label = 'mass',       color = 'black', ls = '--', lw = 3, drawstyle='steps-pre')
    ax[0].plot(x, np.cumsum(Mtot), label = 'cumulative', color = 'black', ls = '-', lw = 3, drawstyle='steps-pre')
    ax[0].set_xlim(rmin, rmax)
    ax[0].semilogy()
    ax[0].set_ylabel(r'M (M$_{\odot}$)')
    ax[0].set_xlabel(r'r (pc)')
    ax[0].legend(loc='best')

    ax[1].plot(x, SD, label = 'Surface Density', color = 'black', lw = 3, drawstyle='steps-pre')

#    rsqr = 1.0 / (centers * centers)
#    rsqr = (np.max(SD) / np.max(rsqr)) * rsqr
#    ax[1].plot(centers, rsqr, lw = 3, color = 'black', ls = '--')


    ax[1].set_xlim(rmin, rmax)
    ax[1].semilogy()
    ax[1].set_xlabel(r'r (pc)')
    ax[1].set_ylabel(r'$\Sigma$ (M$_{\odot}$ pc$^{-2}$')

    fig.set_size_inches(12,6)
    plt.tight_layout()
 
    plt.savefig("particle_IC_profile.png")
    plt.close(fig)

    print np.min(IC.M), np.max(IC.M)
