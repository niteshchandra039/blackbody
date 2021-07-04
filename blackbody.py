import numpy as np
from astropy.constants import astropyconst20 as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rcParams

a = 2*const.h*((const.c**2).to(u.nm**2 / u.s**2))

b = ((const.h*const.c)/(const.k_B)).to(u.K*u.nm)

def planck(lam, T):
    return ( (a/lam**5)*(1/(np.exp(b/(lam*T))-1)) )

wave = np.arange(100,1400,10)


spectrum_1 = (planck(wave*u.nm, 4000*u.K)).to(u.W/(u.m**(2)*u.nm))
spectrum_2 = (planck(wave*u.nm, 5777*u.K)).to(u.W/(u.m**(2)*u.nm))
spectrum_3 = (planck(wave*u.nm, 7000*u.K)).to(u.W/(u.m**(2)*u.nm))




fig, ax = plt.subplots(figsize=(8,6))
rcParams.update({'font.size': 12})


ax.plot(wave,spectrum_3/10**4, c='#8b00ff', label=r'T = 7000 K')
ax.plot(wave,spectrum_2/10**4, c='#00ff00', label=r'T = 5777 K')
ax.plot(wave,spectrum_1/10**4, c='#ff0000', label=r'T = 4000 K')





ax.text(500, 1.5, 'Visible', color='k', fontdict={'fontsize': 20})
ax.annotate('', (400, 1.3), (750, 1.3), arrowprops={'arrowstyle': '<|-|>',
                                                    'color': 'w', 'lw': 2})


# Finally, add some colourful rectangles representing a rainbow in the
# visible region of the spectrum.
# Dictionary mapping of wavelength regions (nm) to approximate RGB values
rainbow_rgb = { (400, 440): '#8b00ff', (440, 460): '#4b0082',
                (460, 500): '#0000ff', (500, 570): '#00ff00',
                (570, 590): '#ffff00', (590, 620): '#ff7f00',
                (620, 750): '#ff0000'}
for wv_range, rgb in rainbow_rgb.items():
    plt.axvspan(*wv_range, color=rgb, ec='none', alpha=1)
    
    
ax.tick_params(direction='in', which='both')
ax.minorticks_on()

ax.set_xlabel(r'Wavelength($\lambda$) (nm)')
ax.set_ylabel(r'$B_\lambda(T)$ $( 10^{4}$ $ W m^{-2} nm^{-1} sr^{-1})$')
plt.legend()

plt.savefig('BlackBody.png', dpi=200, facecolor='w', edgecolor='w',
          orientation='portrait', papertype=None, format=None,
          transparent=False, bbox_inches='tight', pad_inches=0.1, metadata=None)

plt.show()
