[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3667772.svg)](https://doi.org/10.5281/zenodo.3667772)
![Prose Linting](https://github.com/volodymyrss/integral-isgri-rate-meaning/workflows/Prose%20Linting/badge.svg)
![Executing Notebook](https://github.com/volodymyrss/integral-isgri-rate-meaning/workflows/Executing%20Notebook/badge.svg)

# Inter-OSA rate normalization for INTEGRAL ISGRI

The raw count rate produced in ISGRI by a given source strongly depends on time, and on the source position within the FoV.
In principle, appropriate conversion from the physical flux to the reconstructed source count rate can be expressed with the dispersive response provided with the source spectra.

However, the response model and the corrections are designed to allow the count rate to be approximately proportional to the source flux for different source positions and different times. This rate, reported in the standard OSA results, is a reconstructed rate aiming to reproduce the true rate of event detection that would be expected by a given source if 1) the source was on-axis, and 2) the various instrumental efficiency losses (mask support and the detector) did not exist. 

This "true" rate cannot be directly measured, and the meaning of the count rate depends on the instrument model assumed in the reconstruction process in OSA. The difference is especially important between OSA10.2 and OSA11.0, owing to a major progress in the detector understanding: 

OSA11.0 aims to be provide a reconstruction that relies as much as possible on the physical properties of the detector. And even though it would be possible to introduce additional factors to make the rate close to that of OSA10.2, this would inevitably introduce artificial features in the ISGRI spectra. This would be counterproductive given that one of the main goals of OSA11.0 was to actually correct such spectral features introduced by manipulations of the response files, that also made it necessary to have elaborate cross-normalization.

The information about the meaning of the ISGRI rate is contained in the response model expressed in the RMF and ARF structures. An example how to extract this normalization, for a given spectrum, is shown here:


```python
from oda_api.api import DispatcherAPI

disp=DispatcherAPI(host=host)

T1_2osa_utc="2016-01-01T00:00:00"
T2_2osa_utc="2017-01-01T00:00:00"

spec_data_osa10=disp.get_product(instrument='isgri',
                    product='isgri_spectrum',
                    T1=T1_2osa_utc,
                    T2=T2_2osa_utc,
                    query_type='Real',
                    osa_version='OSA10.2',
                    RA=ra,
                    DEC=dec,
                    product_type='Real',
                    selected_catalog=api_cat)

spec_data_osa11=disp.get_product(instrument='isgri',
                    product='isgri_spectrum',
                    T1=T1_2osa_utc,
                    T2=T2_2osa_utc,
                    query_type='Real',
                    osa_version='OSA11.0',
                    RA=ra,
                    DEC=dec,
                    product_type='Real',
                    selected_catalog=api_cat)
spec_data_osa10, spec_data_osa11

# Roques & Jourdain 2018
def crab_ph_cm2_s_kev(en):
    K=7.417e-4
    al1=-1.98
    al2=-2.33
    Ec=500.
    f=K*(en/100)**al1*(np.exp(-en/Ec))
    m=en>Ec*(al1-al2)
    f[m]=(K*((al1-al2)*Ec/100)**(al1-al2)*(en/100)**al2*np.exp(-(al1-al2)))[m]

    return f

def resp_norm(D, e1, e2, plot=False):
    arf=D._4.data_unit[1].data
    rmf_eb=D._5.get_data_unit_by_name('EBOUNDS').data
    rmf_mt=D._5.get_data_unit_by_name('SPECRESP MATRIX').data
    spec=D._3.data_unit[1]
    
    crab_ph_cm2_s_kev    
    ie1=arf['ENERG_LO']
    ie2=arf['ENERG_HI']
    
    source=crab_ph_cm2_s_kev(ie1)
    
    csource=np.outer(arf['SPECRESP']*source*(ie2-ie1),np.ones_like(rmf_eb['E_MIN']))*rmf_mt['MATRIX']
    
    if plot:
        plt.figure()
        plt.plot(
            rmf_eb['E_MIN'],
            csource.sum(0)/(rmf_eb['E_MAX']-rmf_eb['E_MIN'])
        )
        
        plt.plot(
            rmf_eb['E_MIN'],
            spec.data['RATE']/(rmf_eb['E_MAX']-rmf_eb['E_MIN'])
        )
        
        plt.loglog()

    rate_n = spec.data['RATE'][(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)].sum()

    n = csource[:,(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)].sum()
    print("response norm in", e1,e2,"is",n, "rate norm", rate_n)

    return n
    
resp_norm(spec_data_osa10, 30, 100, plot=True)
resp_norm(spec_data_osa11, 30, 100, plot=True)

```

See for a use case:
https://github.com/cdcihub/oda_api_benchmark/blob/master/examples/Crab_lc_longterm.ipynb
