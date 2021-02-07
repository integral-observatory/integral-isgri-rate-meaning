```python
from IPython.core.display import HTML
import requests
import re

import matplotlib.pylab as plt
from astropy import units as u
```


```python
# WebPIMMS request, deduced from web application requests

def webpimms_flux(e1: float, e2: float, slope: float=2):
    t = requests.post("https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl",
                  data={
                      "from": "INTEGRAL/ISGRI Count Rate",
                      "sat": "FLUX",
                      "range": f"{e1}-{e2}",
                      "etype": "kev",
                      "orange": f"{e1}-{e2}",
                      "otype": "kev",
                      "flusso": 1,
                      "nh": 0,
                      "red": "none",
                      "nhi": "none",
                      "model": "Power Law",
                      "gama": slope,
                      "solar": "1.0 Solar Abundance",
                  }
                 ).text      

    return float(re.search("PIMMS predicts a flux .*? of (.*?) ergs/cm/cm/s", t).groups()[0])

#HTML(t)


webpimms_flux(20, 40, 3)
```




    4.39e-11




```python
# query crab spectra and responses

from astroquery.simbad import Simbad
import numpy as np
from astropy.coordinates import SkyCoord

result_table = Simbad.query_object("Crab")

c = SkyCoord(result_table['RA'][0], result_table['DEC'][0], unit=("hourangle", "deg"))



from oda_api.api import DispatcherAPI

products_for_src = lambda o, name: { 
        v.meta_data['product']:v 
            for k,v in o.__dict__.items() 
            if getattr(v,'meta_data',{}).get('src_name') == name
    }



disp=DispatcherAPI(host="https://www.astro.unige.ch/cdci/astrooda/dispatch-data")

T1_2osa_utc="2016-01-01T00:00:00"
T2_2osa_utc="2017-01-01T00:00:00"


spec_data_osa10=disp.get_product(instrument='isgri',
                    product='isgri_spectrum',
                    T1=T1_2osa_utc,
                    T2=T2_2osa_utc,
                    query_type='Real',
                    osa_version='OSA10.2',
                    RA=c.ra.deg,
                    DEC=c.dec.deg,
                    product_type='Real',
                    #selected_catalog=api_cat
                    )

spec_data_osa11=disp.get_product(instrument='isgri',
                    product='isgri_spectrum',
                    T1=T1_2osa_utc,
                    T2=T2_2osa_utc,
                    query_type='Real',
                    osa_version='OSA11.0',
                    RA=c.ra.deg,
                    DEC=c.dec.deg,
                    product_type='Real',
                    #selected_catalog=api_cat,
                    )


```

    - waiting for remote response, please wait run_analysis https://www.astro.unige.ch/cdci/astrooda/dispatch-data
    T1 2016-01-01T00:00:00
    T2 2017-01-01T00:00:00
    query_type Real
    osa_version OSA10.2
    RA 83.63308333333332
    DEC 22.0145
    instrument isgri
    product_type isgri_spectrum
    off_line (False,)
    query_status ('new',)
    verbose (False,)
    session_id NKTHPIK113UMVHOV
    dry_run (False,)
    api True
    oda_api_version 1.0.2
    the job has been submitted on the remote server
     \ the job is working remotely, please wait status=done - job_id=-8431756989072857662  662 
    
    query done succesfully!
    - waiting for remote response, please wait run_analysis https://www.astro.unige.ch/cdci/astrooda/dispatch-data
    T1 2016-01-01T00:00:00
    T2 2017-01-01T00:00:00
    query_type Real
    osa_version OSA11.0
    RA 83.63308333333332
    DEC 22.0145
    instrument isgri
    product_type isgri_spectrum
    off_line (False,)
    query_status ('new',)
    verbose (False,)
    session_id XSN549KGTK0XWWAU
    dry_run (False,)
    api True
    oda_api_version 1.0.2
    the job has been submitted on the remote server
     \ the job is working remotely, please wait status=done - job_id=445217871065656659  659 
    
    query done succesfully!



```python
def extra_ticks(et):
    ax = plt.gca()
    tx = ax.get_xticks()
    ax.set_xticks(
        list(tx)+et
    )
    ax.set_xticklabels(
        list(tx)+et
    )

    
    

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

def powerlaw_ph_cm2_s_kev(en, slope):
    K=1e-5
    f=K*(en/100)**slope
    
    return f

def he_bb_ph_cm2_s_kev(en, T_keV):
    K=1e10
    f=K*(en/100)**3*np.exp(-en/T_keV)
    
    return f


def resp_norm(D, e1, e2, plot=False, model="crab"):
    p = products_for_src(D, 'Crab')
    
    arf=p['isgri_arf'].data_unit[1].data
    rmf_eb=p['isgri_rmf'].get_data_unit_by_name('EBOUNDS').data
    rmf_mt=p['isgri_rmf'].get_data_unit_by_name('SPECRESP MATRIX').data
    spec=p['isgri_spectrum'].data_unit[1]
    
    crab_ph_cm2_s_kev    
    ie1=arf['ENERG_LO']
    ie2=arf['ENERG_HI']
    
    if model == "crab":
        source=crab_ph_cm2_s_kev(ie1)
    else:
        source = model(ie1)
        
    renorm = np.outer(arf['SPECRESP']*source*(ie2-ie1),np.ones_like(rmf_eb['E_MIN']))*rmf_mt['MATRIX']    
    rate_n = np.nansum(spec.data['RATE'][(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)])
    csource=np.outer(arf['SPECRESP']*source*(ie2-ie1),np.ones_like(rmf_eb['E_MIN']))*rmf_mt['MATRIX']
    n = csource[:,(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)].sum()
    source *= rate_n/n
    
    m_flux = (ie1>e1) & (ie2<=e2)
    flux_erg_cm2_s_eband = np.nansum((source*(ie2-ie1)*ie1)[m_flux])*((1.*u.keV).to(u.erg).value)
    
    print("flux_eband", flux_erg_cm2_s_eband)
    
    
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
        
        extra_ticks(
            [20,25,30,35]
        )
        
        plt.xlim([15, 600])
        plt.grid()

    rate_n = np.nansum(spec.data['RATE'][(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)])

    n = csource[:,(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)].sum()
    
    erg_cm2_per_count = flux_erg_cm2_s_eband/n
    
    print("response norm in", e1,e2,"is",n, "rate norm", rate_n, "erg/cm2 per count", erg_cm2_per_count)

    return erg_cm2_per_count


erg_cm2_per_count = []
    
plot=True
for e1,e2 in [
    (20, 40),
    (40, 60),
    (60, 100),
]:
    for origin, c in [
        ('osa10_pl2', resp_norm(spec_data_osa10, e1, e2, plot=plot)),
        ('osa10_pl3.5', resp_norm(spec_data_osa10, e1, e2, plot=plot, model=lambda x: powerlaw_ph_cm2_s_kev(x, -2))),
        ('osa10_pl1.5', resp_norm(spec_data_osa10, e1, e2, plot=plot, model=lambda x: powerlaw_ph_cm2_s_kev(x, -1.5))),        
        ('osa10_bb1k', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: he_bb_ph_cm2_s_kev(x, 1))),  
        ('osa10_bb03k', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: he_bb_ph_cm2_s_kev(x, 0.3))),
        ('osa11_pl2', resp_norm(spec_data_osa11, e1, e2, plot=plot)),
        ('osa11_pl3.5', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: powerlaw_ph_cm2_s_kev(x, -2))),
        ('osa11_pl1.5', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: powerlaw_ph_cm2_s_kev(x, -1.5))),        
        ('osa11_bb1k', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: he_bb_ph_cm2_s_kev(x, 1))),  
        ('osa11_bb03k', resp_norm(spec_data_osa11, e1, e2, plot=plot, model=lambda x: he_bb_ph_cm2_s_kev(x, 0.3))),  
        ('pimms_pl2', webpimms_flux(e1,e2)),
        ('pimms_pl1.5', webpimms_flux(e1,e2,1.5)),
        ('pimms_pl3.5', webpimms_flux(e1,e2,3.5)),
    ]:    
        erg_cm2_per_count.append(dict(
            e1=e1,
            e2=e2,
            origin=origin,
            erg_cm2_in_count=c,
        ))
    
    plot=False

```

    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 7.090497497381998e-09
    response norm in 20 40 is 136.2734 rate norm 136.2734 erg/cm2 per count 5.203140990154163e-11
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 7.09258315555595e-09
    response norm in 20 40 is 136.27339 rate norm 136.2734 erg/cm2 per count 5.2046720682933876e-11
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 7.1830403660945766e-09
    response norm in 20 40 is 136.2734 rate norm 136.2734 erg/cm2 per count 5.271050695181583e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 1.3977833553164586e-08
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 6.663388742371043e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 9.705635738607391e-11
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 4.626784520806423e-13
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 8.318718676701798e-09
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 3.965625729513493e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 8.279860496591275e-09
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 3.9471016028015346e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 7.887348795981727e-09
    response norm in 20 40 is 209.77063 rate norm 209.77065 erg/cm2 per count 3.759987182375322e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 1.3977833553164586e-08
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 6.663388742371043e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 9.705635738607391e-11
    response norm in 20 40 is 209.77065 rate norm 209.77065 erg/cm2 per count 4.626784520806423e-13
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 3.6297541936603313e-09
    response norm in 40 60 is 37.257435 rate norm 37.25743 erg/cm2 per count 9.742362051396846e-11
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 3.6278243868664106e-09
    response norm in 40 60 is 37.25743 rate norm 37.25743 erg/cm2 per count 9.737183392807278e-11
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 3.6243265165552608e-09
    response norm in 40 60 is 37.257427 rate norm 37.25743 erg/cm2 per count 9.727796005839187e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 1.7336830639492434e-11
    response norm in 40 60 is 76.30243 rate norm 76.30242 erg/cm2 per count 2.2721204057904285e-13
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 0.0
    response norm in 40 60 is 76.30242 rate norm 76.30242 erg/cm2 per count 0.0
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.1020438473782535e-09
    response norm in 40 60 is 76.30243 rate norm 76.30242 erg/cm2 per count 5.376033096755265e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.089882855871672e-09
    response norm in 40 60 is 76.30242 rate norm 76.30242 erg/cm2 per count 5.3600957502102694e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.032036782627301e-09
    response norm in 40 60 is 76.30242 rate norm 76.30242 erg/cm2 per count 5.284284167754218e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 1.7336830639492434e-11
    response norm in 40 60 is 76.30243 rate norm 76.30242 erg/cm2 per count 2.2721204057904285e-13
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 0.0
    response norm in 40 60 is 76.30242 rate norm 76.30242 erg/cm2 per count 0.0
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 4.860451041592839e-09
    response norm in 60 100 is 32.83592 rate norm 32.835915 erg/cm2 per count 1.4802238751050807e-10
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 4.874101024960605e-09
    response norm in 60 100 is 32.83592 rate norm 32.835915 erg/cm2 per count 1.4843809031468926e-10
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    --> NAME PRIMARY
    --> NAME GROUPING
    --> NAME SPECRESP MATRIX
    --> NAME EBOUNDS
    flux_eband 4.934124815649976e-09
    response norm in 60 100 is 32.835915 rate norm 32.835915 erg/cm2 per count 1.5026609960407105e-10
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 3.0038440573418865e-12
    response norm in 60 100 is 52.64862 rate norm 52.64862 erg/cm2 per count 5.705456330663807e-14
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 0.0
    response norm in 60 100 is nan rate norm 52.64862 erg/cm2 per count nan
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.949327200999234e-09
    response norm in 60 100 is 52.648617 rate norm 52.64862 erg/cm2 per count 9.400678503422291e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.9538827979557455e-09
    response norm in 60 100 is 52.64862 rate norm 52.64862 erg/cm2 per count 9.409330654792449e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 4.979922966036023e-09
    response norm in 60 100 is 52.64862 rate norm 52.64862 erg/cm2 per count 9.458790959318591e-11
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 3.0038440573418865e-12
    response norm in 60 100 is 52.64862 rate norm 52.64862 erg/cm2 per count 5.705456330663807e-14
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    --> NAME PRIMARY
    --> NAME EBOUNDS
    --> NAME SPECRESP MATRIX
    flux_eband 0.0
    response norm in 60 100 is nan rate norm 52.64862 erg/cm2 per count nan


    <ipython-input-192-255fa5e76f96>:60: RuntimeWarning: divide by zero encountered in float_scalars
      source *= rate_n/n
    <ipython-input-192-255fa5e76f96>:60: RuntimeWarning: invalid value encountered in multiply
      source *= rate_n/n
    <ipython-input-192-255fa5e76f96>:68: RuntimeWarning: invalid value encountered in multiply
      csource=np.outer(arf['SPECRESP']*source*(ie2-ie1),np.ones_like(rmf_eb['E_MIN']))*rmf_mt['MATRIX']



    
![png](count-flux-conversion_files/count-flux-conversion_3_2.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_3.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_4.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_5.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_6.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_7.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_8.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_9.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_10.png)
    



    
![png](count-flux-conversion_files/count-flux-conversion_3_11.png)
    



```python
import json

json.dump(
    erg_cm2_per_count,
    open("conversions.json", "wt")
)

!cat conversions.json | jq
```

    [1;39m[
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.203140990154163e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.2046720682933876e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.271050695181583e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m6.663388742371043e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.626784520806423e-13[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m3.965625729513493e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m3.9471016028015346e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m3.759987182375322e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m6.663388742371043e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.626784520806423e-13[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.203e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.128e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m20[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.502e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.742362051396846e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.737183392807278e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.727796005839187e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m2.2721204057904285e-13[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m0[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.376033096755265e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.3600957502102694e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.284284167754218e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m2.2721204057904285e-13[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m0[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.747e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.766e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m40[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m4.691e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m1.4802238751050807e-10[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m1.4843809031468926e-10[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m1.5026609960407105e-10[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.705456330663807e-14[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa10_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39mnull[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.400678503422291e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.409330654792449e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.458790959318591e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb1k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m5.705456330663807e-14[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"osa11_bb03k"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39mnull[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl2"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m8.813e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl1.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m9.01e-11[0m[1;39m
      [1;39m}[0m[1;39m,
      [1;39m{
        [0m[34;1m"e1"[0m[1;39m: [0m[0;39m60[0m[1;39m,
        [0m[34;1m"e2"[0m[1;39m: [0m[0;39m100[0m[1;39m,
        [0m[34;1m"origin"[0m[1;39m: [0m[0;32m"pimms_pl3.5"[0m[1;39m,
        [0m[34;1m"erg_cm2_in_count"[0m[1;39m: [0m[0;39m8.272e-11[0m[1;39m
      [1;39m}[0m[1;39m
    [1;39m][0m



```python
import pandas as pd

rdiff_pc = lambda x,ax: (100*(ax-x)/ax)

print(f"\033[34m     {'key':20s} {'value':10s} ", end="")
for ref in 'osa10_pl2', 'osa11_pl2', 'pimms_pl2':
    print(f" \033[34m{'vs-'+ref:15s}\033[0m", end="")
print()

for (e1,e2), v in pd.DataFrame(erg_cm2_per_count).groupby(['e1', 'e2']):
    bo = dict(list(v.groupby('origin')))
    bo = {**bo, **hilight_rdata}
    print(f"\033[32m{e1:3d} {e2:3d}\033[0m") 
    for k, v in bo.items():
        try:
            c = v.erg_cm2_in_count.iloc[0]
        except:
            c = v[(e1,e2)]

        print(f"    {k:20s} {c:9.2e}", end="")
            
        for ref in 'osa10_pl2', 'osa11_pl2', 'pimms_pl2':
            dc = rdiff_pc(c, bo[ref].erg_cm2_in_count.iloc[0])

            if abs(dc)>50:
                C = '\033[31m'        
            elif abs(dc)>30:
                C = '\033[33m'
            else:
                C = '\033[32m'
        
            print(f"{C}{dc:15.0f}%\033[0m", end="")
        print(f"")
          
```

    [34m     key                  value       [34mvs-osa10_pl2   [0m [34mvs-osa11_pl2   [0m [34mvs-pimms_pl2   [0m
    [32m 20  40[0m
        osa10_bb03k           4.63e-13[31m             99%[0m[31m             99%[0m[31m             99%[0m
        osa10_bb1k            6.66e-11[32m            -28%[0m[31m            -68%[0m[31m            -59%[0m
        osa10_pl1.5           5.27e-11[32m             -1%[0m[33m            -33%[0m[32m            -25%[0m
        osa10_pl2             5.20e-11[32m              0%[0m[33m            -31%[0m[32m            -24%[0m
        osa10_pl3.5           5.20e-11[32m             -0%[0m[33m            -31%[0m[32m            -24%[0m
        osa11_bb03k           4.63e-13[31m             99%[0m[31m             99%[0m[31m             99%[0m
        osa11_bb1k            6.66e-11[32m            -28%[0m[31m            -68%[0m[31m            -59%[0m
        osa11_pl1.5           3.76e-11[32m             28%[0m[32m              5%[0m[32m             11%[0m
        osa11_pl2             3.97e-11[32m             24%[0m[32m              0%[0m[32m              6%[0m
        osa11_pl3.5           3.95e-11[32m             24%[0m[32m              0%[0m[32m              6%[0m
        pimms_pl1.5           4.13e-11[32m             21%[0m[32m             -4%[0m[32m              2%[0m
        pimms_pl2             4.20e-11[32m             19%[0m[32m             -6%[0m[32m              0%[0m
        pimms_pl3.5           4.50e-11[32m             13%[0m[32m            -14%[0m[32m             -7%[0m
        hili_plaw_1.5         5.96e-11[32m            -15%[0m[31m            -50%[0m[33m            -42%[0m
        hili_plaw_ 2          1.84e-11[31m             65%[0m[31m             54%[0m[31m             56%[0m
        hili_plaw_2.5         2.58e-11[31m             50%[0m[33m             35%[0m[33m             39%[0m
        hili_plaw_ 3          2.02e-11[31m             61%[0m[33m             49%[0m[31m             52%[0m
        hili_plaw_3.5         1.11e-11[31m             79%[0m[31m             72%[0m[31m             73%[0m
        hili_bbody_0.06       1.84e-11[31m             65%[0m[31m             54%[0m[31m             56%[0m
        hili_bbody_0.1        1.84e-11[31m             65%[0m[31m             54%[0m[31m             56%[0m
        hili_bbody_0.3        1.84e-11[31m             65%[0m[31m             54%[0m[31m             56%[0m
        hili_bbody_ 1         1.84e-11[31m             65%[0m[31m             54%[0m[31m             56%[0m
    [32m 40  60[0m
        osa10_bb03k           0.00e+00[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa10_bb1k            2.27e-13[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa10_pl1.5           9.73e-11[32m              0%[0m[31m            -81%[0m[31m           -105%[0m
        osa10_pl2             9.74e-11[32m              0%[0m[31m            -81%[0m[31m           -105%[0m
        osa10_pl3.5           9.74e-11[32m              0%[0m[31m            -81%[0m[31m           -105%[0m
        osa11_bb03k           0.00e+00[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa11_bb1k            2.27e-13[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa11_pl1.5           5.28e-11[33m             46%[0m[32m              2%[0m[32m            -11%[0m
        osa11_pl2             5.38e-11[33m             45%[0m[32m              0%[0m[32m            -13%[0m
        osa11_pl3.5           5.36e-11[33m             45%[0m[32m              0%[0m[32m            -13%[0m
        pimms_pl1.5           4.77e-11[31m             51%[0m[32m             11%[0m[32m             -0%[0m
        pimms_pl2             4.75e-11[31m             51%[0m[32m             12%[0m[32m              0%[0m
        pimms_pl3.5           4.69e-11[31m             52%[0m[32m             13%[0m[32m              1%[0m
        hili_plaw_1.5         1.41e-11[31m             86%[0m[31m             74%[0m[31m             70%[0m
        hili_plaw_ 2          1.29e-11[31m             87%[0m[31m             76%[0m[31m             73%[0m
        hili_plaw_2.5         1.07e-11[31m             89%[0m[31m             80%[0m[31m             78%[0m
        hili_plaw_ 3          8.36e-12[31m             91%[0m[31m             84%[0m[31m             82%[0m
        hili_plaw_3.5         3.73e-12[31m             96%[0m[31m             93%[0m[31m             92%[0m
        hili_bbody_0.06       1.29e-11[31m             87%[0m[31m             76%[0m[31m             73%[0m
        hili_bbody_0.1        1.29e-11[31m             87%[0m[31m             76%[0m[31m             73%[0m
        hili_bbody_0.3        1.29e-11[31m             87%[0m[31m             76%[0m[31m             73%[0m
        hili_bbody_ 1         1.29e-11[31m             87%[0m[31m             76%[0m[31m             73%[0m
    [32m 60 100[0m
        osa10_bb03k                nan[32m            nan%[0m[32m            nan%[0m[32m            nan%[0m
        osa10_bb1k            5.71e-14[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa10_pl1.5           1.50e-10[32m             -2%[0m[31m            -60%[0m[31m            -71%[0m
        osa10_pl2             1.48e-10[32m              0%[0m[31m            -57%[0m[31m            -68%[0m
        osa10_pl3.5           1.48e-10[32m             -0%[0m[31m            -58%[0m[31m            -68%[0m
        osa11_bb03k                nan[32m            nan%[0m[32m            nan%[0m[32m            nan%[0m
        osa11_bb1k            5.71e-14[31m            100%[0m[31m            100%[0m[31m            100%[0m
        osa11_pl1.5           9.46e-11[33m             36%[0m[32m             -1%[0m[32m             -7%[0m
        osa11_pl2             9.40e-11[33m             36%[0m[32m              0%[0m[32m             -7%[0m
        osa11_pl3.5           9.41e-11[33m             36%[0m[32m             -0%[0m[32m             -7%[0m
        pimms_pl1.5           9.01e-11[33m             39%[0m[32m              4%[0m[32m             -2%[0m
        pimms_pl2             8.81e-11[33m             40%[0m[32m              6%[0m[32m              0%[0m
        pimms_pl3.5           8.27e-11[33m             44%[0m[32m             12%[0m[32m              6%[0m
        hili_plaw_1.5         2.24e-11[31m             85%[0m[31m             76%[0m[31m             75%[0m
        hili_plaw_ 2          1.62e-11[31m             89%[0m[31m             83%[0m[31m             82%[0m
        hili_plaw_2.5         1.12e-11[31m             92%[0m[31m             88%[0m[31m             87%[0m
        hili_plaw_ 3          6.34e-12[31m             96%[0m[31m             93%[0m[31m             93%[0m
        hili_plaw_3.5         3.61e-12[31m             98%[0m[31m             96%[0m[31m             96%[0m
        hili_bbody_0.06       1.62e-11[31m             89%[0m[31m             83%[0m[31m             82%[0m
        hili_bbody_0.1        1.62e-11[31m             89%[0m[31m             83%[0m[31m             82%[0m
        hili_bbody_0.3        1.62e-11[31m             89%[0m[31m             83%[0m[31m             82%[0m
        hili_bbody_ 1         1.62e-11[31m             89%[0m[31m             83%[0m[31m             82%[0m



```python

#
# === INTEGRAL ===========================================
#
# Only calculated for plaw and no absorption.
# (Mission, Model, Spectral Index) : (20-40keV, 40-60keV, 60-100keV)

hilight_eband_index=((20,40), (40,60), (60,100))

hilight_data={
    # ('Integral','plaw',0.5): MISSING
    # ('Integral','plaw',1.0): MISSING
    ('Integral','plaw',1.5):(5.959E-11,1.411E-11,2.237E-11),
    # ('Integral','plaw',1.7): MISSING
    ('Integral','plaw',2.0):(1.844E-11,1.290E-11,1.624E-11),
    ('Integral','plaw',2.5):(2.583E-11,1.068E-11,1.116E-11),
    ('Integral','plaw',3.0):(2.024E-11,0.836E-11,0.634E-11),
    ('Integral','plaw',3.5):(1.115E-11,0.373E-11,0.361E-11),

    ('Integral','bbody',0.06):(1.844E-11,1.290E-11,1.624E-11),
    ('Integral','bbody',0.1):(1.844E-11,1.290E-11,1.624E-11),
    ('Integral','bbody',0.3):(1.844E-11,1.290E-11,1.624E-11),
    ('Integral','bbody',1.0):(1.844E-11,1.290E-11,1.624E-11),
}

hilight_rdata = {}

for ik, v in hilight_data.items():
    D={}
    hilight_rdata[f"hili_{ik[1]}_{ik[2]:2g}"] = D
    for (e1,e2), v in zip(hilight_eband_index, hilight_data[ik]):
        print(e1,e2,ik,v)
        D[(e1,e2)]=v
```

    20 40 ('Integral', 'plaw', 1.5) 5.959e-11
    40 60 ('Integral', 'plaw', 1.5) 1.411e-11
    60 100 ('Integral', 'plaw', 1.5) 2.237e-11
    20 40 ('Integral', 'plaw', 2.0) 1.844e-11
    40 60 ('Integral', 'plaw', 2.0) 1.29e-11
    60 100 ('Integral', 'plaw', 2.0) 1.624e-11
    20 40 ('Integral', 'plaw', 2.5) 2.583e-11
    40 60 ('Integral', 'plaw', 2.5) 1.068e-11
    60 100 ('Integral', 'plaw', 2.5) 1.116e-11
    20 40 ('Integral', 'plaw', 3.0) 2.024e-11
    40 60 ('Integral', 'plaw', 3.0) 8.36e-12
    60 100 ('Integral', 'plaw', 3.0) 6.34e-12
    20 40 ('Integral', 'plaw', 3.5) 1.115e-11
    40 60 ('Integral', 'plaw', 3.5) 3.73e-12
    60 100 ('Integral', 'plaw', 3.5) 3.61e-12
    20 40 ('Integral', 'bbody', 0.06) 1.844e-11
    40 60 ('Integral', 'bbody', 0.06) 1.29e-11
    60 100 ('Integral', 'bbody', 0.06) 1.624e-11
    20 40 ('Integral', 'bbody', 0.1) 1.844e-11
    40 60 ('Integral', 'bbody', 0.1) 1.29e-11
    60 100 ('Integral', 'bbody', 0.1) 1.624e-11
    20 40 ('Integral', 'bbody', 0.3) 1.844e-11
    40 60 ('Integral', 'bbody', 0.3) 1.29e-11
    60 100 ('Integral', 'bbody', 0.3) 1.624e-11
    20 40 ('Integral', 'bbody', 1.0) 1.844e-11
    40 60 ('Integral', 'bbody', 1.0) 1.29e-11
    60 100 ('Integral', 'bbody', 1.0) 1.624e-11

