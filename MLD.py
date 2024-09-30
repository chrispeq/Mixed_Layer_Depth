#module file : MLD.py



def zmld_boyer(s, t, p):
    """
    Computes mixed layer depth, based on de Boyer Montégut et al., 2004.
    Parameters
    Downloaded from: https://github.com/pyoceans/python-oceans/blob/master/oceans/sw_extras/sw_extras.py
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [deg(ITS-90)]
    p : array_like
        pressure [db].
    Notes
    -----
    Based on density with fixed threshold criteria
    de Boyer Montégut et al., 2004. Mixed layer depth over the global ocean:
        An examination of profile data and a profile-based climatology.
        doi:10.1029/2004JC002378
    dataset for test and more explanation can be found at:
    http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php
    Codes based on : http://mixedlayer.ucsd.edu/
    """
    import gsw
    import numpy as np
    import numpy.ma as ma
    m = len(s)
    # starti = min(find((pres-10).^2==min((pres-10).^2)));
    starti = np.min(np.where(((p - 10.)**2 == np.min((p - 10.)**2)))[0])
    pres = p[starti:m]
    sal = s[starti:m]
    temp = t[starti:m]
    starti = 0
    m = len(sal)
    pden = gsw.rho(sal, temp,0)-1000

    mldepthdens_mldindex = m
    for i, ppp in enumerate(pden):
        if np.abs(pden[starti] - ppp) > .03:
            mldepthdens_mldindex = i
            break

    # Interpolate to exactly match the potential density threshold.
    presseg = [pres[mldepthdens_mldindex-1], pres[mldepthdens_mldindex]]
    pdenseg = [pden[starti] - pden[mldepthdens_mldindex-1], pden[starti] -
               pden[mldepthdens_mldindex]]
    P = np.polyfit(presseg, pdenseg, 1)
    presinterp = np.linspace(presseg[0], presseg[1], 3)
    pdenthreshold = np.polyval(P, presinterp)

    # The potential density threshold MLD value:
    ix = np.max(np.where(np.abs(pdenthreshold) < 0.03)[0])
    mldepthdens_mldindex = presinterp[ix]

    # Search for the first level that exceeds the temperature threshold.
    mldepthptmp_mldindex = m
    for i, tt in enumerate(temp):
        if np.abs(temp[starti] - tt) > 0.3:
            mldepthptmp_mldindex = i
            break

    # Interpolate to exactly match the temperature threshold.
    presseg = [pres[mldepthptmp_mldindex-1], pres[mldepthptmp_mldindex]]
    tempseg = [temp[starti] - temp[mldepthptmp_mldindex-1],
               temp[starti] - temp[mldepthptmp_mldindex]]
    P = np.polyfit(presseg, tempseg, 1)
    presinterp = np.linspace(presseg[0], presseg[1], 3)
    tempthreshold = np.polyval(P, presinterp)

    # The temperature threshold MLD value:
    ix = np.max(np.where(np.abs(tempthreshold) < 0.3)[0])
    mldepthptemp_mldindex = presinterp[ix]

    return mldepthdens_mldindex, mldepthptemp_mldindex



def zmld_boyervar(s, t, p,denscriteria=0.03,tempcriteria=0.3,pressref=10.):
    """
    Computes mixed layer depth, based on de Boyer Montagut et al., 2004.
    Parameters
    Downloaded from: https://github.com/pyoceans/python-oceans/blob/master/oceans/sw_extras/sw_extras.py
    and adapted to have variable threshold criteria and depth reference
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [deg (ITS-90)]
    p : array_like
        pressure [db].
    Notes
    -----
    Based on density with fixed threshold criteria
    de Boyer Montagut et al., 2004. Mixed layer depth over the global ocean:
        An examination of profile data and a profile-based climatology.
        doi:10.1029/2004JC002378
    dataset for test and more explanation can be found at:
    http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php
    Codes based on : http://mixedlayer.ucsd.edu/
    as it is if it doesn't find the criteria it will default to the depth of deepest density

    """
    import gsw
    import numpy as np
    import numpy.ma as ma

    #m = len(s)
    m = len(np.nonzero(~np.isnan(s))[0])
    if m <= 1:
        mldepthdens_mldindex = 0
        mldepthptemp_mldindex = 0
        return mldepthdens_mldindex, mldepthptemp_mldindex
    else:
        # starti = min(find((pres-10).^2==min((pres-10).^2)));
        starti = np.min(np.where(((p - pressref)**2 == np.min((p - pressref)**2)))[0])
        pres = p[starti:m]
        sal = s[starti:m]
        temp = t[starti:m]
        starti = 0
        m = len(sal)
        pden = gsw.rho(sal, temp,0)-1000

        mldepthdens_mldindex = np.argmax ( pres[ma.where(pden)] )
        for i, ppp in enumerate(pden):
            if np.abs(pden[starti] - ppp) > denscriteria:
                mldepthdens_mldindex = i
                break

        # Interpolate to exactly match the potential density threshold.
        presseg = [pres[mldepthdens_mldindex-1], pres[mldepthdens_mldindex]]
        pdenseg = [pden[starti] - pden[mldepthdens_mldindex-1], pden[starti] -
                pden[mldepthdens_mldindex]]
        P = np.polyfit(presseg, pdenseg, 1)
        presinterp = np.linspace(presseg[0], presseg[1], 3)
        pdenthreshold = np.polyval(P, presinterp)

        # The potential density threshold MLD value:
        ix = np.max(np.where(np.abs(pdenthreshold) <= denscriteria)[0])
        mldepthdens_mldindex = presinterp[ix]

        # Search for the first level that exceeds the temperature threshold.
        mldepthptmp_mldindex = np.argmax ( pres[ma.where(pden)] )
        for i, tt in enumerate(temp):
            if np.abs(temp[starti] - tt) > tempcriteria:
                mldepthptmp_mldindex = i
                break

    # Interpolate to exactly match the temperature threshold.
        presseg = [pres[mldepthptmp_mldindex-1], pres[mldepthptmp_mldindex]]
        tempseg = [temp[starti] - temp[mldepthptmp_mldindex-1],
                temp[starti] - temp[mldepthptmp_mldindex]]
        P = np.polyfit(presseg, tempseg, 1)
        presinterp = np.linspace(presseg[0], presseg[1], 3)
        tempthreshold = np.polyval(P, presinterp)

        # The temperature threshold MLD value:
        ix = np.max(np.where(np.abs(np.ndarray.round(tempthreshold,3)) <= tempcriteria)[0])
        mldepthptemp_mldindex = presinterp[ix]

        return mldepthdens_mldindex, mldepthptemp_mldindex







def zmld_boyervar_densonly(s, t, p,denscriteria=0.03,pressref=10.):
    """
    Computes mixed layer depth, based on de Boyer Montagut et al., 2004.
    Parameters
    Downloaded from: https://github.com/pyoceans/python-oceans/blob/master/oceans/sw_extras/sw_extras.py
    and adapted to have variable threshold criteria and depth reference
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [deg (ITS-90)]
    p : array_like
        pressure [db].
    Notes
    -----
    Based on density with fixed threshold criteria
    de Boyer Montagut et al., 2004. Mixed layer depth over the global ocean:
        An examination of profile data and a profile-based climatology.
        doi:10.1029/2004JC002378
    dataset for test and more explanation can be found at:
    http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php
    Codes based on : http://mixedlayer.ucsd.edu/
    """
    import gsw
    import numpy as np
    import numpy.ma as ma
    m = len(np.nonzero(~np.isnan(s))[0])
    if m <= 1:
        mldepthdens_mldindex = 0
        mldepthptemp_mldindex = 0
        return mldepthdens_mldindex, mldepthptemp_mldindex
    else:
        # starti = min(find((pres-10).^2==min((pres-10).^2)));
        starti = np.min(np.where(((p - pressref)**2 == np.min((p - pressref)**2)))[0])
        pres = p[starti:m]
        sal = s[starti:m]
        temp = t[starti:m]
        starti = 0
        m = len(sal)
        pden = gsw.rho(sal, temp,0)-1000

        mldepthdens_mldindex = np.argmax ( pres[ma.where(pden)] )
        for i, ppp in enumerate(pden):
            if np.abs(pden[starti] - ppp) > denscriteria:
                mldepthdens_mldindex = i
                break

        # Interpolate to exactly match the potential density threshold.
        presseg = [pres[mldepthdens_mldindex-1], pres[mldepthdens_mldindex]]
        pdenseg = [pden[starti] - pden[mldepthdens_mldindex-1], pden[starti] -
                pden[mldepthdens_mldindex]]
        P = np.polyfit(presseg, pdenseg, 1)
        presinterp = np.linspace(presseg[0], presseg[1], 3)
        pdenthreshold = np.polyval(P, presinterp)

        # The potential density threshold MLD value:
        ix = np.max(np.where(np.abs(pdenthreshold) < denscriteria)[0])
        mldepthdens_mldindex = presinterp[ix]
    
        return mldepthdens_mldindex

