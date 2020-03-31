#!/usr/bin/env python
######################################################################
# 
# DESCRIPTION
# ===========
#
# Library for reading ECMWF (and COSMO) data in grib1-format and
# extracting (and plotting) data fields
#
# ####################################################################


# load libraries -----------------------------------------------------
import pygrib
import numpy as np
from scipy import interpolate
# ====================================================================




def grib_cont_list(infile): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Creates a dictionary of the grib file content


    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file


    Returns
    --------
    L : dict
        dictionary of grib content
         * short names are used as keys
         * list of long names and number of occurence are given as list
           for each key
    '''

    print(infile)
    print('===========================================================')
# open input file ....................................................
    gf=pygrib.open(infile)
    gf.seek(0)


# create a list of variables .........................................
    L={};
    for dat in gf:
        na=dat.name
        sna=dat.shortName
        
# take short name as identifier
# append long name and number of levels 

        if sna in L:
            L[sna][1]=L.get(sna,0)[1]+1
        else:
            L[sna]=[na,0]
        
        
    return L # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################



def get_grib_field_lll(infile,var_name): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Extracts longitude, latitude, level and one 2d or 3d-field from grib file 


    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file

    var_name : str
        short name of variable in infile


    Returns
    --------
    lon : numpy array
        longitude

    lat : numpy array
        latitude
    
    lev : list
        level list

    var : numpy array
        field of variable var_name
    '''


# open and read from COSMOS file -------------------------------------
    Indx = pygrib.index(infile,'shortName')


# get certain variable ...............................................
    try:
        Cont=Indx.select(shortName=var_name)
    except:
        print(('ERROR: variable '+var_name + ' not in ' + infile))
        return


    Nx=Cont[0].Ni
    Ny=Cont[0].Nj
    Nz=np.shape(Cont)[0]

    lat,lon=Cont[0].latlons()
    
    lat=np.transpose(lat)
    lon=np.transpose(lon)

    if Nz > 1:
        var=np.zeros([Nx,Ny,Nz])
        lev=np.zeros([Nz])

        n=0
        for layer in Cont:
            var[:,:,n] = np.transpose(layer.values)
            lev[n] = layer.level
            n = n + 1

        # sort fields in vertical
        lev =  np.array(lev)
        ils = lev.argsort()  

        l = lev[ils]
        v = var[..., ils]

    else:
        layer = Cont[0]
        v = np.transpose(layer.values)
        l = [layer.level,]



    return lon,lat,l,v # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################


def get_grib_field(infile,var_name): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Extracts one 2d or 3d-field from grib file 


    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file

    var_name : str
        short name of variable in infile


    Returns
    --------
    var : numpy array
        field of variable var_name
    '''


# open and read from COSMOS file -------------------------------------
    Indx = pygrib.index(infile,'shortName')


# get certain variable ...............................................
    Cont=Indx.select(shortName=var_name)


    Nx=Cont[0].Ni
    Ny=Cont[0].Nj
    Nz=np.shape(Cont)[0]


    if Nz > 1:
        var = np.zeros([Nx,Ny,Nz])

        lev=np.zeros([Nz])
        n=0
        for layer in Cont:
            var[:,:,n] = np.transpose(layer.values)
            lev[n] = layer.level
            n = n + 1

        # sort fields in vertical
        lev =  np.array(lev)
        ils = lev.argsort()  

        l = lev[ils]
        v = var[..., ils]


    else:

        layer = Cont[0]
        v = np.transpose(layer.values)
        l = [layer.level,]


    return v # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################


def get_field_from_indx(indx, var_name, lonlat = False):

    '''
    Reads field from grib file given the grib index.


    Parameters
    ----------
    indx : pygrib index object 
        index for assessing variables in grib file (for 'shortName')

    var_name : str
        name of variable (shortName)

    lonlat : bool, optional, default = False
        option if longitude and latitude are also retrieved (default: False)


    Returns
    --------
    var : numpy array
        field of variable
    '''


    # get full data object with variable name ........................
    Cdat = indx.select(shortName = var_name)
    
    v = []
    lev = []
    
    # extract the layers of the field ................................
    for c in Cdat:
        v.append(c.values.T)
        lev.append(c.level)

    var = np.dstack(v)
    lev = np.array(lev)

    # sort fields in vertical
    ils = lev.argsort()  
    
    l = lev[ils]
    
    if np.ndim(var) == 3:
        v = var[..., ils]


    if lonlat:
        # obtain longitude and latitude ..............................
        lat,lon = c.latlons()

        lat = np.transpose(lat)
        lon = np.transpose(lon)
  
        return lon, lat, v.squeeze()
    else:
        return v.squeeze() 

######################################################################
######################################################################

        
def get_fc_fields(cfile, var_list, lonlat = True):

    '''
    Reads a list of variables from grib file.

    Parameters
    ----------
    cfile : str
        grib filename (full path)
    
    var_list : list
        list of variable names (shortName)

    lonlat : bool, optional, default = False
        option if longitude and latitude are also retrieved (default: False)


    Returns
    --------
    var : dict of numpy arrays
        dictionary of variable fields

    '''


    indx = pygrib.index(cfile,'shortName')

    var = {}
     
    for vname in var_list:
        print('... reading %s' % vname )

        try:
            if lonlat and 'lon' not in var:
                var['lon'], var['lat'], var[vname] = get_field_from_indx(indx, vname, lonlat=True)
            else:
                var[vname] = get_field_from_indx(indx, vname, lonlat=False)

        except:
            print('%s not in %s' %(vname, cfile))

    return var 

######################################################################
######################################################################

def get_full_grib_field_info(infile,var_name): # LLLLLLLLLLLLLLLLLLLLL

    '''
    Prints all keys and their corresponding values for a grib1 field.


    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file

    var_name : str
        short name of variable in infile


    Returns
    --------
    grb_info : dict
        dictionary of field keys and their values
    '''


# open and read from COSMOS file -------------------------------------
    Indx = pygrib.index(infile,'shortName')

    
# get certain variable ...............................................
    cont = Indx.select(shortName=var_name)

    # set first field 
    f = cont[0]

    # start with info dict
    grb_info={}

    for key in list(f.keys()):
        grb_info[key] = f[key]
        print((key,'  ',f[key]))
    
    return grb_info # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################



def get_grib_cont(infile): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Extracts of content of a whole grib file.

    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file


    Returns
    --------
    var : dict
        dictionary of field data
          * short names are used as keys
          * the fields are attributed to each key
    '''

    L = grib_cont_list(infile)

    cont={}
    for var in L:
        cont[var] = get_grib_field(infile,str(var))

    return cont # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################


def mlev2pres(infile,var): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


    '''
    Calculates the pressure of a ECMWF hybrid model level given as
    
                      p = a + b * p_surf


    Parameters
    ----------
    infile : str
        input filename of the concerned grib1 file

    var : str
        short name of variable in infile


    Returns
    --------
    pres : numpy array
        pressure levels
    '''


# open and read from COSMOS file -------------------------------------
    Indx = pygrib.index(infile,'shortName')



# get surface pressure ...............................................
    try:
        sp_cont = Indx.select(shortName='sp')
        sp = np.transpose(sp_cont[0].values)
    except:
        pass # no sp variable

    try:
        sp_cont = Indx.select(shortName='lnsp')
        lnsp = np.transpose(sp_cont[0].values)
        sp = np.exp(lnsp)

    except:
        pass # no lnsp variable
 

    try:
        sp_shape = sp.shape
    except:
        raise



# get level information for variable .................................
    var_cont = Indx.select(shortName=var)
    
    Nx = var_cont[0].Ni
    Ny = var_cont[0].Nj
    Nz = np.shape(var_cont)[0]
    lev = {}

    n = 0
    for layer in var_cont:
        lev[n] = layer.level
        n = n + 1

# get level parameters ...............................................
    ab = var_cont[0].pv
    Nlev_half = ab.shape[0] / 2


    # these are the corefficients of the half levels
    a, b = ab[0:Nlev_half], ab[Nlev_half:]

    # little check
    if a[0] != 0. or a[-1] != 0. or b[0] != 0. or b[-1] != 1.:
        print('Warning: Coefficient Check failed!')
        print(('a = ', a))
        print(('b = ', b))
        
    pres = np.zeros([Nx,Ny,Nz])

    if Nz != Nlev_half:
        # assume the variable is given on full levels
        a = 0.5 * (a[1:] + a[:-1])
        b = 0.5 * (b[1:] + b[:-1])


# pressure on half levels ............................................
    for n in range(Nz):
        
        pres[:,:,n] = a[n] + b[n] * sp[:,:]
    
    return pres


######################################################################
######################################################################

def calc_COSMO_turning(lon, lat, lonN = -170., latN = 40.):

    '''
    Calculate the angle which 2d vectors in COSMO have to be turned.

    This is needed as COSMO has a rotated coordinate system and saves
    and e.g. horizontal wind vector with respect to this and not to 
    the geographic reference frame. 


    Parameters
    ----------
    lon : numpy array
        longitude of grid point(s)

    lat : numpy array
        latitude of grid point(s)


    Returns
    --------
    delta : numpy array
        angle used to rotate 2d vectors
    '''

    # deg 2 rad tranformation
    lon  = lon / 180. * np.pi
    lat  = lat / 180. * np.pi

    lonN  = lonN / 180. * np.pi
    latN  = latN / 180. * np.pi
    

    arg =  np.cos(latN) * np.sin(lonN - lon) / \
        (np.cos(lat) * np.sin(latN) - np.sin(lat) * np.cos(latN) * np.cos(lonN-lon))
    delta = np.arctan( arg )
    
    return delta

######################################################################
######################################################################


def turn_COSMO_winds(delta, u, v):

    '''
    Rotates the components of the wind vectors given in the
    reference frame of COSMO to the geographical reference frame.


    Parameters
    ----------
    delta : numpy array
        rotation angle

    u : numpy arary
        wind in COSMOs x-direction 

    v : numpy array
        wind in COSMOs y-direction


    Returns
    --------
    ug : numpy arary
       zonal wind

    vg : numpy arary
       meridional wind
    '''

    # initialize wind fields
    ug = np.zeros(u.shape)
    vg = np.zeros(v.shape)

    # set sine and cosine of delta
    cos_delta = np.cos(delta)
    sin_delta = np.sin(delta)
    
    # start rotation in all levels
    for l in range(u.shape[2]):
        ug[:,:,l]  =   u[:,:,l] * cos_delta  +  v[:,:,l] * sin_delta
        vg[:,:,l]  = - u[:,:,l] * sin_delta  +  v[:,:,l] * cos_delta

    return ug, vg

######################################################################
######################################################################


def interpolate_to_plev(plev, p, var):

    '''
    Linearly interpolates a variable given on some vertical level type to
    a set of pressure levels.


    Parameters
    ----------
    plev : numpy array or list
        sorted vector of pressure levels

    p : numpy array
        pressure field in the coordinate system of the variable

    var : numpy array
        variable which should be interpolated


    Returns
    --------
    varlev : numpy array
        variable interpolated linearly to given pressure levels


    Notes
    ------
    This routine can also be used for the interpolation to arbitrary
    vertical level types, e.g. isentropic levels, when p and plev are
    replaced accordingly.
    '''

    # field dimensions
    Nx, Ny = var.shape[0:2]
    Np = plev.shape[0]

    # initialize interpolated field
    varlev = np.zeros((Nx,Ny,Np))
    

    # loop over vertical columns
    for i in range(Nx):
        for j in range(Ny):
            pz = p[i,j,:]
            varz = var[i,j,:]

            
            f = interpolate.interp1d(pz,varz,kind='linear',bounds_error=False)
            
            varlev[i,j,:] = f(plev)
            
    
    return varlev

######################################################################
######################################################################
    

