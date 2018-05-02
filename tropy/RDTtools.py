#!/usr/bin/python

# Libraries
# ==========
import sys, os
import numpy as np
import h5py
import datetime
from scipy.interpolate import griddata
   
#####################################################################
#####################################################################

class Cell(object):

    def __init__(self):
        self.lon = np.array([])
        self.lat = np.array([])
        return
    # ................................................................


    def wind(self):
        al = np.pi * self.direction / 180.
        U = self.speed

        # zonal and meridional wind
        self.u = -U * np.sin(al)
        self.v = -U * np.cos(al)

        self.p = self.pressure
        return

    # ................................................................

    def __getitem__(self, attr):
        return getattr(self, attr) 

    # ................................................................

    def keys(self):
        return  self.__dict__.keys()

    # ................................................................

    def nkeys(self):
        
        #  key set
        ks = self.keys()

        # only keys for numpy fields
        nk = []

        tnp = type(np.array([]))

                
        for k in ks:
            var = self.__dict__[k]

            if type(var) == tnp:
                nk.append(k)

        return nk

    # ................................................................

    def cell_connected(self, next_cell):

        # get id of cell
        cell_id = self.__dict__['id[32]']

        # get previous id of next cell
        next_cell_prev_id = next_cell['prev_id[32]']

        if cell_id == next_cell_prev_id:
            return True
        else:
            return False
    # ................................................................

    def get_cs_prop(self, prop):

        # check for type of the hfile argument
        thf = type(self.hfile)


        if thf == type('str'):
            hfile = self.hfile
        elif thf == type([]):
            hfile = self.hfile[0]

        fprod = h5py.File(hfile,'r')
        
        cs = self.cloud_systems[0]

        p = fprod[cs][prop].value

        fprod.close()    

        return p

    # ................................................................

    def contour(self):
        return self.get_cs_prop('CONTOUR')

    # ................................................................


    
    def area(self):
        return self.get_cs_prop('AREA')

    # ................................................................


######################################################################
######################################################################


class CellTrajectory(Cell):
    
    def __init__(self, name, start_date, first_cell):

        self.name = name

        self.start_date = start_date

        self.hfile = [first_cell.hfile]

        ks =  first_cell.keys()
        ks.remove('hfile')

        for k in ks:
            self.__dict__[k] = first_cell[k]


        return

    # ................................................................

    def append(self, next_cell):
        
        nk = self.nkeys() # get key for numpy sets

        for k in nk:

            # get old variables
            old_vars = self.__dict__[k]
            
            # get the new variable
            new_var = next_cell[k]
            
            # append the new one
            self.__dict__[k] = np.append(old_vars, new_var)

        self['hfile'].append(next_cell['hfile'])
            
        return


    # ................................................................

    def item(self, n):
        
        dset = self
        
        # get key set
        ks = dset.__dict__.keys()

        # data type
        tnp = type(np.array([]))

        ditem = Cell()

        for k in ks:
            var = dset.__dict__[k]
            
            if type(var) == tnp:
                ditem.__dict__[k] = np.array([var[n]])

            elif type(var) == type([]):
                ditem.__dict__[k] = var[n]

            else:
                ditem.__dict__[k] = var
        
        return ditem
        
    # ................................................................

    def last(self):
                
        return self.item(-1)
    # ................................................................


    def contour(self,n):
        
        cell = self.item(n)

        return cell.contour()


    # ................................................................


    
    def area(self,n):
        
        cell = self.item(n)


        return cell.area()

    # ................................................................
    
    def times(self):
        t = []
        for hfile in self.hfile:
            fname = os.path.basename(hfile)
            # ATTENTION  !!! : filename nomenclature assumed !
            time_str = fname[17:29]
            try:
                time = datetime.datetime.strptime(time_str, '%Y%m%d%H%M')
                t.append(time)
            except:
                print 'Warning: filename nomenclature might have changed!'

        return t
        
    
######################################################################
######################################################################

class RDT(Cell):

    # ................................................................

    def __init__(self, hfile):
        self.hfile = hfile
        return

    # ................................................................

    def read(self):

        # check if file exists
        if not h5py.is_hdf5(self.hfile):
            print 'No hdf-file ',self.hfile ,' available!'
            raise IOError


        
        fprod = h5py.File(self.hfile,'r')

        # get header of rdt list
        rdt_head=fprod['RDT'][0].dtype.names

        # get values
        dset = fprod['RDT'].value

        # get length of full data set 
        Nlen = len(dset['lon'])

        # get names of cloud systems
        cs = []
        for n in range(Nlen):
            cs.append('RDT_CLOUD_SYSTEM_%06d' % n)

        self.__dict__['cloud_systems'] = np.array(cs)

        # get properties of cloud systems
        for name in rdt_head:
            self.__dict__[name] = np.array(dset[name])

        fprod.close()    

        return


    # ................................................................

    def mask(self, m):

        nk = self.nkeys()
        
        
        for k in nk:
            var = self.__dict__[k]
            self.__dict__[k] = var[m==True]

        return
 
    # ................................................................


    def item(self, n):
        
        dset = self
        
        # get key set
        ks = dset.__dict__.keys()

        # data type
        tnp = type(np.array([]))

        ditem = Cell()

        for k in ks:
            var = dset.__dict__[k]
            
            if type(var) == tnp:
                ditem.__dict__[k] = np.array([var[n]])
            else:
                ditem.__dict__[k] = var
        
        return ditem
        
    # ................................................................



######################################################################
######################################################################


class RDT_time_series(dict):

    # ................................................................

    def __init__(self, sat_type, arch):

        # set time interval
        if sat_type in ('msg2','meteosat9'):
            self.time_increment = datetime.timedelta(minutes = 15)
            self.sat_type = 'MSG2'
            self.scan_type = 'hrs'

        elif sat_type in ('msg1',  'meteosat8'):
            self.time_increment = datetime.timedelta(minutes = 5)
            self.sat_type = 'MSG1'
            self.scan_type = 'rss'

        else:
            self.sat_type = 'UNKNOWN'

        self.arch_dir = arch + '/'

        return

    # ................................................................

    def read(self, date1, date2, region = 'eu'):

        if self.sat_type == 'UNKNOWN':
            print 'unknown satellite spec'
            return

        # create region string .......................................
        tmpl_str = '________'
        region_str = region + tmpl_str[len(region):]
        

        # do loop over all time slots ................................
        date = date1
        
        while (date <= date2):

            date_string = date.strftime('%Y%m%d%H%M')            
            date_string_nice = date.strftime('%Y-%m-%d_%H:%M')            
 
            print 'reading RDT for %s' % date_string_nice
            # create filename
            fname = self.arch_dir + 'SAFNWC_%s_RDT__%s_%s-%s.buf.h5' \
                % (self.sat_type, date_string, self.scan_type, region_str)
            
            try:
                self[date_string_nice] = RDT(fname)
                self[date_string_nice].read()

            except IOError:
                del self[date_string_nice]
                print '... jump to next file'
 
            date += self.time_increment

        return

    # ................................................................

    def geo_mask(self, lon, lat):
        
        for date in self:

            # take data set
            h = self[date]


            # check if numpy arrays are there at all ...
            if not h.nkeys():
                continue

            # make mask
            i = 5
            mask = griddata((lon[::i,::i].flatten(),lat[::i,::i].flatten()),lat[::i,::i].flatten()*0,zip(h.lon,h.lat),method='linear',fill_value=1) == 0
            
            h.mask(mask)
        return


    # ----------------------------------------------------------------
    def mask_with_region(self, reg):

        for date in self:

            # take data set
            h = self[date]


            # check if numpy arrays are there at all ...
            if not h.nkeys():
                continue

            # make mask
            
            (lon1, lon2), (lat1, lat2) = reg

            lon_mask = np.logical_and(h.lon >= lon1, h.lon <=lon2)
            lat_mask = np.logical_and(h.lat >= lat1, h.lat <=lat2)

            mask = np.logical_and(lon_mask, lat_mask)

            
            h.mask(mask)

        return
 
    # ................................................................

    def connectivity(self):

        Rset = self

        # try to build a connection tree
        dates = Rset.keys()
        dates.sort()
        
        # 1st date
        date1 = dates[0]

        # all dates without the 1st
        ndates = dates[1:]

        # start connectivity dict
        connectivity = {}

        # loop over all dates wo the 1st
        for date in ndates:
            
            date0 = date1
            date1 = date
	
            # IDs in the last RDT data set
            id1 = list(Rset[date0]['id[32]'])

            # IDs in the current data set
            # which show the connection to the previous one
            id2 = list(Rset[date1]['prev_id[32]'])

            # index set
            iset = []

            for n, idval in enumerate(id2):

                if idval in id1:
                    i = id1.index(idval)
                else:
                    i = None
                    

                iset.append(i)


            connectivity[date] = iset

        return connectivity 

######################################################################
######################################################################


class RDT_traj(dict):
    
    def __init__(self, Rset):


        # get connectivity 
        con = Rset.connectivity()


        # get all the dates in the right order
        dates  = Rset.keys()
        dates.sort()

        date1 = dates[0]
        ndates = dates[1:]

        # number of all trajectoires
        Nall = 0

        # start with initial date ------------------------------------
        cur_traj = []
        r = Rset[date1]

        d = datetime.datetime.strptime(date1,'%Y-%m-%d_%H:%M')

        for n, lon in enumerate(r.lon):

            Nall += 1
            name = 'traj%04d' % Nall
            self[name] = CellTrajectory(name, d, r.item(n))

            cur_traj.append(name)
        # ============================================================


        # proceed ----------------------------------------------------
        for date in ndates:
            
            r = Rset[date]
            
            prev_traj = cur_traj
            cur_traj = []

            d = datetime.datetime.strptime(date1,'%Y-%m-%d_%H:%M')
            

            # get index set from the connectivity data
            iset = con[date]
            
            
            # check for connectivity for each cell
            for n, i in enumerate(iset):

                # append cell to trajectory if connection is given
                if i or i == 0:
                    name = prev_traj[i]

                    self[name].append(r.item(n))
                    
                # start a new trajectory if no connection is there
                else:
                    Nall+=1
                    name = 'traj%04d' % Nall
                    self[name] = CellTrajectory(name, d, r.item(n))

                cur_traj.append(name)
    # ================================================================


    def remove_short(self, min_len = 2):
            
        # get list of trajectory names ...........................
        traj_name_list = sorted(self.keys())


        for traj_name in traj_name_list:

            traj = self[traj_name]
                
            length = len(traj.lon)
            
            if length < min_len:
                del self[traj_name]
                

        return
    

######################################################################
######################################################################



def draw_RDT_contour(con, ecol, m):

     # transform contour coordinates to map projection
    x, y = m(con['longitude'],con['latitude'])
   
        
    # plug it to cuurent figure
    m.plot(x,y,color='black', ls = '-', lw=3)
    
    # plug it to cuurent figure
    m.plot(x,y,color=ecol, ls = '-')
    
    return

######################################################################
######################################################################


def select_color_for_cs(h, i, crit = 'btmin'):
    
    if crit == 'btmin':
        min_bt = h['min_bt'][i] - 273.15


        if min_bt>-20:
            ecol='GreenYellow'
        elif min_bt > -40:
            ecol='red'
        elif min_bt<-40:
            ecol='magenta'
            
    elif crit == 'nature':
        nat = h['nature'][i]

        if nat == 0:
            ecol = 'red'
#        elif nat == 2:
#            ecol = 'orange'
        else: 
            ecol = 'black'

    return ecol

######################################################################
######################################################################


def read_RDT_traj(date1, date2, rdt_dir, 
                  sat_type = 'msg2',
                  domain = ((6.,15.), (47., 55.)),
                  region = 'eu',
                  remove_single = True,
                  traj_min_length = 2):
    
    '''
    Reads a set of RDT trajectories.


    USAGE:
    =====
    trajectories = read_RDT_traj(date1, date2, rdt_dir, 
                                 sat_type = 'msg2',
                                 domain = ((6.,15.), (47., 55.)),
                                 region = 'eu',
                                 remove_single = True):


    INPUT:
    ======
    date1: datetime object of start time
    date2: datetime object of end time
    rdt_dir: name of full path to rdt files
    domain: tuble of min/max of lon and lat (OPTIONAL)
            in format ((lon_min, lon_max), (lat_min, lat_max))
    region: string of region used in NWCSAF  (OPTIONAL)
            (last entry in rdt filename)
    remove_single: flag if to remove trajectories with just 
                   one single entry  (OPTIONAL)

    


    OUTPUT:
    =======
    trajectories: dictionary of rdt cell trajectories
    '''


    # initialize RDT set .............................................
    Rset = RDT_time_series(sat_type, rdt_dir)

    # read data ......................................................
    Rset.read(date1,date2, region = region)
    
    # mask data ......................................................
    Rset.mask_with_region(domain)

    # sort data ......................................................
    trajectories = RDT_traj(Rset) 

    if remove_single:
        trajectories.remove_short(min_len = traj_min_length)

    return trajectories


######################################################################
######################################################################

    
    
                  

if __name__ == '__main__':

    IS_TESTED = True

    if not IS_TESTED:
        pass
    else:
        sys.path.append('/u1/home/fabian/lib/tropy/')
        import ECMWFtools
    # for testing
        cfile='/u1/home/fabian/data/COSMO/forecast/2012/06/30/06/lfff00110000'
    
        lon,lat,lev, qi3d = ECMWFtools.get_grib_field_lll(cfile,'QI')
        

        rdt_dir = '/u1/home/fabian/data/NWCSAF/2011/06/06/RDT_Tw-40_Tc-40_fake_06/'
        
        date1 = datetime.datetime(2011,6,6,9,0)
        date2 = datetime.datetime(2011,6,6,18,0)

        # rdt_dir = '/u1/home/fabian/data/NWCSAF/2012/07/05/RDT_Tw-40_Tc-40_fake_00/'

        # date1 = datetime.datetime(2012,7,5,9,0)
        # date2 = datetime.datetime(2012,7,5,18,0)


        traj = read_RDT_traj(date1, date2, rdt_dir, region = 'eu')

        
            
        
