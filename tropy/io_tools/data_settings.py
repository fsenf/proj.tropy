#!/usr/bin/env python
#-*- coding:utf-8 -*-


######################################################################
######################################################################

def bt_settings():

    sett = {}

    sett['add_offset'] = 0
    sett['scale_factor'] = 0.01
    sett['dtype'] = '|u2'
    sett['_FillValue'] = 0
    sett['longname'] = 'Brightness Temperature'
    sett['unit'] = 'K'

    return sett

######################################################################
######################################################################

def refl_settings():

    sett = {}

    sett['add_offset'] = - 0.2
    sett['scale_factor'] = 0.0001
    sett['dtype'] = '|u2'
    sett['_FillValue'] = 0
    sett['longname'] = 'Reflectivity'
    sett['unit'] = ''

    return sett

######################################################################
######################################################################


def histogram_settings():
    
    unit = 'counts'
    longname = 'histogram of absolute occurence frequency'

    return standard_int_setting(longname = longname, unit = unit)

 
######################################################################
######################################################################


def Dge_settings():
    
    unit = 'um'
    longname = 'generalized effective diameter'

    return standard_real_setting(longname = longname,  unit = unit)

 
######################################################################
######################################################################


def standard_real_setting(longname = 'Unknown', unit = 'Unknown'):
    
    sett = {}
    sett ['add_offset'] = 0
    sett['scale_factor'] = 1
    sett['dtype'] = '|f4'
    sett['_FillValue'] = - 99.
    sett['unit'] = unit
    sett['longname'] = longname

    return sett

 
######################################################################
######################################################################

def standard_int_setting(longname = 'Unknown', unit = 'Unknown'):
    
    sett = {}
    sett ['add_offset'] = 0
    sett['scale_factor'] = 1
    sett['dtype'] = '|u4'
    sett['_FillValue'] = -1
    sett['unit'] = unit
    sett['longname'] = longname

    return sett

 
######################################################################
######################################################################


def settings(vartype = 'bt'):

    if vartype == 'bt':
        return bt_settings()

    elif vartype == 'hist':
        return histogram_settings() 

    elif vartype == 'Dge':
        return Dge_settings() 
       

######################################################################
######################################################################

