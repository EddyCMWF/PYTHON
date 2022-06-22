#!/usr/bin/python
#
# Python module containing UK site info
#

def get_site_dict_all():
    #
    site_dict_all = { 'AliceHolt':  { 'fname':'aliceholt', 'years':(2000,2005),      \
                                      'Jindex':0,'PFT_index':0,'Site_int':1800.,     \
                                      'CF_n':'FC','Rd_n':'Rn','H_n':'H','LE_n':'LE', \
                                      'CO2_range':(-10,30),'CO2_close_range':(-5,10), \
                                      'CO2_smF':4 },      \
                      'Auchencorth':{ 'fname':'amoss',  'years':(2002,2013),          \
                                      'Jindex':1,'PFT_index':2,'Site_int':1800.,      \
                                      'CF_n':'NEE','Rd_n':None,'H_n':'H','LE_n':'LE', \
                                      'CO2_range':(-10,20),'CO2_close_range':(-5,10), \
                                      'CO2_smF':6   },      \
                      'Cardington': {'fname':'cardington', 'years':(2005,2011),      \
                                     'Jindex':2,'PFT_index':2, 'Site_int':1800.,     \
                                     'CF_n':'NEE','Rd_n':'Rg','H_n':'H','LE_n':'LE', \
                                     'CO2_range':(-10,20),'CO2_close_range':(-5,5), \
                                     'CO2_smF':4   },      \
                      'EastSaltoun': {'fname':'eastsaltoun','years':(2003,2006),      \
                                      'Jindex':3,'PFT_index':2,'Site_int':3600.,      \
                                      'CF_n':'NEE','Rd_n':'Rn','H_n':'H','LE_n':'LE', \
                                      'CO2_range':(-10,20),'CO2_close_range':(-5,10), \
                                      'CO2_smF':8   },     \
                      'EasterBush': {'fname':'easterbush','years':(2003,2013),       \
                                     'Jindex':4,'PFT_index':2,'Site_int':1800.,      \
                                     'CF_n':'NEE','Rd_n':None,'H_n':'H','LE_n':'LE', \
                                     'CO2_range':(-10,20),'CO2_close_range':(-5,10), \
                                     'CO2_smF':4   },       \
                      'GriffinForrest': {'fname':'griffin','years':(1997,2010),          \
                                         'Jindex':5,'PFT_index':1,'Site_int':3600.,      \
                                         'CF_n':'NEE','Rd_n':'Rg','H_n':'H','LE_n':'LE', \
                                         'CO2_range':(-10,20),'CO2_close_range':(-5,10), \
                                         'CO2_smF':6    }   \
                  }
    #
    return site_dict_all
#

def get_site_dict(site_name):
    #
    site_dict_all=get_site_dict_all()
    site_dict=site_dict_all[site_name]
    #
    return site_dict

