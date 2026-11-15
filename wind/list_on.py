# List of Sentinel-1 files and some key parameters to read 

#Default parameters:
pn = '/nfs1/store/satellite/sentinel-1/'
pol = 'hh'
resolution = 500    # in meters


# Section extraction
distance = 10e4
every_x = 50

# The large dictionar where S1 names are keys, and xs, ys - (xstart, ystart) describe 
# where from the profiles should start
fn_list_on={}

#Greenland, 1 march 
# on ice wind, need other cross-section analysis
fn_list_on['S1B_EW_GRDM_1SDH_20170301T072943_20170301T073043_004512_007DAB_A371.zip'] ={"area_x" :[0, 500],
                                                                                        "area_y" :[400, 700],
                                                                                        "xs" : [250, 300, 350],
                                                                                        "ys" : [150, 200, 250]
                                                                                        }

# +
#Greenland, 3 march
# on ice
fn_list_on['S1A_EW_GRDM_1SDH_20170303T174659_20170303T174803_015531_019850_434A.zip']={"xs": [330, 380, 430],
                                                                                       "ys": [350, 500, 650]
                                                                                       }

fn_list_on['S1A_EW_GRDM_1SDH_20170305T141328_20170305T141428_015558_019922_5072.zip']={'area_x' :[500, -1],
                                                                                    'area_y' :[300, 600],
                                                                                    'ys' :[100, 50, 0],
                                                                                    'xs' :[100, 125, 150]
                                                                                       }

### Central Arctic 

## Arctic Zn, on-ice conditions, beautiful ice displacement
fn_list_on['S1B_EW_GRDM_1SDH_20160920T194910_20160920T195015_002157_0039E2_21CD.zip']={"xs" :[400, 600, 800],
                                                                                        "ys" :[800, 600, 400]}

