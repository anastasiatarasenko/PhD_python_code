pn = '/nfs1/store/satellite/sentinel-1/'
pol = 'hh'
resolution = 500    # in meters


# Section extraction
distance = 10e4
every_x = 50

# The large dictionar where S1 names are keys, and xs, ys - (xstart, ystart) describe 
# where from the profiles should start

fn_list_off={}


#### Baffin sea, off ice:
# # +
fn_list_off['S1A_EW_GRDM_1SDH_20170306T212543_20170306T212643_015577_0199B7_AD06.zip'] = {'area_x' :[200, 800],
                                                                                      'area_y' :[200, 600],
                                                                                      'xs' :[300, 400, 500],
                                                                                      'ys' :[100, 150, 200]
                                                                                       }




# +
#off-ice

fn_list_off['S1A_EW_GRDM_1SDH_20170311T065637_20170311T065737_015641_019B9F_8EBE.zip']={"distance":22e4,
                                                                                           "area_x": [0, -1],
                                                                                           "area_y": [500, -1],
                                                                                           "xs": [200, 400, 600],
                                                                                           "ys": [250, 250, 250]
                                                                                          }
#xs = [200, 400, 600]
#ys = [750, 750, 750]




#off-ice, 12 March 2017
fn_list_off['S1A_EW_GRDM_1SDH_20170312T073652_20170312T073756_015656_019C11_21FD.zip']={"area_x": [300, 800],
                                                                                    "area_y": [0, 400],
                                                                                    "xs": [50, 200, 350], 
                                                                                    "ys": [350, 350, 350]
                                                                                   }


 #off ice
fn_list_off['S1A_EW_GRDM_1SDH_20170330T050948_20170330T051048_015917_01A3D7_3D6F.zip'] ={"area_x" : [0, 400],
                                                                                         "area_y" : [0, 400],
                                                                                         "xs" : [200, 300, 300],
                                                                                         "ys" : [200, 250, 300],
                                                                                         "every_x" : 20}
# xs = [450, 450, 450]
# ys = [400, 500, 600]


fn_list_off['S1A_EW_GRDM_1SDH_20160818T061339_20160818T061439_012651_013DC8_FBB4.zip'] = {'area_x' :[600, -1],
                                                                                      'area_y' :[0, 300],
                                                                                      'ys' :[200, 225, 250],
                                                                                      'xs' :[25, 50, 75], 
                                                                                      'no filter':True
                                                                                       }


fn_list_off['S1A_EW_GRDM_1SDH_20160818T061339_20160818T061439_012651_013DC8_FBB4.zip']={'area_x' :[500, -1],
                                                                                    'area_y' :[0, 400],
                                                                                    "ys": [210, 230, 260],
                                                                                    "xs": [100, 150, 200]
                                                                                       }

