import sys
sys.path.append("../")
import easyio as eio
import pylab as pl

def main():

#    mz_arr = pl.array([100.0,
#        101.0,
#        200.0,
#        201.0,
#        300.0,
#        301.0])
#    int_arr = pl.array([100,
#        100,
#        200,
#        200,
#        300,
#        300])
#    rt_arr = pl.array([1.0,10.0,20.0])
#    scan_index = pl.array([1,3,5]).astype(int)
#    act_scan_num =pl.array([1,3,5]).astype(int) 
#    tot_int_arr = pl.array([int_arr[0]+int_arr[1],int_arr[2]+int_arr[3],int_arr[4]])
#
#    points_in_scan_arr = pl.array([2,2,2]).astype(int)
    mz_arr = [100.0,
        101.0,
        200.0,
        201.0,
        300.0,
        301.0]
    int_arr = [100.0,
        100.0,
        200.0,
        200.0,
        300.0,
        300.0]
    rt_arr = [1.0,10.0,20.0]
    #scan_index = [0,2,4]
    act_scan_num =[1,3,5]
    tot_int_arr = [int_arr[0]+int_arr[1],int_arr[2]+int_arr[3],int_arr[4]]

    points_in_scan_arr = [2,2,2]


    num_scans = len(points_in_scan_arr )

    dfwrite = eio.DataFileWriter('test_out',num_scans)


    dfwrite.set_write_points_in_scan_arr(points_in_scan_arr )
    dfwrite.set_mz_values(mz_arr)
    dfwrite.set_intensity_values(int_arr)
    #dfwrite.set_scan_index(scan_index)
    dfwrite.set_scan_index()
    dfwrite.set_rt_values(rt_arr)
    dfwrite.set_tot_intensity_arr(tot_int_arr)
    #dfwrite.set_act_scan_num(act_scan_num)

    dfwrite.closeWriter()



if __name__ == "__main__":
    main()
