from easyIOmassspec import easyio
import matrix
import os

def main():

    parameters = {'absolute_intensity_thresh': 100.0,
                  'mz_factor': 10000}

    input_dir = '/Users/jasonzhou/Desktop/adap-3d/project/Data/'

    for file_name in os.listdir(input_dir):
        if file_name == '.DS_Store':
            continue
        df_str = input_dir + file_name
        dfr = easyio.DataFileReader(df_str, False)
        matrix_variables = matrix.Matrix(dfr, parameters)



        stop = 3


if __name__ == "__main__":
    main()