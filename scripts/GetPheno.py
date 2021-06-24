## This is not finished!!!!
# for testing purposes only

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="Input File")
parser.add_argument("-o", "--output", dest="output", help="Output Prefix")
args = parser.parse_args()

samples_file = args.input
out_file = args.output

def convert_time(cell,shift):
    if type(cell) == str:
        cell = cell.strip('[RuSi]')
        time = int(cell.split('-')[0])+(float(cell.split('-')[1])/60)
        if time < 12:
            time += 24
        return(time - shift)
    pass

def lunar_day(day,shift):
    lunar_day = Converter.Solar2Lunar(pd.to_datetime(day)).day
    if lunar_day < shift:
        lunar_day += 30
    return(lunar_day)

samples = pd.read_csv(samples_file,header = None)
samples = pd.concat([pd.DataFrame([0]*len(samples)),samples,pd.DataFrame(samples[0].str.split('_').tolist())], axis = 1)
samples.columns = ['FID','IID','strain','date','time','indiv']
samples['date'] = pd.to_datetime((samples['date'] + '-2019'),dayfirst = True)
samples['time'] = samples['time'].apply(lambda x: convert_time(x,12))
samples.to_csv(out_file + '_all.txt',sep = '\t', index = False)


#strain file
strain = samples[['FID','IID','strain']]
strain.to_csv(out_file + '_strain.txt',sep = '\t', index = False)
#strain covariates
#samples[['FID','IID','time']].to_csv('strain_cov.txt',sep='\t',index = False)

#lunar day file
