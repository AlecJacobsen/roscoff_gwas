## This is not finished!!!!
# for testing purposes only

import argparse
import pandas as pd
from lunarcalendar import Converter
import statistics

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
samples.drop("indiv", axis = 1, inplace = True)
samples['b-strain'] = samples['strain'].apply(lambda x: 1 if x == 'FM' else 0)
samples['date'] = pd.to_datetime((samples['date'] + '-2019'),dayfirst = True)
samples['lunar_day'] = samples['date'].apply(lambda x: lunar_day(x,5))
samples['time'] = samples['time'].apply(lambda x: convert_time(x,12))
samples['status'] = [1]*len(samples['time'])
samples['dia_break'] = samples['date'].apply(lambda x: 1 if x > pd.to_datetime('2019-04-10') else 0)
peak_day = []
for i in range(len(samples['date'])):
    month = samples['dia_break'][i]
    strain = samples['b-strain'][i]
    peak_day.append((samples['date'][i] - min(samples[(samples['b-strain'] == strain) & (samples['dia_break'] == month)]['date'])).total_seconds() / (3600*24))
samples['peak_day'] = peak_day
off_peakness = []
for i in range(len(samples['date'])):
    month = samples['dia_break'][i]
    strain = samples['b-strain'][i]
    med = statistics.median(samples[(samples['b-strain'] == strain) & (samples['dia_break'] == month)]['peak_day'])
    off_peakness.append(abs(med - samples['peak_day'][i]))
samples['off_peakness'] = off_peakness
samples.to_csv(out_file + "_all.txt",sep = '\t', index = False)

#strain file
strain = samples[['FID','IID','strain','b-strain']]
strain.to_csv(out_file + '_strain.txt',sep = '\t', index = False)
