from cProfile import Profile
from pstats import Stats
from io import StringIO
import pandas as pd

def print_stats(pr):
	s = StringIO()
	sortby = 'cumulative'
	ps = Stats(pr, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue()) 


pr = Profile()
pr.enable()

# Your code is here
def open_xls():
    for i in range(10):
        d =pd.ExcelFile('~/Documents/scripts/PTM/first/dataset.xls')
        
def open_csv():
    for i in range(10):
        d =pd.read_table('~/Documents/scripts/PTM/first/dataset.csv')        
open_csv()
open_xls()


pr.disable()    
print_stats(pr)
