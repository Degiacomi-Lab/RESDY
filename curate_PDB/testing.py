import csv
import pandas as pd
list1 = ['dog', 'cat', 'fox']

df = pd.DataFrame(list1)

df.to_csv('animals.csv', index=False)