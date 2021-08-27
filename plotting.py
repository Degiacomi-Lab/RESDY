import pandas as pd
from os import name
import numpy
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import inquirer
import os


def mark_af():

    csv_name = input('Name of csv:')
    path = 'Output/' + csv_name
    all_pka_sasa_res = pd.read_csv(path)
    all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
    df_af = all_pka_sasa_res((all_pka_sasa_res['PDB Code'])[:3] == '/AF')


    return

def plot_mean_with_sd ():

    csv_name = input('Name of csv:')
    path = 'Output/' + csv_name
    all_pka_sasa_res = pd.read_csv(path)
    all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
    plt.errorbar(x, y, e, linestyle='None', marker='^')

    return

