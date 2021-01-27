from scipy.stats import chi2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Calculate and return lambda (inflation)
def __get_inflation__(observed):
    obsv_median = np.median(observed)
    Chi = chi2.ppf(1.0 - obsv_median, 1)
    lmbd = Chi / chi2.ppf(0.5, 1)
    return lmbd

# Returns: (fig, ax, lambda)
# fig and ax for more custermizations
def qqplot(filename,
           output='output.png',
           column_title = 'P',
           title='Q-Q plot',
           xlabel='Expected –log10 P-values',
           ylabel='Observed –log10 P-values',
           dpi=300):
   
    fh = open(filename, 'r')
    line = fh.readline()
    tmp = line.strip().split(',')

    # Read data into datafram
    df = pd.DataFrame()
    # Check delimiter type, "," or " "
    if len(tmp) == 1:
        df = pd.read_csv(filename, delim_whitespace=True)
    else:
        df = pd.read_csv(filename)

    Pobsv = df.loc[:, column_title]  # Observed p values
    Pobsv = np.sort(Pobsv.dropna().values)
    Pexpt = np.arange(1, len(Pobsv) + 1, 1) / (len(Pobsv) + 1)  # Expected
    logPobsv = -np.log10(Pobsv)
    logPexpt = -np.log10(Pexpt)

    # Calculate lambda
    infl = __get_inflation__(Pobsv)

    fig, ax = plt.subplots(dpi=dpi)
    ax.plot(logPexpt, logPobsv, linestyle='', marker='o', markersize=3, markeredgewidth=0.5,
            fillstyle='none', color='k', alpha=0.8)
    ax.plot(logPexpt, logPexpt, color='r', linewidth=0.4)

    # Label x and y axis
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    annotation = "λ = " + str("{0:.4f}".format(infl))
    ax.annotate(annotation, xy=(0.7, 0.2), xycoords='axes fraction')
    fig.savefig(output)

    return (fig, ax, infl)  # Return for more custermizations