from matplotlib import pyplot as plt
import numpy as np
from glob import glob
from collections import Counter
import re

ex_folder = 'fline_data/'
experiments = ['no_1000', 'add_20_1000', 'add_50_1000', 'replace_20_1000', 'replace_50_1000']
for i in range(len(experiments)):
    experiments[i] = ex_folder + experiments[i]
    files = glob(experiments[i] + '/*')

    x = np.arange(1, 18)
    y = np.zeros([10, 17])
    for k, file in enumerate(files):
        with open(file) as f:
            line = f.readline()
            line = line.lstrip('[')
            line = line.rstrip(']')
            line = line.split(', ')
            line = np.sort(np.array(line, dtype=int))
            count = Counter(line)
            labels, values = zip(*count.items())
            for j, l in enumerate(labels):
                y[k][l-1] += values[j]
    y_mean = np.mean(y, axis=0)
    y_std = np.std(y, axis=0)

    plt.errorbar(x, y_mean, yerr=y_std, linestyle='--', marker='o', barsabove=True, capsize=4)
    plt.plot(x, y_mean, linestyle='--', marker='o')

    plt.yscale('log')
    plt.title('Average number of broken lines over 10 runs')
    plt.xlabel('Number of lines')
    plt.ylabel('Number of events')
    plt.legend(['no solar panels', '20% addition', '50% addition', '20% replacement', '50% replacement'])
plt.show()
