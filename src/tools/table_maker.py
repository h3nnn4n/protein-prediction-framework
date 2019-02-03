import os
import re


names = os.listdir()

targets = ['all_best_fxn.dat', 'all_rmsd_fxn.dat']

data = {}
keys = []

ncols = 0

for name in names:
    if os.path.isdir(name) and len(name) == 4:
        os.chdir(name)
        for target in targets:
            key = name + '_' + target.split('.')[0]

            if key not in data.keys():
                data[key] = []
                keys.append(key)

            with open(target) as f:
                lines = f.readlines()

                for line in lines:
                    tokens = re.sub("\s+", " ", line.strip()).split(' ')

                    if ncols == 0:
                        ncols = len(tokens)

                    data[key].append(tokens)
        os.chdir('..')

w = ncols - 2

header = """\\begin{table}
  \centering
    \\begin{tabular}{l|""" + 'l|' * w + """l}
      \hline"""

tail = """  \end{tabular}
  \caption{%s}
  \label{result}
\end{table}"""


for key in keys:
    print(header)
    for d in data[key]:
        for k, v in enumerate(d):
            if k < len(d) - 1:
                print(v, end=' & ')
            else:
                print(v, end=' \\\\ \hline\n')
    caption = 'rmsd' if 'rmsd' in key else 'best'
    print(tail % caption)
    print()
