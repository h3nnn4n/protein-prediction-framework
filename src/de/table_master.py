import os
import re


names = os.listdir()

data = {}

proteins = []

ncols = 0

for name in names:
    if os.path.isdir(name) and len(name) == 4:
        os.chdir(name)
        proteins.append(name)
        target = 'all_data.dat'
        with open(target) as f:
            lines = f.readlines()

            for line in lines:
                tokens = re.sub("\s+", " ", line.strip()).split(' ')

                test = tokens[0]

                if test not in data.keys():
                    data[test] = {}

                mode = tokens[1]

                data[test][mode] = {}

                data[test][mode]['best'] = float(tokens[2])
                data[test][mode]['mean'] = float(tokens[3])
                data[test][mode]['stddev'] = float(tokens[4])
                data[test][mode]['median'] = float(tokens[5])
        os.chdir('..')

w = 3

header = """\\begin{table}
  \centering
    \\begin{tabular}{l|""" + 'l|' * w + """l}
      \hline"""

tail = """  \end{tabular}
  \caption{%s}
  \label{result}
\end{table}"""

p1zdd = """ 1ZDD    & DE$_{C_1}                 $& 54.27               & 7.67              & $82.97\pm15.49     $  \\\ \hline
        & DE$_{C_2}                 $& 65.77               & 9.42              & $82.76\pm9.22      $  \\\ \hline
        & GA  			     & -40.40              & 10.9              & $-36.20\pm2.60     $  \\\ \hline
        & MSA			     & -62.99        	   & 2.62              & $-48.96\pm7.77     $  \\\ \hline """

p1crn = """ 1CRN    & DE$_{C_1}                 $& 82.86               & 21.56             & $126.95\pm25.98    $  \\\ \hline
        & DE$_{C_2}                 $& 72.48               & 15.44             & $109.08\pm22.96    $  \\\ \hline
        & GA   		             & -22.70              & 5.8               & $-18.20\pm2.9      $  \\\ \hline
        & MSA		             & -76.93        	   & 6.96              & $-54.01\pm 17.30   $  \\\ \hline """

p1enh = """ 1ENH    & DE$_{C_1}                 $& 294.25              & 14.72             & $372.11\pm52.05    $  \\\ \hline
        & DE$_{C_2}                 $& 255.54              & 19.28             & $320.38\pm41.06    $  \\\ \hline
        & GA                         & -56.08              & 14.99             & $-51.52\pm1.94     $  \\\ \hline
        & MSA                        & -95.86              & 5.70              & $-80.75\pm8.48     $  \\\ \hline """

p1ail = """ 1AIL    & DE$_{C_1}                 $& 357.84              & 25.00             & $440.63\pm58.11    $  \\\ \hline
        & DE$_{C_2}                 $& 332.54              & 16.88             & $411.81\pm56.84    $  \\\ \hline
        & GA                         & -75.07              & 12.34             & $-71.08\pm3.35     $  \\\ \hline
        & MSA			     & -128.55             & 8.27              & $-117.54\pm10.28   $  \\\ \hline """

transl = {'1zdd': p1zdd, '1crn': p1crn, '1enh': p1enh, '1ail': p1ail}


for protein in proteins:
    br = float('inf')
    bs = float('inf')
    bm = float('inf')

    brn = None
    bsn = None
    bmn = None

    print(header)
    print(transl[protein])

    for k, v in data.items():
        if protein in k:
            if v['best_fxn']['best'] < bs:
                bs = v['best_fxn']['best']
                bsn = k

            if v['best_fxn']['mean'] < bm:
                bm = v['best_fxn']['mean']
                bmn = k

            if v['rmsd_fxn']['best'] < br:
                br = v['rmsd_fxn']['best']
                brn = k

    for k, v in data.items():
        if protein in k:
            ss = [("& %s &", k[5:].replace('_', '-')),
                  ("_%.2f#&", v['best_fxn']['best']),
                  ("_%.2f#&", v['rmsd_fxn']['best']),
                  ("$_%.2f \pm %.2f#$ \\\ \hline", (v['best_fxn']['mean'], v['best_fxn']['stddev']))]

            print(ss[0][0] % ss[0][1], end='')

            if k == bsn:
                w = ss[1][0] % ss[1][1]
                w = re.sub('_', '\\\\textbf{', w)
                w = re.sub('#', '}', w)
                print(w, end='')
            else:
                w = ss[1][0] % ss[1][1]
                w = re.sub('_', ' ', w)
                w = re.sub('#', ' ', w)
                print(w, end='')

            if k == brn:
                w = ss[2][0] % ss[2][1]
                w = re.sub('_', '\\\\textbf{', w)
                w = re.sub('#', '}', w)
                print(w, end='')
            else:
                w = ss[2][0] % ss[2][1]
                w = re.sub('_', ' ', w)
                w = re.sub('#', ' ', w)
                print(w, end='')

            if k == bmn:
                w = ss[3][0] % ss[3][1]
                w = re.sub('_', '\\\\bm{', w)
                w = re.sub('#', '}', w)
                print(w, end='')
            else:
                w = ss[3][0] % ss[3][1]
                w = re.sub('_', ' ', w)
                w = re.sub('#', ' ', w)
                print(w, end='')

            print()


            # print(" & %s & %.2f & %.2f & $%.2f \pm %.2f$ \\\ \hline" %
                  # (k[5:].replace('_', '-'), v['best_fxn']['best'], v['rmsd_fxn']['best'], v['best_fxn']['mean'],  v['best_fxn']['stddev']))

    print(tail % protein.upper())
    print()
