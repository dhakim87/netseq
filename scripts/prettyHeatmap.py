import seaborn as sns
import csv
import matplotlib.pylab as plt
from collections import defaultdict;

counter = defaultdict(int)

header = None
colA = "score"
colB = "mstScore"

Xs = []
Ys = []

with open("group_105.cdt") as inputFile:
    reader = csv.reader(inputFile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for row in reader:
        if header is None:
            header = row
            colA = header.index(colA)
            colB = header.index(colB)
            continue
        
        counter[(int(float(row[colA])), int(row[colB]))] += 1
        Xs.append((int(float(row[colA]))))
        Ys.append(int(row[colB]))

heatmap = []
for y in range(50,130): #mst score
    row = []
    for x in range(120,190):
        row.append(counter[(x,y)])
    heatmap.append(row)



ax = sns.heatmap(heatmap, linewidth=0.5, xticklabels=range(120,190), yticklabels=range(50,130), square=True, annot=False)
plt.title("Observed")
plt.show()

plt.scatter(Xs, Ys, alpha=0.01)
plt.axvline(x=175, c='k', linestyle='dashed')
plt.axhline(y=80, c='k', linestyle='dashed')
plt.xlabel("Score")
plt.ylabel("MST Score")
plt.title("Candidate Network Score Distribution")
plt.figtext(0.5, 0.01, "Dotted Axes Identify \"True\" Network", wrap=True, horizontalalignment='center', fontsize=12)
plt.show();
#Score 175, Pair Score: 105, MST Score 80
