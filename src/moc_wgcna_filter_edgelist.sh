# 2019-11-22
# Samuel Lee

# Filter wgcna edge list file 
# for minimum value of 0.1

cd "/c/Users/slee/projects/mocNets"


awk '{ if ($3 > 0.01) { print } }' "results/moc-wgcna-10th-power-MOC-adjacency-3column.txt" > \
"results/moc-wgcna-10th-power-MOC-adjacency-3column-filterDot01.txt"

awk '{ if ($3 > 0.01) { print } }' "results/moc-wgcna-10th-power-BEN-adjacency-3column.txt" > \
"results/moc-wgcna-10th-power-BEN-adjacency-3column-filterDot01.txt"
