from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sciezka = "/home/markowar/bios/1prr.pdb"
p = PDBParser()
strukura=p.get_strukura("1MBO",sciezka)
print(strukura)

for item in strukura.get_models():
    print(item)
model=strukura[0]
dssp = DSSP(model,sciezka,dssp='/usr/bin/dssp')



kolory = {'H': 'red','B':'orange','E':'yellow', 'G':'blue','I':'violet','T': 'pink','S':'green', '-': 'black'}
for amino in list(dssp.keys()):
    plt.scatter(dssp[amino][4], dssp[amino][5], c=kolory[dssp[amino][2]],alpha=0.7)
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Alfa helisa'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='Mostek beta'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='yellow', markersize=10, label='Nić'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='3-10 helisa'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='violet', markersize=10, label='Pi helisa'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='pink', markersize=10, label='Turn'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='Bend'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='None')
]

plt.legend(handles=legend_elements)
plt.title('Wykres Ramachandrana')
plt.xlabel('Wartość Phi')
plt.ylabel('Wartość Psi')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks([-180,-120,-60,0,60,120, 180]) 
plt.yticks([-180,-120,-60,0,60,120, 180]) 
plt.savefig('wykres.png')
