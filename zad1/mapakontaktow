#!pip install biopython
#!pip install matplotlib

import Bio.PDB as PDB
import matplotlib.pyplot as plt
import argparse

#nie wiemy jak uruchomić plik skryptowy z argumentami w google colab więc nalezy
#zmienic wartosc zmiennej kod

#parserI = argparse.ArgumentParser()
#parserI.add_argument("--pdb", help="Nazwa pliku PDB")
#args = parserI.parse_args()
#kod = args.pdb

kod = "2hhb"

pdbl = PDB.PDBList()
plik = pdbl.retrieve_pdb_file(kod, pdir=".", file_format="pdb")

parser = PDB.PDBParser()
struktura = parser.get_structure("RNA", plik)

# c alfa
atomy = []
for model in struktura:
    for lancuch in model:
        for reszta in lancuch:
            if reszta.has_id("CA"):
                atomy.append(reszta["CA"])

# obliczanie dystansow miedzy atomami c alfa
dystansy = []
for i in range(len(atomy)):
    for j in range(i+1, len(atomy)):
        dystans = atomy[i] - atomy[j]
        if dystans <= 8:
            dystansy.append((i, j))

# mapa kontaktow
mapa_kontaktow = [[0] * len(atomy) for _ in range(len(atomy))]
for i, j in dystansy:
    mapa_kontaktow[i][j] = 1
    mapa_kontaktow[j][i] = 1

# Wizualizacja mapy kontaktow
plt.imshow(mapa_kontaktow, cmap="binary")
plt.xlabel("Index reszty")
plt.ylabel("Index reszty")
plt.title("Mapa kontaktow")
plt.show()
