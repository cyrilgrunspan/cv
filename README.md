# cv

CV de Cyril Grunspan (`index.html`, `cv.pdf`).

`cagnotte.html` : moteur de calcul en ligne de la cagnotte critique d'une course avec investissement perdu
(C. Grunspan et R. Pérez-Marco, *Is the game worth the candle? Optimal stopping in a race with sunk investment*,
définition 3.17 et algorithme 1). Entrées : l'état (a, h), la part q du challenger et la longueur z ;
sortie : un encadrement de la cagnotte critique v*_z(a, h). Même algorithme en Python : `cagnotte.py`
(`python3 cagnotte.py a h q z`). Calcul en virgule flottante, non certifié.
