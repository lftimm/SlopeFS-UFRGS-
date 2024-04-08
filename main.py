import SlopeFs as sfs
import scipy as sp
import numpy as np

def minimize_fel(soil):
    def fel(circ):
        dcirc = {
            'xc': circ[0],
            'yc': circ[1],
            'R': circ[2]
        }
        soil = sfs.SoilSpace(circle=dcirc)

        props = soil.properties
        model = sfs.Model(soil)

        c = props['c']
        gam = props['gam']
        are = model.polys_A
        alp = model.alphas
        phi = props['phi']
        slen = props['slope_len']
        n = props['num_slice']

        fellen = sfs.SoilFs.fellenius(c, gam, are, alp, phi, slen, n)
        return fellen

    return sp.optimize.minimize(fel, list(soil.properties['Circle'].values()))

def minimize_bis(soil):
    def bis(circ):
        dcirc = {
            'xc': circ[0],
            'yc': circ[1],
            'R': circ[2]
        }
        soil = sfs.SoilSpace(circle=dcirc)

        props = soil.properties
        model = sfs.Model(soil)

        c = props['c']
        gam = props['gam']
        are = model.polys_A
        alp = model.alphas
        phi = props['phi']
        dxs = model.dxs

        bishop = sfs.SoilFs.bishop(c, gam, are, alp, phi, dxs)

        return bishop

    return sp.optimize.minimize(bis, list(soil.properties['Circle'].values()))


def main():

if __name__ == '__main__':
    main()
