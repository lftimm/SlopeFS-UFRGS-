from typing import Optional, Dict, Tuple, List
import math
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

class SoilSpace:
    """
        SoilSpace Class:
            It represents the entire space of the problem.
            - The physical characteristics of the slope
            - The mechnanical characteristics of the soil
        Uses:
            Used standalone to define the soil and slope before sending it to Model and Soil_fs
            Default values present to simplify analysing behaviour
    """

    def __init__(self, c=20, phi=30, gam=18.5, alp=45, h=15, num_slice=50, circle: Dict[str, float] = None):
        self.properties = {
            'c': c,
            'phi': math.radians(phi),
            'gam': gam,
            'alp': math.radians(alp),
            'h': h,
            'num_slice': num_slice
        }

        self.update_slope_len()

        if not circle:
            self.update_circle()
        else:
            assert ['xc', 'yc', 'R'] == list(circle.keys()), f'Wrong dictionary keys: {list(circle.keys())}'
            self.update_circle(xc=circle['xc'], yc=circle['yc'], r=circle['R'])

    def update_slope_len(self):
        self.properties['slope_len'] = self.properties['h'] / math.tan(self.properties['alp'])

    def update_circle(self, xc=None, yc=None, r=None):
        if xc is None and yc is None and r is None:
            xc = 0.5 * self.properties['h'] / math.tan(self.properties['alp'])
            yc = (4 / 3) * self.properties['h']
            r = math.sqrt(xc ** 2 + yc ** 2)

        circle = {
            'xc': xc,
            'yc': yc,
            'R': r
        }

        self.properties['Circle'] = circle

    def __str__(self):
        return f'{self.properties}'


class Model:
    """
    Model class:
        Responsible for finding the values associated with the critical circle.
    Uses:
        Used for finding the variables needed to calculate the FS of the slope.
    """

    def __init__(self, soil: Optional[SoilSpace] = None):
        if soil is None:
            soil = SoilSpace()

        self.soil = soil
        self.sl: Dict[str, float] = soil.properties

        self.circle = self.sl['Circle']

        self.points: Tuple[Tuple, Tuple] = self.intersec()
        self.c_points: List[Tuple] = self.split_geometry()
        self.polys: List[np.array] = self.mk_polys()

        self.dxs, self.alphas = self.calc_alphas()
        self.polys_A: List[float] = self.calc_areas()

    def intersec(self) -> Tuple[Tuple, Tuple]:
        """
            Calculates the points of intersection of the slope and the circle.
            It uses second degree equations to find the intersections.
            It returns the points in order from left to right.
        """
        t = math.tan(self.sl['alp'])
        h = self.sl['h']
        l = self.sl['slope_len']
        r, xc, yc = self.circle['R'], self.circle['xc'], self.circle['yc']

        a = 1 + t ** 2
        b1 = -2 * xc
        b = b1 - 2 * t * yc
        c = xc ** 2 + yc ** 2 - r ** 2

        delta = b ** 2 - 4 * a * c


        assert delta > 0, 'Math error, delta <= 0, Circle doesn\'t intersect slope.'

        def f(x):
            if 0 <= x < l:
                return x, t * x
            elif x < 0:
                return xc - math.sqrt(r ** 2 - yc ** 2), 0
            else:
                return xc + math.sqrt(r ** 2 + 2 * h * yc - h ** 2 - yc ** 2), h

        x1 = (-b + math.sqrt(delta)) / (2 * a)
        x2 = (-b - math.sqrt(delta)) / (2 * a)
        p1 = f(x1)
        p2 = f(x2)

        p_l = p1 if p1[0] == min(p1[0], p2[0]) else p2
        p_r = p2 if p_l[0] == p1[0] else p1

        return p_l, p_r

    def split_geometry(self) -> List[Tuple]:
        """
            It splits the circle into equal parts based on the number of slices given.
            Together there is the total_angle method, it measures the total angle of the intersection points.
            It returns a list of tuples containing the points.
            One thing might be removed, in the definition of f() there is a rounding done with the map() function.
        """
        # a = math.tan(self.sl['alp'])
        p_l, p_r = self.points
        r, xc, yc = self.circle['R'], self.circle['xc'], self.circle['yc']
        ns = self.sl['num_slice']

        v_p_r = np.array(p_r)
        v_c = np.array([xc, yc])

        def total_angle(p1, p2, xc, yc):
            dist = lambda p1, p2: math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
            ct = (xc, yc)
            a, b, c = dist(p1, p2), dist(p1, ct), dist(p2, ct)
            tot_angle = math.acos(-(a ** 2 - (b ** 2 + c ** 2)) / (2 * b * c))

            return tot_angle

        tot_a = total_angle(p_l, p_r, xc, yc)
        alp = tot_a / ns
        gam = math.atan(abs((v_c[1] - v_p_r[1]) / (v_c[0] - v_p_r[0])))

        f = lambda n: map(lambda x: x,
                          r * np.array([math.cos(-(gam + n * alp)), math.sin(-(gam + n * alp))]) + v_c)
        return [tuple(f(n)) for n in range(ns + 1)]

    def mk_polys(self) -> np.array:
        """
            This method creates the polygons whose areas are going to be calculated.
            It takes the list of points in the circle, reflects them into the corresponding part of the surface.
            Together with it there is the pair_points method, which takes the previous points and orders them counter-clockwise.
            It returns an array of 4x1 arrays with the points .
        """
        c_parts = self.c_points
        pts_x, pts_y = zip(*c_parts)
        a = math.tan(self.sl['alp'])
        h = self.sl['h']
        l = self.sl['slope_len']

        def f(x):
            if 0 <= x <= l:
                return a * x
            elif x > l:
                return h
            else:
                return 0

        up_c_parts = [(x, f(x)) for x in pts_x]

        full_points = [c_parts, up_c_parts]

        def pair_points(points):
            polys = []
            for i in range(0, len(points[0])):
                try:
                    polys.append([points[0][i], points[1][i], points[1][i + 1], points[0][i + 1]])
                except IndexError:
                    pass

            return np.array(polys)

        a = pair_points(full_points)

        return a

    def calc_alphas(self):
        """
            Alpha angle of each slice of the circle (each polygon).
            Utilizes the arctan(x).
        """

        polys = self.polys

        alp = lambda dy, dx: math.atan(dy / dx)
        n_polys = polys.shape[0]

        # Potential performance lost
        alphas = []
        dxs = []
        for i in range(n_polys):
            dy = polys[i][0][1] - polys[i][-1][1]
            dx = polys[i][0][0] - polys[i][-1][0]

            dxs.append(dx)
            alpha = alp(dy,dx)
            alphas.append(alpha)

        return dxs, alphas

    def calc_areas(self) -> List[float]:
        """
            It calculates the areas of the polygons.
            It uses the shoelace formula for calculating the area.
            It returns an array containing the areas of each of the polygons.
        """
        p = self.polys

        areas = []

        for poly in p:
            n = len(poly)
            a = 0
            for i in range(n):
                x1, y1 = poly[i]
                x2, y2 = poly[(i + 1) % n]
                a += (x1 * y2 - x2 * y1) / 2
            areas.append(a)

        return areas


class SoilFs:
    """
        SoilFs:
            Wraps all the code together, calling the model, calculating FS and calling the view.
        Uses:
            To be used by the end user when dealing with the program.
    """

    def __init__(self, soil: Optional[SoilSpace] = None, methods=['Fellenius', 'OSM'], minimize=True):
        self.model = Model(soil)
        self.sl = self.model.soil.properties
        self.c0 = self.sl['Circle']

        self.results = self.end_results(methods,minimize)
        self.view = lambda: View.plot(self)

    def end_results(self, methods, minimize) -> Dict[str, float]:
        """
            Finalizes everything by gathering all the previous steps and calculating the FS.
            Bishop is an implicit equation, its roots are found using Newton's method (via Scipy.optimize)
            It returns both values in a dictionary format, according to the format used by the SoilSpace class.
        """
        gam = self.sl['gam']
        c = self.sl['c']
        phi = self.sl['phi']

        are = self.model.polys_A
        alp = self.model.alphas
        dxs = self.model.dxs

        results = {'c0': list(self.sl['Circle'].values())}
        if 'Fellenius' in methods:
            fel = self.fel0(c, gam, are, alp, phi, dxs)
            results['fel0'] = fel

            if minimize:
                min_fel = self.min_fel(self.model.soil)
                results['min_fel_c0'] = list(min_fel.x)
                results['min_fel'] = min_fel.fun
                self.sl['Circle'] = self.c0

        if 'OSM' in methods:
            osm0 = self.bis0(c, gam, are, alp, phi, dxs)
            results['osm0'] = osm0
            if minimize:
                min_osm = self.min_bis(self.model.soil)
                results['min_osm_c0'] = list(min_osm.x)
                results['min_osm'] = min_osm.fun
                self.sl['Circle'] = self.c0

        return results

    @staticmethod
    def fel0(c, gam, are, alp, phi, dxs):
        size = len(are)
        fell1 = sum(
            [c * dxs[i] / math.cos(alp[i]) for i in range(size)]
        )
        fell2 = sum(
            [math.tan(phi) * gam * are[i] * math.cos(alp[i]) for i in range(size)]
        )
        fell3 = sum(
            [gam * are[i] * math.sin(alp[i]) for i in range(size)]
        )

        fs = (fell1 + fell2) / fell3
        return fs

    @staticmethod
    def bis0(c, gam, are, alp, phi, dxs):
        def bishop_calc(fs):
            size = len(are)
            bip1 = (sum(
                [gam * are[i] * math.sin(alp[i]) for i in range(size)])) ** -1
            bip2 = sum(
                [(c * dxs[i] + gam * are[i] * math.tan(phi)) /
                 (math.cos(alp[i]) + math.sin(alp[i]) * math.tan(phi) / fs) for i in range(size)])
            return fs - bip1 * bip2

        return optimize.newton(bishop_calc, x0=2)

    @staticmethod
    def min_fel(slope):
        def f(c):
            slope.update_circle(c[0], c[1], c[2])
            try:
                model = Model(slope)
                c = model.sl['c']
                gam = model.sl['gam']
                are = model.polys_A
                alp = model.alphas
                phi = model.sl['phi']
                dxs = model.dxs
                fs = SoilFs.fel0(c, gam, are, alp, phi, dxs)
                if np.isnan(fs):
                    return np.inf
            except AssertionError:
                return 9999
            return fs

        c0 = list(slope.properties['Circle'].values())
        fun = optimize.minimize(f, c0, method='SLSQP')
        return fun

    @staticmethod
    def min_bis(slope):
        def f(c):
            slope.update_circle(c[0], c[1], c[2])
            try:
                model = Model(slope)
                if any(dx == 0 for dx in model.dxs):
                    return np.inf
                c = model.sl['c']
                gam = model.sl['gam']
                are = model.polys_A
                alp = model.alphas
                phi = model.sl['phi']
                dxs = model.dxs
                fs = SoilFs.bis0(c, gam, are, alp, phi, dxs)
                if np.isnan(fs):
                    return np.inf
            except AssertionError:
                return 9999
            return fs

        c0 = list(slope.properties['Circle'].values())
        fun = optimize.minimize(f, c0, method='SLSQP')
        return fun

    def __str__(self):
        return '\n'.join(f'{key}: {value}' for key, value in self.results.items())


class View:
    """
    View Class:
        It presents all the information to the user.
    Uses:
        To be called by SoilFS.
    """
    @staticmethod
    def plot(full: SoilFs):
        pprt = full.sl
        pprt2 = full.model
        xcord,ycord = zip(*pprt2.c_points[::-1])
        h = pprt['h']
        len = pprt['slope_len']
        xc, yc, R = list(pprt['Circle'].values())

        plt.xlim([-len, 2*len])
        plt.ylim([-len, 2*len])
        plt.plot([-2*len, 0], [0, 0], color='#000000')
        plt.plot([0, len], [0, h],  color='#000000')
        plt.plot([len, 2*len], [h, h], color='#000000')

        plt.plot(xcord, ycord, color='#3200ff')
        plt.plot([xc, xcord[-1]], [yc, ycord[-1]],color='#3200ff')
        plt.plot([xcord[0], xc], [ycord[0], yc],color='#3200ff')

        plt.plot(xc,yc,color='#3200ff')

        title = ''
        to_plot = full.results.keys()

        if 'fel0' in to_plot:
            fel = '{:.2f}'.format(full.results['fel0'])
            title += f'Fellenius={fel} '
        if 'osm0' in to_plot:
            osm = '{:.2f}'.format(full.results['osm0'])
            title += f'OSM={osm}'

        xc,yc,R = map(lambda x: round(x,2), [xc,yc,R])
        title += f'\n c=(xc:{xc},yc:{yc},R:{R})'

        plt.title(title)

        polys = pprt2.polys[::-1]
        for poly in polys:
            a, b = zip(*poly)
            plt.plot(a, b, color='#ff3200')

        plt.show()
        plt.close()

        to_look = ['min_fel', 'min_osm']
        found = [_ in to_plot for _ in to_look]
        size_of_plot = sum(found)
        i = 1
        if any(found):
            plt.figure(figsize=(13,5))
            if found[0]:
                plt.subplot(1,size_of_plot,i)
                plt.xlim([-len, 2 * len])
                plt.ylim([-len, 2 * len])
                plt.plot([-2 * len, 0], [0, 0], color='#000000')
                plt.plot([0, len], [0, h], color='#000000')
                plt.plot([len, 2 * len], [h, h], color='#000000')

                fel_c0 = full.results['min_fel_c0']

                nc = {'xc':fel_c0[0],'yc':fel_c0[1],'R':fel_c0[2]}
                soiln = pprt.copy()
                del soiln['Circle']
                del soiln['slope_len']
                soiln['circle'] = nc
                soiln['alp'] = math.degrees(soiln['alp'])

                s1 = SoilSpace(**soiln)
                s2 = Model(s1)
                xcord,ycord = zip(*s2.c_points[::-1])
                xc, yc, R = list(s1.properties['Circle'].values())

                plt.plot([xcord[0],xc],[ycord[0],yc],color='#3200ff')
                plt.plot([xc, xcord[-1]], [yc, ycord[-1]], color='#3200ff')
                plt.plot(xcord, ycord, color='#3200ff')

                polys = s2.polys[::-1]
                for poly in polys:
                    a, b = zip(*poly)
                    plt.plot(a, b, color='#ff3200')

                ff = '{:.2f}'.format(full.results['min_fel'])

                xcn,ycn,Rn= map(lambda x:round(x,2),[xc,yc,R])

                plt.title('FS_{Fellenius}=' + f'{ff}' + '\n' + f'c=(xc:{xcn},yc:{ycn},R:{Rn})')
                i += 1

            if found[1]:
                plt.subplot(1, size_of_plot, i)
                plt.xlim([-len, 2 * len])
                plt.ylim([-len, 2 * len])
                plt.plot([-2 * len, 0], [0, 0], color='#000000')
                plt.plot([0, len], [0, h], color='#000000')
                plt.plot([len, 2 * len], [h, h], color='#000000')

                osm_c0 = full.results['min_osm_c0']
                nc = {'xc':osm_c0[0],'yc':osm_c0[1],'R':osm_c0[2]}
                soiln = pprt.copy()
                del soiln['Circle']
                del soiln['slope_len']
                soiln['circle'] = nc
                soiln['alp'] = math.degrees(soiln['alp'])
                s1 = SoilSpace(**soiln)
                s2 = Model(s1)
                xcord,ycord = zip(*s2.c_points[::-1])
                xc, yc, R = list(s1.properties['Circle'].values())

                plt.plot(xcord, ycord, color='#3200ff')
                plt.plot([xc, xcord[-1]], [yc, ycord[-1]], color='#3200ff')
                plt.plot([xcord[0], xc], [ycord[0], yc], color='#3200ff')

                polys = s2.polys[::-1]
                for poly in polys:
                    a, b = zip(*poly)
                    plt.plot(a, b, color='#ff3200')

                osm = '{:.2f}'.format(full.results['min_osm'])

                xcn,ycn,Rn= map(lambda x:round(x,2),[xc,yc,R])
                plt.title('FS_{OSM}=' + f'{osm}' + '\n' + f'c=(xc:{xcn},yc:{ycn},R:{Rn})')

                i += 1

            plt.show()


