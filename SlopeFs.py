from typing import Optional
import os

import math
import numpy as np
from scipy import optimize

from copy import deepcopy
import csv


class SoilSpace:
    """
        SoilSpace Class:
            It represents the entire space of the problem.
            - The physical characteristics of the slope
            - The mechnanical characteristics of the soil
        Uses:
            Used standalone to define the soil before sending it to Model and Soil_fs
            Default values present to simplify analysing behaviour
    """

    def __init__(self, c=20, phi=25, gam=18.5, alp=30, h=15, num_slice=5):
        self.properties = {
            'c': c,
            'phi': math.radians(phi),
            'gam': gam,
            'alp': math.radians(alp),
            'h': h,
            'num_slice': num_slice
        }

        self.update_slope_len()
        self.update_circle()

    def update_slope_len(self):
        self.properties['slope_len'] = self.properties['h'] / math.tan(self.properties['alp'])

    def update_circle(self):
        circle = {
            'xc': 0.5 * self.properties['h'] / math.tan(self.properties['alp']),
            'yc': 1.67 * self.properties['h']
        }

        circle['R'] = 1.5 * math.sqrt(circle['xc'] ** 2 + circle['yc'] ** 2)

        self.properties['Circle'] = circle

    def __str__(self):
        return f'{self.properties}'


class Model:
    """
    Model class:
    Responsible for calculating the factor of safety of the slope.
    """

    def __init__(self, soil: Optional[SoilSpace] = None):
        if soil is None:
            soil = SoilSpace()

        self.sl = soil.properties

        self.circle = self.sl['Circle']

        self.points = self.intersec()
        self.c_points = self.splitgeometry()
        self.polys = self.mk_polys()

        self.dxs, self.alphas = self.calc_alphas()
        self.polys_A = self.calc_areas()

        self.results = self.end_results()

    def intersec(self):
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

        if delta > 0:
            def f(x):
                if 0 < x < l:
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

    def splitgeometry(self):
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

        tot_a = self.total_angle(p_l, p_r, xc, yc)
        alp = tot_a / ns
        gam = math.atan(abs((v_c[1] - v_p_r[1]) / (v_c[0] - v_p_r[0])))

        f = lambda n: map(lambda x: round(x, 2),
                          r * np.array([math.cos(-(gam + n * alp)), math.sin(-(gam + n * alp))]) + v_c)
        return [tuple(f(n)) for n in range(ns + 1)]

    def total_angle(self, p1, p2, xc, yc):
        dist = lambda p1, p2: math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
        ct = (xc, yc)
        a, b, c = dist(p1, p2), dist(p1, ct), dist(p2, ct)
        tot_angle = math.acos(-(a ** 2 - (b ** 2 + c ** 2)) / (2 * b * c))

        return tot_angle

    def mk_polys(self):
        """
            This method creates the polygons whose areas are going to be calculated.
            It takes the list of points in the circle, reflects them into the corresponding part of the surface.
            Together with it there is the pair_points method, which takes the previous points and orders them counter-clockwise.
            It returns an array of 4x1 arrays with the points .
        """
        c_parts = self.c_points
        pts_x, pts_y = zip(*c_parts)
        a = self.sl['alp']
        h = self.sl['h']
        l = self.sl['slope_len']

        def f(x):
            if 0 <= x <= l:
                return round(a * x, 2)
            elif x > l:
                return h
            else:
                return 0

        up_c_parts = [(x, f(x)) for x in pts_x]

        full_points = [c_parts, up_c_parts]

        return self.pair_points(full_points)

    def pair_points(self, points):
        polys = []
        for i in range(0, len(points[0])):
            try:
                polys.append([points[0][i], points[1][i], points[1][i + 1], points[0][i + 1]])
            except IndexError:
                pass

        return np.array(polys)

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
            alphas.append(alp(dy, dx))

        return dxs, alphas

    def calc_areas(self):
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

    def end_results(self):
        """
            Finalizes everything by gathering all the previous steps and calculating the FS.
            Bishop is an implicit equation, its roots are found using Newton's method (via Scipy.optimize)
            It returns both values in a dictionary format, according to the format used by the SoilSpace class.
        """
        gam = self.sl['gam']
        c = self.sl['c']
        phi = self.sl['phi']
        slen = self.sl['slope_len']
        n = self.sl['num_slice']

        are = self.polys_A
        alp = self.alphas
        dxs = self.dxs

        u = 0
        size = len(are)

        def fellenius():
            fell1 = sum(
                [c * slen / n + (gam * are[i] * math.cos(alp[i]) - u * slen) * math.tan(phi) for i in range(size)])
            fell2 = sum([gam * are[i] * math.sin(alp[i]) for i in range(size)])

            fs = fell1 / fell2
            return fs

        def bishop(fs):
            bip1 = (sum([gam * are[i] * math.sin(alp[i]) for i in range(size)])) ** -1
            bip2 = sum([(c * dxs[i] + gam * are[i] * math.tan(phi)) / (
                    math.cos(alp[i]) + math.sin(alp[i]) * math.tan(phi) / fs) for i in range(size)])
            return fs - bip1 * bip2

        bip = optimize.newton(bishop, x0=2)
        fel = fellenius()

        return {'fel': fel, 'bis': bip}


class SoilFs:
    """
        SoilFs:
        Serves as a wrapper for everything.
        Houses methods for showing and analyzing data.
    """

    def __init__(self, soil: Optional[SoilSpace] = None):
        self.soil = soil
        self.model = Model(self.soil)

    def variate_soil(self, var: str, end: float, step: float, start=None, base: Optional[SoilSpace] = None) -> list[SoilSpace]:
        if base is None:
            base = self.soil

        assert var in base.properties, f'Invalid key: Property \'{var}\' doesn\'t exist.'
        assert type(base.properties[var]) is int or float, f'Invalid key \'{var}\': must be numeric.'
        assert step > 0, f'Invalid Range: step must be greater than 0.'

        def asserts(start, end, step):
            assert end > start, f'Invalid Range: end > start.'
            n_steps = (end - start) / step
            assert n_steps > 1, f'Invalid Range: Step too big. {n_steps}'

        copy = deepcopy(base)
        samples = []

        if var in ['alp', 'phi']:
            if start is None:
                start = round(math.degrees(copy.properties[var]))

            asserts(start, end, step)

            while start <= end:
                copy.properties[var] = math.tan(math.radians(start))
                copy.update_slope_len()
                copy.update_circle()
                samples.append(deepcopy(copy))
                start += step
        else:
            if start is None:
                start = copy.properties[var]

            asserts(start, end, step)

            while start <= end:
                copy.properties[var] = start
                copy.update_slope_len()
                copy.update_circle()
                samples.append(deepcopy(copy))
                start += step

        return samples

    def change_to_csv(self, var: str, end: float, step: float, start=None, filename='default.txt'):
        """
            Method for analysis.
            One variable can be isolated for study, given a range from start to finish.
            Its output is written in a csv file that can be imported by Excel.
        """
        if os.path.exists(filename):
            counter = 2
            name, extension = os.path.splitext(filename)
            while True:
                new_name = f'{name}({counter}).{extension}'
                if not os.path.exists(new_name):
                    filename = new_name
                    break
                counter += 1

        fields = ['c', 'phi', 'gam', 'alp', 'h', 'num_slice', 'slope_len', 'fel', 'bis']

        copy = deepcopy(self.soil)
        variations = self.variate_soil(var, end, step, start, base=copy)
        normal_rows = [soil.properties for soil in variations]
        result_rows = [Model(soil).results for soil in variations]
        final_rows = []
        for row1, row2 in zip(normal_rows, result_rows):
            updated_row = {**row1, **row2}
            del updated_row['Circle']
            final_rows.append(list(updated_row.values()))

        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(fields)
            writer.writerows(final_rows)
