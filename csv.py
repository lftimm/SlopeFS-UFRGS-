'''
def variate_soil(self, var: str, end: float, step: float,
                 start=None, base: Optional[SoilSpace] = None) -> list[SoilSpace]:
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

def change_to_csv(self, var: str, end: float, step: float, start=None, filename='default.txt') -> None:
    """
        Method for analysis.
        One variable can be isolated for study, given a range from start to finish.
        Its output is written in a csv file that can be imported by Excel.
    """
    if os.path.exists(filename):
        counter = 2
        name, extension = os.path.splitext(filename)
        while True:
            new_name = f'{name}({counter}){extension}'
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
'''