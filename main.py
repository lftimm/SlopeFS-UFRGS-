import SlopeFs


def main():
    soil = SlopeFs.SoilSpace(h=10, num_slice=100)
    model = SlopeFs.SoilFs(soil)
    model.change_to_csv('h', 100, 10, filename='test_h.txt')


if __name__ == '__main__':
    main()
