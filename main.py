from SlopeFs import SoilSpace, SoilFs


def main():
    soil = SoilSpace(h=5, num_slice=500)
    model = SoilFs(soil)
    model.change_to_csv('h', 100, 1, filename='test_h.txt')


if __name__ == '__main__':
    main()