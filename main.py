import SlopeFs as sfs

def main():
    soil = sfs.SoilSpace(h=12)
    fs = sfs.SoilFs(soil)
    print(fs)

if __name__ == '__main__':
    main()
